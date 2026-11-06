function disp_struct = get_display(P,xname,out_names)
% DISP_STRUCT = GET_DISPLAY(P,XNAME,OUT_NAMES) constructs the terms of the
% expression
%   int ... int K(s,s') x(s') ds' ... ds' + ...
% of the PI operator P acting on a function x, for use in the 'display'
% method of 'sopvar', and in the 'get_display' methods of the classes of
% operators that have 'sopvar' objects as factors or blocks.
%
% INPUTS
% - P:          'sopvar' object representing a PI operator;
% - xname:      (optional) char, the name of the function that the operator
%               acts on. Defaults to x, unless a variable is also named x
%               (ignoring a number at the end of the name), in which case
%               it is the first of u, w, v that is not the name of a
%               variable. If the operator acts on R, the name '1' indicates
%               the constant 1, which is not displayed;
% - out_names:  (optional) 1 x M cellstr, the names to display for the
%               output variables P.vars.out. Defaults to the names given by
%               'varname2display';
%
% OUTPUTS
% - disp_struct:      struct with the fields
%       terms:  1 x n cell of char, the terms of the sum, without their
%               sign, see 'display_term';
%       signs:  1 x n array of 1 and -1, the sign of each term;
%       use_K:  logical, true if a kernel is displayed by name. The kernel
%               K_i can then be extracted by calling P.kernels{i};
%       ok:     logical, true if the operator is displayed explicitly (no
%               kernel by name), with at most 6 terms of at most 100
%               characters, so that it may be part of a larger display;
%       xname:  char, the name of the function that the operator acts on;
%       is_id:  logical, true if P is an identity operator, in which case
%               the only term is the function with its arguments, x(s);
%       rows:   1 x M struct array with fields terms and signs, as above,
%               for each element of the output of a matrix-valued operator
%               P, e.g. the ith element of P x. Empty unless P is
%               matrix-valued, and each element is displayed explicitly on
%               one line. If P maps from a vector-valued space, the terms
%               act on the elements x_1, x_2, ... of the function x;
%
% NOTES
% - The terms correspond to the parameters P.params, see 'display_term'.
%   They are numbered like the elements of P.params, so a term with a zero
%   kernel is not displayed, but does take an index.
% - A kernel with at most 4 nonzero coefficients is displayed explicitly,
%   if P is scalar-valued (P.dims=[1,1]). The kernels of a matrix-valued
%   operator are displayed by name in terms, and explicitly in rows if the
%   element of every row has at most 4 nonzero coefficients in each kernel,
%   and the row has at most 100 characters.
% - An identity operator is displayed as x(s), see 'is_id_op'.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PIETOOLS - get_display
%
% Copyright (C) 2026 PIETOOLS Team
%
% This program is free software; you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation; either version 2 of the License, or
% (at your option) any later version.
%
% This program is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
%
% You should have received a copy of the GNU General Public License
% along with this program; if not, write to the Free Software
% Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% If you modify this code, document all changes carefully and include date
% authorship, and a brief description of modifications
%
% DJ, 11/06/2026: Initial coding;

% Default settings
max_terms = 6;          % maximal number of terms of an operator displayed explicitly
max_chars = 100;        % maximal number of characters of each such term
max_row_chars = 100;    % maximal number of characters of each row of a matrix-valued
                        % operator displayed explicitly

if nargin<2 || isempty(xname)
    % x, unless a variable is also named x (ignoring a number at the end of
    % the name), in which case the first of u, w, v that is not
    var_bases = regexprep([reshape(P.vars.out,1,[]),reshape(P.vars.in,1,[])],'_?\d+$','');
    for xname = 'xuwv'
        if ~ismember(xname,var_bases)
            break
        end
    end
end
if nargin<3
    out_names = [];
end
disp_struct = struct('terms',{cell(1,0)},'signs',zeros(1,0),'use_K',false,'ok',false,'xname',xname,'is_id',false,'rows',struct('terms',{},'signs',{}));

if is_id_op(P)
    if isempty(out_names)
        out_names = cellfun(@varname2display,reshape(P.vars.out,1,[]),'UniformOutput',false);
    end
    disp_struct.terms = {xname};
    if ~isempty(out_names)
        disp_struct.terms = {[xname,'(',strjoin(reshape(out_names,1,[]),','),')']};
    end
    disp_struct.signs = 1;        disp_struct.ok = true;        disp_struct.is_id = true;
    return
end

for k=1:numel(P.params)
    if nnz(P.params{k})==0
        continue
    end
    [term_str,term_sign,use_K] = display_term(P,k,xname,out_names);
    disp_struct.terms{end+1} = term_str;      disp_struct.signs(end+1) = term_sign;
    disp_struct.use_K = disp_struct.use_K || use_K;
end
disp_struct.ok = ~disp_struct.use_K && numel(disp_struct.terms)<=max_terms && all(cellfun(@numel,disp_struct.terms)<=max_chars);

if P.dims(1)*P.dims(2)>1
    disp_struct.rows = display_rows(P,xname,out_names,max_row_chars);
end

end

function rows = display_rows(P,xname,out_names,max_chars)
% ROWS = DISPLAY_ROWS(P,XNAME,OUT_NAMES,MAX_CHARS) constructs the terms of
% each element of the output of the matrix-valued operator P acting on x,
% for 'get_display'.
%
% INPUTS
% - P:          'sopvar' object, with P.dims(1)*P.dims(2)>1;
% - xname, out_names: as for 'get_display';
% - max_chars:  maximal number of characters of the terms of a row;
%
% OUTPUTS
% - rows:       1 x m struct array with fields terms and signs, the sum
%               that defines element i of P x, with one term for each
%               element K_ij of the kernel of each parameter, acting on x_j
%               (or x if P.dims(2)=1). Empty if a kernel element has more
%               than 4 nonzero coefficients, or if the terms of a row take
%               more than max_chars characters;

m = P.dims(1);      n = P.dims(2);
rows = struct('terms',{},'signs',{});

x_names = repmat({xname},1,n);
if n>1
    for j=1:n
        % A comma separates the index from a subscript at the end of the name
        sep = repmat(',',1,xname(end)>=8320 && xname(end)<=8329);
        x_names{j} = [xname,sep,sprintf('%d',j)-'0'+8320];
    end
end

row_terms = repmat({cell(1,0)},1,m);       row_signs = repmat({zeros(1,0)},1,m);
for k=1:numel(P.params)
    if nnz(P.params{k})==0
        continue
    end
    [i,j] = decode_params(P,k);
    entry_idcs = unique([i,j],'rows');
    for t=1:size(entry_idcs,1)
        [term_str,term_sign,use_K] = display_term(P,k,x_names{entry_idcs(t,2)},out_names,entry_idcs(t,:));
        if use_K
            return
        end
        row_terms{entry_idcs(t,1)}{end+1} = term_str;      row_signs{entry_idcs(t,1)}(end+1) = term_sign;
    end
end
for i=1:m
    if isempty(row_terms{i})
        row_terms{i} = {'0'};      row_signs{i} = 1;
    end
    if sum(cellfun(@numel,row_terms{i}))+3*(numel(row_terms{i})-1)>max_chars
        return
    end
end
rows = struct('terms',row_terms,'signs',row_signs);

end
