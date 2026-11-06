function display(Rop,name,width_max,expand,call_name)
% DISPLAY(ROP,NAME,WIDTH_MAX,EXPAND,CALL_NAME) displays the tensor-PI
% operator ROP in the command window as
%   (NAME x)(s_1,...,s_d) = (R_1 x_1)(s_1) o ... o (R_d x_d)(s_d),
% for o the tensor product, with the operators R_j displayed as the actual
% operators if they are small enough, or as a sum of such products if ROP
% has multiple terms. Called automatically when displaying a 'tensopvar'
% object without semicolon.
%
% INPUTS
% - Rop:    'tensopvar' object;
% - name:   (optional) char, the name to display for the operator. Defaults
%           to the name of the variable passed to display, or to 'Rop' if
%           there is none. The name 'ans' is replaced by 'Rop';
% - width_max:  (optional) maximal number of characters of a line with
%           multiple terms. Defaults to the width of the command window, or
%           to 100 if this cannot be determined;
% - expand: (optional) logical, false to always display the operators R_ij
%           as such, instead of as the actual operators if they are small
%           enough. Defaults to true;
% - call_name:  (optional) char, the expression with which the user can
%           call the operator, used in the call to extract an operator
%           R_ij. Defaults to name, also if that is ans, and is set by the
%           display of a container that has Rop as its only block, e.g.
%           'Rop.ops{1,1}';
%
% OUTPUTS
% Display of the operator in the command window, with the terms constructed
% by 'get_display'.
%
% NOTES
% - The product of the factors is displayed as the circled times symbol for
%   a tensor product (Rop.otype(1)=2), and as .* for an elementwise product
%   (Rop.otype(1)=1). The same holds for the inputs x_1,...,x_d, with
%   Rop.otype(2).
% - The output variables are those in Rop.vars.out, with a subscript to
%   distinguish different copies of the same variable, numbered over the
%   factors in order, as in 'apply_op'. A number at the end of the name of
%   a variable is displayed as a subscript, followed by the subscript of
%   the copy, separated by a comma, e.g. s1 as s_1 and its second copy as
%   s_1,2, see 'varname2display'.
% - The operators R_ij, with i the term and j the factor, are displayed as
%   the actual operators, acting on x_j, if all of them in the terms that
%   are displayed are small enough, see 'get_display'. Otherwise, they are
%   displayed as R_ij, dropping i if there is only one term, and j if
%   there is only one factor, and can be extracted by calling Rop.ops{i,j}.
% - Up to three terms are displayed. For more terms, the first and last
%   term are displayed, and the other terms are replaced by an ellipsis.
%   The terms are placed on the same line, unless the line would have more
%   characters than fit in the command window, in which case the whole next
%   term is moved to the next line, so that no term is split by the command
%   window, with the + aligned under the equality sign.
% - An identity operator is displayed as (NAME x)(s) = x(s), with the
%   output variables.
% - The input is named x, unless this is the name of a variable (ignoring a
%   number at the end of the name), in which case it is named u, w, or v,
%   whichever is the first that is not the name of a variable.
% - If there are multiple factors, the input space is indicated by the line
%   "where x(s_1,...,s_d) = x_1(s_1) o ... o x_d(s_d)", with the input
%   variables in Rop.vars.in, and their copies numbered as the output
%   variables. With more than four factors, the factors x_2,...,x_{d-1} are
%   replaced by an ellipsis. For a pointwise product of the inputs
%   (Rop.otype(2)=1), the factors are set equal, x_1(s) = x_2(s) = ...
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PIETOOLS - display
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
max_factors = 4;        % maximal number of factors of x displayed; with more
                        % factors, the middle factors are replaced by an ellipsis

if nargin<2 || isempty(name)
    name = inputname(1);
    if isempty(name)
        name = 'Rop';
    end
end
if nargin<5 || isempty(call_name)
    call_name = name;                   % the name that the user can call
end
if strcmp(name,'ans')
    % Display an unnamed object under a default name
    name = 'Rop';
    fprintf('\n  ans = %s, where\n',name);
end
if nargin<3
    width_max = [];
end
if nargin<4
    expand = true;
end

if isempty(Rop.ops)
    fprintf('\n  %s is an empty tensor-PI operator.\n\n',name);
    return
end
if ~Rop.chkval
    try
        Rop = initialize(Rop);
    catch
        builtin('disp',Rop);
        return
    end
end
[K,d] = size(Rop.ops);

% Construct the terms: up to three, or the first and last, and an ellipsis
disp_struct = get_display(Rop,[],expand);

% An identity operator is displayed as (Rop x)(s) = x(s)
if disp_struct.is_id
    fprintf('\n');
    lhs = [name,' ',disp_struct.xsym];
    if ~isempty(disp_struct.out_args)
        lhs = ['(',lhs,')(',strjoin(disp_struct.out_args,','),')'];
    end
    print_terms([lhs,' ='],disp_struct.terms,1,width_max);
    fprintf('\n');
    return
end

lhs = ['(',name,' ',disp_struct.xsym,')'];
if ~isempty(disp_struct.out_args)
    lhs = [lhs,'(',strjoin(disp_struct.out_args,','),')'];
end
lhs = [lhs,' ='];
fprintf('\n');
print_terms(lhs,disp_struct.terms,disp_struct.signs,width_max);
if d>1 || ~disp_struct.expanded
    fprintf('\n');
end

% Indicate the input space, x(s_1,...,s_d) = x_1(s_1) o ... o x_d(s_d)
if d>1
    x_fctrs = cell(1,d);
    for j=1:d
        x_fctrs{j} = disp_struct.xnames{j};
        if ~isempty(disp_struct.fargs_in{j})
            x_fctrs{j} = [x_fctrs{j},'(',strjoin(disp_struct.fargs_in{j},','),')'];
        end
    end
    if d>max_factors
        x_fctrs = {x_fctrs{1},'...',x_fctrs{d}};
    end
    x_lhs = disp_struct.xsym;
    if ~isempty(disp_struct.in_args)
        x_lhs = [x_lhs,'(',strjoin(disp_struct.in_args,','),')'];
    end
    fprintf('  where %s = %s.\n',x_lhs,strjoin(x_fctrs,disp_struct.in_sym));
end

if ~disp_struct.expanded
    op_name = 'R';
    if K>1
        op_name = [op_name,char(7522)];   % subscript i
    end
    if d>1
        op_name = [op_name,char(11388)];  % subscript j
    end
    if K>1
        iname = 'i';
    else
        iname = '1';
    end
    if d>1
        jname = 'j';
    else
        jname = '1';
    end
    fprintf('  Call "%s.ops{%s,%s}" to extract operator %s.\n',call_name,iname,jname,op_name);
end
fprintf('\n');

end
