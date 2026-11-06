function disp_struct = get_display(Rop,xnames,expand)
% DISP_STRUCT = GET_DISPLAY(ROP,XNAMES,EXPAND) constructs the terms of the
% expression
%   (R_11 x_1)(s) o ... o (R_1d x_d)(s) + ...
% of the tensor-PI operator ROP acting on the functions x_1,...,x_d, with
% o the tensor product, or the actual operators R_ij acting on x_j, if they
% are small enough, for use in the 'display' method of 'tensopvar', and in
% the 'get_display' methods of the classes of operators that have
% 'tensopvar' objects as blocks.
%
% INPUTS
% - Rop:    'tensopvar' object, with K terms and d factors in each term;
% - xnames: (optional) char, or 1 x d cellstr, the names of the functions
%           on which the factors act. A char is used for all factors.
%           Defaults to x_j for factor j, or to x if d=1, with x replaced
%           by u, w, or v, the first that is not the name of a variable, if
%           a variable is named x (ignoring a number at the end of the
%           name);
% - expand: (optional) logical, true if the operators R_ij may be displayed
%           as the actual operators if they are small enough. If false,
%           they are always displayed as R_ij. Defaults to true;
%
% OUTPUTS
% - disp_struct:      struct with the fields
%       terms:  1 x n cell of char, the terms of the sum, without their
%               sign. If there are more than 3 terms, and the operators
%               R_ij are not displayed, only the first and last term are
%               constructed, with '...' in between;
%       signs:  1 x n array of 1 and -1, the sign of each term (1 for
%               '...');
%       expanded:   logical, true if the operators R_ij are displayed as
%               the actual operators, and false if they are displayed as
%               R_ij, in which case they can be extracted by calling
%               Rop.ops{i,j};
%       xsym:   char, the name of the function x that Rop acts on, with
%               x = x_1 o ... o x_d, if d>1;
%       xnames: 1 x d cellstr, the names of the functions on which the
%               factors act;
%       out_args:   cellstr, the names of the output variables of Rop, with
%               a subscript for the copy of a variable if it has multiple
%               copies;
%       in_args:    cellstr, the same for the input variables;
%       fargs_in:   1 x d cell of cellstr, the names of the input variables
%               of each of the factors;
%       in_sym:     char, the symbol that joins the factors x_j of x, which
%               is o for a tensor product, and = for a pointwise product;
%       is_id:  logical, true if Rop is an identity operator, see
%               'is_id_op'. The only term is then the function with its
%               spatial variables, e.g. x(s), and expanded is true;
%
% NOTES
% - The operators R_ij are displayed as the actual operators if there are
%   at most 3 terms, all R_ij can be displayed explicitly by 'get_display'
%   (which requires them to be scalar-valued), and the sum has at
%   most 6 terms of at most 100 characters. If not, none of the operators
%   is displayed.
% - The output variables of the factor R_ij, which can be given names with
%   a subscript for the copy of the variable, are displayed as such in the
%   display of R_ij. The dummy variables of the integrals in R_ij are
%   displayed as gamma.

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
% DJ, 11/06/2026: Get the subscripts from 'varname2display';

% Default settings
max_ops = 3;            % maximal number of terms for which the operators R_ij are
                        % displayed as the actual operators, and above which only
                        % the first and last term are displayed
max_terms = 6;          % maximal number of terms displayed with the actual operators
max_chars = 100;        % maximal number of characters of each such term

if nargin<3
    expand = true;
end
if ~Rop.chkval
    Rop = initialize(Rop);
end
[K,d] = size(Rop.ops);

otimes_sym = char(8855);
if Rop.otype(1)==2
    out_sym = [' ',otimes_sym,' '];
else
    out_sym = ' .* ';
end
if Rop.otype(2)==2
    in_sym = [' ',otimes_sym,' '];
else
    in_sym = ' = ';
end
[~,sub_dig] = varname2display('0123456789');                 % subscript digits
sub_str = @(n) sub_dig(sprintf('%d',n)-'0'+1);

% Name of the input, x, unless this is also the name of a variable (ignoring a
% number at the end of the name), in which case the first of u, w, v that is not
var_bases = regexprep([Rop.vars.out,Rop.vars.in],'_?\d+$','');
for xsym = 'xuwv'
    if ~ismember(xsym,var_bases)
        break
    end
end
if nargin>=2 && ~isempty(xnames)
    if ischar(xnames)
        xnames = repmat({xnames},1,d);
    end
else
    xnames = cell(1,d);
    for j=1:d
        xnames{j} = xsym;
        if d>1
            xnames{j} = [xsym,sub_str(j)];
        end
    end
end

% Determine the arguments of the output (side 1) and the input (side 2) of
% each factor: the variables it maps to or from, with a subscript for the
% copy if the variable has multiple copies. With a tensor product, the
% copies are numbered over the factors in order
var_names = {Rop.vars.out,Rop.vars.in};
depmat_sides = {Rop.depmat.out,Rop.depmat.in};       % d x M, d x N: copies of each variable in each factor
deg_sides = {Rop.deps.out,Rop.deps.in};           % 1 x M, 1 x N: copies of each variable
fctr_args_side = cell(1,2);     args_side = cell(1,2);      copy_offs_side = cell(1,2);
for s=1:2
    nvars_side = numel(var_names{s});
    dep_side = depmat_sides{s};       deg_side = deg_sides{s};
    if Rop.otype(s)==2
        copy_offs = cumsum(dep_side,1) - dep_side;
    else
        copy_offs = zeros(d,nvars_side);
    end
    copy_offs_side{s} = copy_offs;
    fctr_args = cell(1,d);
    for j=1:d
        fctr_args{j} = cell(1,0);
        for k=1:nvars_side
            for c=1:dep_side(j,k)
                if deg_side(k)>1
                    op_var_names = varname2display(var_names{s}{k},copy_offs(j,k)+c);
                else
                    op_var_names = varname2display(var_names{s}{k});
                end
                fctr_args{j}{end+1} = op_var_names;
            end
        end
    end
    fctr_args_side{s} = fctr_args;
    if Rop.otype(s)==2
        args_side{s} = [fctr_args{:}];
    else
        % Pointwise product: all factors use the same copies of each variable
        args_side{s} = cell(1,0);
        for k=1:nvars_side
            for c=1:deg_side(k)
                if deg_side(k)>1
                    op_var_names = varname2display(var_names{s}{k},c);
                else
                    op_var_names = varname2display(var_names{s}{k});
                end
                args_side{s}{end+1} = op_var_names;
            end
        end
    end
end
disp_struct = struct('terms',{cell(1,0)},'signs',zeros(1,0),'expanded',false,'xsym',xsym,'xnames',{xnames}, ...
    'out_args',{args_side{1}},'in_args',{args_side{2}},'fargs_in',{fctr_args_side{2}},'in_sym',in_sym,'is_id',false);
disp_struct.is_id = K==1 && is_id_op(Rop);
if disp_struct.is_id
    % The function itself, with the spatial variables, x(s) (x for d>1: the
    % product of the x_j)
    x_id_str = xnames{1};
    if d>1
        x_id_str = xsym;
    end
    if ~isempty(args_side{1})
        x_id_str = [x_id_str,'(',strjoin(args_side{1},','),')'];
    end
    disp_struct.terms = {x_id_str};        disp_struct.signs = 1;        disp_struct.expanded = true;
    return
end

% Display the operators R_ij as the actual operators, if there are at most 3
% terms, and these are small enough
terms = cell(1,0);      signs = zeros(1,0);
expanded = expand && K<=max_ops;
for i=1:K*expanded
    fctr_strs = cell(1,d);
    for j=1:d
        op = Rop.ops{i,j};
        if ~isa(op,'sopvar')
            expanded = false;
            break
        end
        % The names of the output variables of the factor, with subscripts for
        % the copies of the variables. The input variables are not displayed
        out_var_names = cell(1,numel(op.vars.out));
        for r=1:numel(op.vars.out)
            k = find(strcmp(var_names{1},op.vars.out{r}),1);
            if ~isempty(k) && deg_sides{1}(k)>1
                out_var_names{r} = varname2display(op.vars.out{r},copy_offs_side{1}(j,k)+1);
            else
                out_var_names{r} = varname2display(op.vars.out{r});
            end
        end
        fctr_disp = get_display(op,xnames{j},out_var_names);
        if ~fctr_disp.ok
            expanded = false;
            break
        end
        if isempty(fctr_disp.terms)
            fctr_strs{j} = '0';
        elseif numel(fctr_disp.terms)>1 || fctr_disp.signs(1)<0
            % A sum, or a negative term, is enclosed in brackets
            fctr_strs{j} = ['(',repmat('- ',1,fctr_disp.signs(1)<0),fctr_disp.terms{1}];
            for r=2:numel(fctr_disp.terms)
                fctr_strs{j} = [fctr_strs{j},' ',repmat('+',1,fctr_disp.signs(r)>0),repmat('-',1,fctr_disp.signs(r)<0),' ',fctr_disp.terms{r}];
            end
            fctr_strs{j} = [fctr_strs{j},')'];
        else
            fctr_strs{j} = fctr_disp.terms{1};
        end
    end
    if ~expanded
        break
    elseif d==1
        % A single factor: its terms are the terms of the sum
        terms = [terms,fctr_disp.terms];               signs = [signs,fctr_disp.signs];       %#ok<AGROW>
    else
        terms{end+1} = strjoin(fctr_strs,out_sym);   signs(end+1) = 1;               %#ok<AGROW>
    end
end
if expanded && (numel(terms)>max_terms || any(cellfun(@numel,terms)>max_chars))
    expanded = false;
end

% Otherwise, display the operators as R_ij acting on x_j, with the output
% variables, for the first and last term only if there are more than 3
if ~expanded
    terms_shown = 1:K;
    if K>max_ops
        terms_shown = [1,K];
    end
    terms = cell(1,numel(terms_shown));
    for ii=1:numel(terms_shown)
        i = terms_shown(ii);
        fctr_strs = cell(1,d);
        for j=1:d
            op_name = 'R';
            if K>1
                op_name = [op_name,sub_str(i)];
            end
            if K>1 && d>1 && max(i,j)>9
                op_name = [op_name,','];
            end
            if d>1
                op_name = [op_name,sub_str(j)];
            end
            fctr_strs{j} = ['(',op_name,' ',xnames{j},')'];
            if ~isempty(fctr_args_side{1}{j})
                fctr_strs{j} = [fctr_strs{j},'(',strjoin(fctr_args_side{1}{j},','),')'];
            end
        end
        terms{ii} = strjoin(fctr_strs,out_sym);
    end
    signs = ones(1,numel(terms));
    if K>max_ops
        terms = [terms(1),{'...'},terms(2)];
        signs = ones(1,3);
    end
end
disp_struct.terms = terms;
disp_struct.signs = signs;
disp_struct.expanded = expanded;

end
