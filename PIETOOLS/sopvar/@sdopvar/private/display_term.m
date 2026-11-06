function [term_str,term_sign,use_K] = display_term(P,k,xsym,out_names,entry_idcs)
% [TERM_STR,TERM_SIGN,USE_K] = DISPLAY_TERM(P,K,XSYM,OUT_NAMES,ENTRY_IDCS)
% builds a char representing the term of the operator P defined by the
% parameters P.PARAMS.A{K} and P.PARAMS.B{K}, for use in the 'display' and
% 'get_display' methods of 'sdopvar'.
%
% INPUTS
% - P:      'sdopvar' object, with M output variables and N input
%           variables;
% - k:      integer, the linear index of the parameter in P.params.A, used
%           to name the kernel K_k;
% - xsym:   (optional) char, the name of the function that the operator
%           acts on, as in K(s,s') x(s'). Defaults to 'x'. If the operator
%           acts on R, and xsym is '1', the function is not displayed;
% - out_names:  (optional) 1 x M cellstr, the names to display for the
%               output variables. Defaults to those given by
%               'varname2display';
% - entry_idcs: (optional) 1 x 2 array, the row and column of the element
%               of the matrix-valued kernel. If given, the term is that of
%               the scalar kernel K_ij, acting on the function XSYM (the
%               jth element of the function that P acts on), e.g. for the
%               ith row of P x;
%
% OUTPUTS
% - term_str:   1 x n char representing the term, without its sign, e.g.
%                   0_int^s K2(s,g) x(g) dg
%               (with the integral sign, subscripts and superscripts as
%               unicode characters). The integrals are those in the
%               variables P.vars.in, in order, with the outermost first.
%               The dummy variables are displayed as gamma, with the number
%               of the integrated variable as subscript if there are
%               multiple. The limits are displayed as subscripts and
%               superscripts if they are integers or variables, and as {x}_
%               and ^{x} otherwise, e.g. for 0.5;
% - term_sign:  1 or -1, the sign of the term. Only -1 if the kernel is
%               displayed explicitly as a single term with a negative
%               coefficient;
% - use_K:      logical, true if the kernel is displayed as K_k(...), and
%               false if it is displayed explicitly;
%
% NOTES
% - A variable that is also an output variable has the parameter indexed by
%   1, 2 or 3, see 'decode_params'. For 1, the term is a multiplier, so
%   that there is no integral, and the function is evaluated at the output
%   variable, which is also used as the limit of the integrals for 2 and 3.
%   For 2, the integral is taken from the lower boundary of the domain to
%   the output variable, and for 3, from the output variable to the upper
%   boundary. A variable that is only an input variable is integrated over
%   its domain.
% - The kernel depends on the output variables, on the dummy variables that
%   are integrated, and on the decision variables P.Zd, which are
%   displayed by their names. It is displayed explicitly if it has at most
%   4 nonzero coefficients in P.params.A{k} and P.params.B{k} together,
%   and P is scalar-valued (P.dims=[1,1]), or if ENTRY_IDCS is given, and
%   the element K_ij has at most 4 nonzero coefficients. The cost does not
%   depend on the number of decision variables.
% - The characters are declared by their code, to keep the file ASCII.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PIETOOLS - display_term
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
% DJ, 11/06/2026: Get the names of dummy variables, the limits, and the
%                   subscripts from 'varname2display';

% Default settings
max_nnz = 4;            % maximal number of nonzero coefficients of a kernel (element)
                        % displayed explicitly

if nargin<3 || isempty(xsym)
    xsym = 'x';
end
if nargin<5
    entry_idcs = [];
end

out_vars = reshape(P.vars.out,1,[]);
in_vars = reshape(P.vars.in,1,[]);
nvars_in = numel(in_vars);
m = P.dims(1);
n = P.dims(2);
if nargin<4 || isempty(out_names)
    out_names = cellfun(@varname2display,out_vars,'UniformOutput',false);
end
out_disp = reshape(out_names,1,[]);

% Type of each input variable: 0 if it is only an input variable, and 1, 2
% or 3 if it is also an output variable, as for the parameter
shared_vars = intersect(in_vars,out_vars);
multi_idx = ones(1,numel(shared_vars));
if ~isempty(shared_vars)
    multi_idx = gamma_of_cell(k,numel(shared_vars));
end
var_kind = zeros(1,nvars_in);
[is_shared,shared_pos] = ismember(in_vars,shared_vars);
var_kind(is_shared) = multi_idx(shared_pos(is_shared));
[~,out_pos] = ismember(in_vars,out_vars);        % position of each input variable among the output variables

% Dummy variables, numbered by the kth integrated variable if there are
% multiple (the dummy variable of a multiplier is not displayed)
dum_disp = cellfun(@(v) varname2display([v,'_dum']),in_vars,'UniformOutput',false);
if sum(var_kind~=1)>1
    int_rank = cumsum(var_kind~=1);
    for q=1:nvars_in
        dum_disp{q} = varname2display(sprintf('%s_dum%d',in_vars{q},int_rank(q)));
    end
end

% Declare the integrals, outermost first, the differentials, innermost
% first, the arguments of the function, and of the kernel
int_strs = cell(1,0);       diff_strs = cell(1,0);
x_args = cell(1,nvars_in);      kern_args = out_disp;
for q=1:nvars_in
    if var_kind(q)==1
        % Multiplier: the function is evaluated at the output variable
        x_args{q} = out_disp{out_pos(q)};
        continue
    end
    x_args{q} = dum_disp{q};
    kern_args{end+1} = dum_disp{q};                                %#ok<AGROW>
    % Limits: the boundaries of the domain, or the output variable, as
    % displayed. A limit is a subscript or superscript, or in braces, see
    % 'varname2display'
    [~,lim_lo] = varname2display(sprintf('%g',P.dom.in(q,1)));
    [~,~,lim_hi] = varname2display(sprintf('%g',P.dom.in(q,2)));
    if var_kind(q)==2
        [~,~,lim_hi] = varname2display(out_disp{out_pos(q)});
    elseif var_kind(q)==3
        [~,lim_lo] = varname2display(out_disp{out_pos(q)});
    end
    int_strs{end+1} = [lim_lo,char(8747),lim_hi];                       %#ok<AGROW>
    diff_strs = [{['d',dum_disp{q}]},diff_strs];                       %#ok<AGROW>
end

% Display the kernel, explicitly if it is small. The terms are those of the
% constant part A and of the part B that multiplies the decision variables
nnz_k = nnz(P.params.A{k}) + nnz(P.params.B{k});
if ~isempty(entry_idcs)
    [row_idcs,col_idcs,dvar_idcs,coefs,degmat] = decode_params(P,k);
    entry_sel = row_idcs==entry_idcs(1) & col_idcs==entry_idcs(2);              % terms of the element K_ij
end
if nnz_k<=max_nnz && m*n==1
    [row_idcs,col_idcs,dvar_idcs,coefs,degmat] = decode_params(P,k);
    [kern_str,term_sign] = explicit_kernel(row_idcs,col_idcs,dvar_idcs,coefs,degmat,[out_disp,dum_disp],P.Zd,m,n);
    use_K = false;
elseif ~isempty(entry_idcs) && nnz(entry_sel)<=max_nnz
    % Element K_ij of a matrix-valued kernel, as a scalar kernel
    nsel = nnz(entry_sel);
    [kern_str,term_sign] = explicit_kernel(ones(nsel,1),ones(nsel,1),dvar_idcs(entry_sel),coefs(entry_sel),degmat(entry_sel,:),[out_disp,dum_disp],P.Zd,1,1);
    use_K = false;
else
    if isscalar(P.params.A)
        kern_str = 'K';
    else
        kern_str = varname2display(sprintf('K%d',k));
    end
    if ~isempty(kern_args)
        kern_str = [kern_str,'(',strjoin(kern_args,','),')'];
    end
    term_sign = 1;
    use_K = true;
end

% Combine the integrals, kernel, function, and differentials
if nvars_in>0
    x_str = [xsym,'(',strjoin(x_args,','),')'];
elseif strcmp(xsym,'1')
    x_str = '';                  % the constant 1 is not displayed
else
    x_str = xsym;
end
term_parts = [int_strs,{kern_str,x_str},diff_strs];
term_parts = term_parts(~cellfun(@isempty,term_parts));
term_str = strjoin(term_parts,' ');

end
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %



%% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %%
function [kern_str,ksign] = explicit_kernel(i,j,l,coefs,degmat,var_names,dvar_names,m,n)
% [KERN_STR,KSIGN] = EXPLICIT_KERNEL(I,J,L,COEFS,DEGMAT,VAR_NAMES,
%   DVAR_NAMES,M,N) builds a char representing the kernel defined by the
% terms coefs(t)*(decision variable l(t))*(monomial degmat(t,:)) in
% element (I(t),J(t)) of the m x n kernel.
%
% INPUTS
% - i,j:        nt x 1 arrays, the row and column of each term;
% - l:          nt x 1 array, the index of the decision variable of each
%               term, or 0 for a constant;
% - coefs:      nt x 1 array, the coefficient of each term;
% - degmat:     nt x D array, the degrees of each term in the variables;
% - var_names:  1 x D cellstr, the names of the variables as displayed;
% - dvar_names: cell, the names of the decision variables;
% - m,n:        integers, the dimensions of the kernel;
%
% OUTPUTS
% - kern_str:  char representing the kernel. Matrix-valued kernels are
%               displayed as [K11, K12; K21, K22]. A scalar kernel with
%               more than one term is enclosed in brackets, and the
%               constant 1 is displayed as an empty char. A decision
%               variable is displayed before the monomial, e.g. d1*s*g;
% - ksign:  1 or -1, the sign of a scalar kernel with a single term, which
%           is not included in kern_str, and 1 otherwise;

pm_chars = '+-';                                    % pm_chars(1+(x<0)) is the sign of x

% Add the terms of equal monomials, which may arise from multipliers
[term_keys,~,key_idx] = unique([i,j,l,degmat],'rows');
coef_sums = accumarray(key_idx,coefs);
is_nonzero = coef_sums~=0;
term_keys = term_keys(is_nonzero,:);        coef_sums = coef_sums(is_nonzero);

% Display the terms of each element
elem_strs = repmat({'0'},m,n);
elem_nterms = zeros(m,n);
elem_signs = ones(m,n);
for t=1:numel(coef_sums)
    row_t = term_keys(t,1);         col_t = term_keys(t,2);
    dvar_t = term_keys(t,3);        mono_degs = term_keys(t,4:end);
    mono_facs = cell(1,0);
    if dvar_t>0
        mono_facs{end+1} = char(dvar_names{dvar_t});                 %#ok<AGROW>
    end
    for g=find(mono_degs>0)
        mono_facs{end+1} = var_names{g};                             %#ok<AGROW>
        if mono_degs(g)>1
            [~,~,pow_str] = varname2display(sprintf('%d',mono_degs(g)));       % superscript
            mono_facs{end} = [var_names{g},pow_str];
        end
    end
    if abs(coef_sums(t))~=1 || isempty(mono_facs)
        mono_facs = [{sprintf('%.4g',abs(coef_sums(t)))},mono_facs];
    end
    mono_str = strjoin(mono_facs,'*');
    if elem_nterms(row_t,col_t)==0
        elem_strs{row_t,col_t} = mono_str;
        elem_signs(row_t,col_t) = sign(coef_sums(t));
    else
        elem_strs{row_t,col_t} = [elem_strs{row_t,col_t},' ',pm_chars(1+(coef_sums(t)<0)),' ',mono_str];
    end
    elem_nterms(row_t,col_t) = elem_nterms(row_t,col_t)+1;
end

if m*n==1
    kern_str = elem_strs{1};
    ksign = elem_signs(1);
    if elem_nterms(1)>1
        % Sum of terms: the first term keeps its sign
        kern_str = ['(',repmat('-',1,elem_signs(1)<0),kern_str,')'];
        ksign = 1;
    elseif strcmp(kern_str,'1')
        kern_str = '';
    end
else
    ksign = 1;
    row_strs = cell(1,m);
    for row_t=1:m
        for col_t=1:n
            if elem_nterms(row_t,col_t)>0
                elem_strs{row_t,col_t} = [repmat('-',1,elem_signs(row_t,col_t)<0),elem_strs{row_t,col_t}];
            end
        end
        row_strs{row_t} = strjoin(elem_strs(row_t,:),', ');
    end
    kern_str = ['[',strjoin(row_strs,'; '),']'];
end

end
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
