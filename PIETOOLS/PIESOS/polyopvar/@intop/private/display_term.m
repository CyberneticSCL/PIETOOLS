function [term_str,term_sign,use_K] = display_term(Kop,sub_idcs,k,xsym)
% [TERM_STR,TERM_SIGN,USE_K] = DISPLAY_TERM(KOP,SUB_IDCS,K,XSYM)
% builds a char representing the term of the functional KOP defined by the
% kernel KOP.PARAMS{K}, for use in the 'display' method of 'intop'.
%
% INPUTS
% - Kop:        'intop' object with N spatial variables;
% - sub_idcs:   1 x N array of integers, the subscripts of the kernel in
%               Kop.params (empty if N=0), so that
%               Kop.intmat{i}(sub_idcs(i),:) specifies the integrals
%               along direction i;
% - k:          integer, the linear index of the kernel in Kop.params,
%               used to name the kernel K_k;
% - xsym:       (optional) char, the name of the function that the
%               functional acts on, as in K(s) x(s). Defaults to 'x';
%
% OUTPUTS
% - term_str:   1 x n char representing the term, without its sign, e.g.
%                   0_int^1 0_int^s2 K3(s1,s2) x(s1,s2) ds1 ds2
%               (with the integral sign, subscripts and superscripts as
%               unicode characters). The integrals are ordered with the
%               outermost first, using get_int_order, with the limits
%               following from intmat. The lower limit precedes the
%               integral symbol, and the upper limit follows it. The limits
%               are displayed as subscripts and superscripts if they are
%               integers or variables, and as {x}_ and ^{x} otherwise,
%               e.g. for 0.5. In the limits, the comma between the first
%               and second index of a variable is dropped. A dummy variable
%               gamma is displayed as the subscript or superscript gamma in
%               a limit, followed by its number;
% - term_sign:  1 or -1, the sign of the term. Only -1 if the kernel is
%               displayed explicitly as a single term with a negative
%               coefficient;
% - use_K:      logical, true if the kernel is displayed as K_k(...), and
%               false if it is displayed explicitly;
%
% NOTES
% - The dummy variables are displayed as gamma, with the number of the
%   integrated copy as subscript if there are multiple (copies with equal
%   limits are multipliers, and not integrated). The domains are those of
%   Kop.dom.
% - Along a direction i, a variable sj with equal limits sp,sp is a
%   multiplier, evaluated at sp: it is not integrated, and sp is displayed
%   wherever sj would appear, also in the limits of other variables.
% - The kernel is displayed explicitly if it has at most 4 nonzeros in
%   Kop.params{k}, and Kop is scalar-valued (m=n=1). Decision variables are
%   displayed by their names Kop.Zd. The cost does not depend on their
%   number: only the nonzeros of Kop.params{k} are accessed.
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
max_nnz = 4;            % maximal number of nonzero coefficients of a kernel
                        % displayed explicitly

if nargin<4
    xsym = 'x';
end

nvars = numel(Kop.vars.in);
ncopies = Kop.depmat.in;
ndum = sum(ncopies);
dum_pos = [0,cumsum(ncopies)];               % copies of variable i: dum_pos(i)+1:dum_pos(i+1)
dom = Kop.dom.in;
[m,n] = size(Kop);

dum_names = cell(1,ndum);

% Resolve the multipliers: copy copy_idx along direction i is evaluated at the
% integrated copy rep_idx(copy_idx), or at the boundary value
% rep_val(copy_idx) if rep_idx(copy_idx)=0.
% A multiplier may be evaluated at another multiplier, so follow the chain.
rep_idx = zeros(1,ndum);       rep_val = zeros(1,ndum);
lim_lo_idx_all = cell(1,nvars);       lim_hi_idx_all = cell(1,nvars);
for i=1:nvars
    intmat_row = Kop.intmat{i}(sub_idcs(i),:);
    lim_lo_idx = intmat_row(1:2:end);   lim_hi_idx = intmat_row(2:2:end);
    lim_lo_idx_all{i} = lim_lo_idx;        lim_hi_idx_all{i} = lim_hi_idx;
    for copy_idx=1:ncopies(i)
        eval_idx = copy_idx;
        for it=1:ncopies(i)+1
            if lim_lo_idx(eval_idx)~=lim_hi_idx(eval_idx)
                break
            elseif lim_lo_idx(eval_idx)==0 || lim_lo_idx(eval_idx)==ncopies(i)+1
                rep_val(dum_pos(i)+copy_idx) = dom(i,1+(lim_lo_idx(eval_idx)>0));
                eval_idx = 0;
                break
            end
            eval_idx = lim_lo_idx(eval_idx);
        end
        if eval_idx>0
            rep_idx(dum_pos(i)+copy_idx) = dum_pos(i)+eval_idx;
        end
    end
end
% Names of the dummy variables, numbered by the kth integrated copy if there
% are multiple (the others are multipliers, evaluated at an integrated copy)
dum_vars = repelem(reshape(Kop.vars.in,1,[]),reshape(ncopies,1,[]));      % variable of each copy
is_int = rep_idx==(1:ndum);
nint = sum(is_int);
for dum_idx=find(is_int)
    dum_names{dum_idx} = varname2display([dum_vars{dum_idx},'_dum']);
    if nint>1
        dum_names{dum_idx} = varname2display(sprintf('%s_dum%d',dum_vars{dum_idx},sum(is_int(1:dum_idx))));
    end
end
eff_names = dum_names;           % name or value that each copy is displayed as
for dum_idx=1:ndum
    if rep_idx(dum_idx)==0
        eff_names{dum_idx} = sprintf('%g',rep_val(dum_idx));
    else
        eff_names{dum_idx} = dum_names{rep_idx(dum_idx)};
    end
end

int_sym = char(8747);
int_strs = cell(1,nvars);       diff_strs = cell(1,nvars);
for i=1:nvars
    lim_lo_idx = lim_lo_idx_all{i};        lim_hi_idx = lim_hi_idx_all{i};
    int_strs{i} = cell(1,0);        diff_strs{i} = cell(1,0);
    if ncopies(i)==0
        continue
    end
    lim_names = [{sprintf('%g',dom(i,1))},eff_names(dum_pos(i)+1:dum_pos(i+1)),{sprintf('%g',dom(i,2))}];
    % Integrate over the variables that are not multipliers, ordered
    % inner to outer by get_int_order, using the index to break ties
    int_order = get_int_order(Kop.intmat{i}(sub_idcs(i),:));
    int_copies = find(lim_lo_idx~=lim_hi_idx);
    [~,sort_idcs] = sortrows([int_order(int_copies).',int_copies.']);
    int_copies = int_copies(sort_idcs);
    nint_i = numel(int_copies);
    for j=1:nint_i
        copy_idx = int_copies(j);
        % The limits as subscript and superscript, or in braces, see
        % 'varname2display'. The lower limit precedes the integral symbol,
        % the upper limit follows it
        [~,lim_lo] = varname2display(lim_names{lim_lo_idx(copy_idx)+1});
        [~,~,lim_hi] = varname2display(lim_names{lim_hi_idx(copy_idx)+1});
        int_strs{i}{nint_i+1-j} = [lim_lo,int_sym,lim_hi];
        diff_strs{i}{j} = ['d',dum_names{dum_pos(i)+copy_idx}];
    end
end
int_strs_all = cell(1,0);       diff_strs_all = cell(1,0);
for i=1:nvars
    % The first variable is integrated innermost, so that the differentials
    % appear in order
    int_strs_all = [int_strs{i},int_strs_all];
    diff_strs_all = [diff_strs_all,diff_strs{i}];
end

param_k = Kop.params{k};
kern_args = strjoin(eff_names,',');
if nnz(param_k)<=max_nnz && m*n==1
    [kern_str,term_sign] = explicit_kernel(Kop,param_k,rep_idx,rep_val,dum_names);
    use_K = false;
else
    if isscalar(Kop.params)
        kern_str = 'K';
    else
        kern_str = varname2display(sprintf('K%d',k));
    end
    if nvars>0
        kern_str = [kern_str,'(',kern_args,')'];
    end
    term_sign = 1;
    use_K = true;
end

if nvars>0
    x_str = [xsym,'(',kern_args,')'];
else
    x_str = xsym;
end
term_parts = [int_strs_all,{kern_str,x_str},diff_strs_all];
term_parts = term_parts(~cellfun(@isempty,term_parts));
term_str = strjoin(term_parts,' ');

end
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %



%% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %%
function [kern_str,ksign] = explicit_kernel(Kop,param_k,rep_idx,rep_val,dum_names)
% [KERN_STR,KSIGN] = EXPLICIT_KERNEL(KOP,PARAM_K,REP_IDX,REP_VAL,
%   DUM_NAMES) builds a char representing the kernel defined by the
% sparse coefficients PARAM_K, as a polynomial in the dummy variables,
% with each multiplier replaced by the variable (or value) it is
% evaluated at.
%
% INPUTS
% - Kop:        'intop' object;
% - param_k:    m*(nd+1) x n*nZ sparse array of coefficients, as in
%               Kop.params;
% - rep_idx:    1 x D array of integers, with rep_idx(g)>0 the index of
%               the dummy variable that variable g is evaluated at, and
%               rep_idx(g)=0 if it is evaluated at the boundary value
%               rep_val(g);
% - rep_val:    1 x D array of doubles, see rep_idx;
% - dum_names:  1 x D cellstr, the names of the dummy variables;
%
% OUTPUTS
% - kern_str:   char representing the kernel. Matrix-valued kernels are
%               displayed as [K11, K12; K21, K22]. A scalar kernel with
%               more than one term is enclosed in brackets, and the
%               constant 1 is displayed as an empty char;
% - ksign:      1 or -1, the sign of a scalar kernel with a single term,
%               which is not included in kern_str, and 1 otherwise;

nvars = numel(Kop.Z.in);
ndum = numel(dum_names);
ndvars = numel(Kop.Zd);
[m,n] = size(Kop);
[nmons,strides,nZ] = Zdims(Kop.Z.in);
pm_chars = '+-';

% Determine the degree of each term in each (merged) dummy variable
[row_pos,col_pos,coefs] = find(param_k);
row_pos = row_pos(:);   col_pos = col_pos(:);   coefs = coefs(:);       % row vectors if P has one row
nt = numel(coefs);
blk_idx = ceil(row_pos/(ndvars+1));       dvar_idx = row_pos - (blk_idx-1)*(ndvars+1);    % dl=1: constant, dl=l+1: Zd{l}
col_blk = ceil(col_pos/nZ);            mono_pos = col_pos - (col_blk-1)*nZ;
degmat = zeros(nt,ndum);
for t=1:nt
    mono_idcs = mod(floor((mono_pos(t)-1)./strides),nmons) + 1;
    degs = zeros(1,0);
    for i=1:nvars
        degs = [degs,Kop.Z.in{i}(mono_idcs(i),:)];
    end
    for g=1:ndum
        if rep_idx(g)>0
            degmat(t,rep_idx(g)) = degmat(t,rep_idx(g)) + degs(g);
        else
            coefs(t) = coefs(t)*rep_val(g)^degs(g);
        end
    end
end
% Add the terms of equal monomials, which may arise from the merging
[term_keys,~,key_idx] = unique([blk_idx,col_blk,dvar_idx,degmat],'rows');
coef_sums = accumarray(key_idx,coefs);
is_nonzero = coef_sums~=0;
term_keys = term_keys(is_nonzero,:);        coef_sums = coef_sums(is_nonzero);

elem_strs = repmat({'0'},m,n);
elem_nterms = zeros(m,n);
elem_signs = ones(m,n);
for t=1:numel(coef_sums)
    row_t = term_keys(t,1);         col_t = term_keys(t,2);     dvar_t = term_keys(t,3);
    mono_degs = term_keys(t,4:end);
    mono_facs = cell(1,0);
    for g=find(mono_degs>0)
        mono_facs{end+1} = dum_names{g};                             %#ok<AGROW>
        if mono_degs(g)>1
            [~,~,pow_str] = varname2display(sprintf('%d',mono_degs(g)));       % superscript
            mono_facs{end} = [dum_names{g},pow_str];
        end
    end
    if dvar_t>1
        mono_facs = [Kop.Zd(dvar_t-1).',mono_facs];
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
        kern_str = ['(',repmat('-',1,elem_signs(1)<0),kern_str,')'];
        ksign = 1;
    elseif strcmp(kern_str,'1')
        kern_str = '';
    end
else
    ksign = 1;
    row_strs = cell(1,m);
    for i=1:m
        for j=1:n
            if elem_nterms(i,j)>0
                elem_strs{i,j} = [repmat('-',1,elem_signs(i,j)<0),elem_strs{i,j}];
            end
        end
        row_strs{i} = strjoin(elem_strs(i,:),', ');
    end
    kern_str = ['[',strjoin(row_strs,'; '),']'];
end

end
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
