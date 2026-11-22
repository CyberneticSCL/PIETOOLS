function KCop = int(Cop,W,vars,dom)
% KCOP = INT(COP,W,VARS,DOM) takes the composition of an integral
% (functional) with a tensor-PI operator COP,
%   KCop * x = int_{Omega} W*Cop*x dOmega
% where Omega is the output domain of COP, or the domain DOM of the
% variables VARS.
% KCOP = INT(COP,W,INTMAT) integrates over the region INTMAT of the copies
% of the output variables, given as the intmat of an 'intop'.
%
% INPUTS
% - Cop:    p x n 'tensopvar' object representing a tensor-PI operator
%           acting on a single degree-d distributed monomial;
% - W:      (optional) q x p 'double', 'polynomial' or 'dpvar' weight
%           inside the integral, defaults to the identity. If W depends on
%           the output variables of Cop, it is integrated along with Cop,
%           and may depend only on the variables
%             'si_j'  for the jth copy of vars.out{i}, j=1,...,deps.out(i),
%                     for a tensor product of the
%                     outputs (Cop.otype(1)==2, or degree above 1);
%             'si'    also for the only copy, if the ith output variable
%                     vars.out{i} has degree 1 (deps.out(i)==1) and
%                     Cop.otype(1)==1 (a pointwise product);
%           where si is the name vars.out{i}. Other variables are rejected,
%           as is a weight depending on both 'si' and 'si_1', unless they
%           are in VARS;
% - vars:   (optional) 1 x N 'cellstr' or 'polynomial' (pvar) array, the
%           variables to integrate over: each copy of the output variables
%           of Cop, named as for W, once, and any other variables. W is
%           integrated over (or evaluated at, see dom) the other variables
%           first. Defaults to the copies of the output variables;
% - dom:    (optional) N x 2 'double', 'polynomial' or 'cell' array, the
%           limits [L,U] of each variable in VARS. A limit is a number, or
%           a variable in VARS (a pvar, or its name in a cell, named as in
%           VARS). For a copy of an output variable, each limit is the
%           boundary of its domain (L=a, U=b) or another copy of the same
%           output variable; for the other variables, any number or
%           variable. Equal limits that are a variable, [t,t], evaluate at
%           t (a Dirac delta), e.g. [s_1,s_1]; for a copy, t must be
%           another copy of the same output variable. Equal numbers are an
%           interval of length zero, not an evaluation. Defaults to the
%           domain of Cop, so it must be given if VARS has other variables.
%           See NOTES;
% - intmat: (in the form int(Cop,W,intmat)) 1 x M cell, for each output
%           variable si of Cop a 1 x 2*Di array, the region of the copies
%           si_1,...,si_Di as a row of the intmat of an 'intop' (0 for ai,
%           Di+1 for bi, j for si_j, [j,j] for evaluation at si_j). The
%           first deps.out(i) copies are those of the output of Cop, and
%           lie between ai, bi and each other; W is integrated over the
%           other copies first, as over the other variables of VARS;
%
% OUTPUTS
% - KCop:   q x n 'intop' object defining a functional representing the
%           composition of an integral operator and Cop;
%
% NOTES
% - The variables "vars" must include all variables on which the output of
% Cop depends, integration with respect to only a subset of the variables
% is not supported.
% - Variable limits define an iterated integral, e.g. for the copies s_1,
%   s_2 of s in [a,b],
%     int(Cop,W,[s_1;s_2],[a,b;a,s_1])
%       = int_a^b int_a^{s_1} W(s_1,s_2)*(Cop*x)(s_1,s_2) ds_2 ds_1,
%   the integral over the triangle s_2<=s_1, an 'intop' with intmat
%   [0,3,0,1]. With dom = [a,b;s_1,s_1], it is the integral over the
%   diagonal, int_a^b W(s_1,s_1)*(Cop*x)(s_1,s_1) ds_1 (intmat [0,3,1,1]).
%   A variable may not appear in its own limits, also not through the
%   limits of the variables in its limits. The limits of the copies may
%   not be other variables, and a copy between two copies must have these
%   ordered by the other limits (see intop.canonicalize_intmat), so that
%   the iterated integral is the integral over the region.
%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PIETOOLS - int
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
% For support, contact M. Peet, Arizona State University at mpeet@asu.edu
% DJ, 10/15/2026: Initial coding;
% DJ, 11/22/2026: Variable limits of integration, equal limits for
%                   evaluation (Dirac delta), and regions as an intmat;

if nargin<2
    W = [];
end
% Check that the tensor-PI operator is properly declared
if ~Cop.chkval
    Cop = initialize(Cop);
end

% Variables of integration on which the output of Cop does not depend:
% integrate W over them, so that KCop integrates over all of VARS. The
% copies of the output variables are all integrated over below, over the
% domain of Cop
% A region given as an intmat (a cell of numeric rows) is read directly,
% with the same outputs as dom2intmat
is_intmat = nargin==3 && iscell(vars) && all(cellfun(@isnumeric,vars));
if is_intmat
    [vars,lims,intmat_lims,ord_W] = intmat2lims(vars,Cop.vars.out,Cop.dom.out,Cop.deps.out);
elseif nargin>=3
    if nargin<4
        dom = [];
    end
    % Copies of the output variables, named as for W (see below), with the
    % index of the copy each name refers to and its domain
    names_cp = cell(1,0);       idcs_cp = zeros(1,0);       dom_cp = zeros(0,2);
    ncp = 0;
    for i=1:numel(Cop.vars.out)
        d_i = Cop.deps.out(i);
        names_i = get_tensorvars(Cop.vars.out{i},d_i,'full');
        idcs_i = ncp+(1:d_i);
        if d_i==1 && Cop.otype(1)==1
            names_i = [names_i,get_tensorvars(Cop.vars.out{i},1,{'bare','full'})];
            idcs_i = [idcs_i,ncp+1];
        end
        names_cp = [names_cp,names_i];
        idcs_cp = [idcs_cp,idcs_i];
        dom_cp = [dom_cp;repmat(Cop.dom.out(i,:),numel(idcs_i),1)];
        ncp = ncp+d_i;
    end
    % The limits of the copies set the intmat of Kop (below). W is
    % integrated over the other variables, or evaluated at equal variable
    % limits, each before the variables in its limits (order ord_W)
    [vars,lims,intmat_lims,ord_W] = dom2intmat(vars,dom,names_cp,idcs_cp,dom_cp,Cop.deps.out); % DJ, 11/19/2026
end
if nargin>=3
    if isempty(W) && ~isempty(ord_W)
        W = speye(size(Cop,1));
    end
    for k=ord_W
        [L,U] = lims{k,:};
        if ~isnumeric(L) && isequal(L,U)                                    % DJ, 11/19/2026
            % Dirac delta: W at vars{k}=L, a constant W is unchanged
            if ~isa(W,'double')
                W = subs(W,polynomial(vars(k)),L);
            end
        elseif isa(W,'double') && isnumeric(L) && isnumeric(U)
            W = (U-L)*W;
        else
            % W a polynomial or dpvar, or a variable limit: integrate
            W = int(W,polynomial(vars(k)),L,U);
        end
    end
end

% Extract the output function space of Cop
m = size(Cop,1);
M = numel(Cop.vars.out);
vars_out = Cop.vars.out;
dom_out = Cop.dom.out;
degs_out = Cop.deps.out;

% Declare an integral on the output function space
intmat = cell(1,M);
for i=1:M
    intmat{i} = repmat([0,degs_out(i)+1],1,degs_out(i));
end
if nargin>=3
    % With the limits of dom, checked to be ordered by dom2intmat, so that
    % intop need not check them
    intmat = intmat_lims;                                                   % DJ, 11/19/2026
end
params = {speye(m)};
Kop = intop(params,intmat,vars_out,dom_out,true);

% For weights depending on variables: replace the identity by W, with its
% variables renamed to the dummy variables of integration
has_spatial_W = ~isempty(W) && ~isa(W,'double') && ~isempty(W.varname);
if has_spatial_W
    % Variables allowed in W, and the dummy variables of Kop they are
    % renamed to: the copies <name>_<j>, as named by apply_op for a tensor
    % product, and, for the only copy of a pointwise output, also the bare
    % name
    dum_vars = Kop.varnames.in_dum;
    names_W = cell(1,0);        names_dum = cell(1,0);
    for i=1:M
        names_i = get_tensorvars(vars_out{i},degs_out(i),'full');
        dum_i = reshape(dum_vars{i},1,[]);
        if degs_out(i)==1 && Cop.otype(1)==1
            names_i = [names_i,get_tensorvars(vars_out{i},1,{'bare','full'})];
            dum_i = [dum_i,dum_vars{i}(1)];
        end
        names_W = [names_W,names_i];
        names_dum = [names_dum,dum_i];
    end
    [is_ok,idcs] = ismember(W.varname,names_W);
    if ~all(is_ok)
        error("The weight may only depend on the variables 'si_j' (copy j of the output variable si) or, for one copy of a pointwise output, 'si'.")
    end
    % Rename all variables at once: a copy may not be renamed twice
    new_names = names_dum(idcs);
    if numel(unique(new_names))<numel(new_names)
        error("The weight may not depend on both 'si' and 'si_1', the same variable.")
    end
    W.varname = reshape(new_names,size(W.varname));
    Kop = intop({W},intmat,vars_out,dom_out,true);
end

% Take the composition of the integral with Cop
KCop = Kop*Cop;
if ~isempty(W) && ~has_spatial_W
    KCop = W*KCop;
end

end
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %



%% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %%
function [vars,lims,intmat,ord_W] = dom2intmat(vars,dom,names_cp,idcs_cp,dom_cp,deps_out)
% [VARS,LIMS,INTMAT,ORD_W] = DOM2INTMAT(VARS,DOM,NAMES_CP,IDCS_CP,DOM_CP,
% DEPS_OUT) reads the variables and limits of integration of 'int', where
% a limit may be a variable of integration, and returns the intervals of
% the copies of the output variables as the intmat of an 'intop', and the
% order in which W is integrated over the other variables.
%
% INPUTS
% - vars:       variables of integration, as in 'int';
% - dom:        limits of integration, as in 'int', or [] for the domain of
%               the copies;
% - names_cp:   1 x P 'cellstr', the names of the copies of the output
%               variables (and the bare name of a single pointwise copy);
% - idcs_cp:    1 x P array, the copy (1,...,ncp) each name refers to, the
%               copies ordered by output variable and then by copy;
% - dom_cp:     P x 2 array, the domain of each name;
% - deps_out:   1 x M array, the number of copies of each output variable;
%
% OUTPUTS
% - vars:       1 x N 'cellstr', the names of the variables of integration;
% - lims:       N x 2 cell, the limits of each variable: a 'double', or the
%               'polynomial' pvar of a variable of VARS;
% - intmat:     1 x M cell, the intervals of the copies of each output
%               variable, as the field intmat of an 'intop' (0 and d+1 for
%               the boundaries, j for copy j, [j,j] for evaluation at j),
%               with ordered rows (see intop.canonicalize_intmat);
% - ord_W:      1 x K array, the indices in VARS of the other variables, in
%               the order in which W is integrated over them: each before
%               the variables in its limits;

% Limits that are variables (a pvar, or a name) are taken out, and replaced
% by 0 for parse_int_inputs, which reads the numbers
if isa(dom,'polynomial')
    dom = arrayfun(@(k) dom(k),reshape(1:numel(dom),size(dom)),'UniformOutput',false);
elseif isnumeric(dom)
    dom = num2cell(dom);
end
names_lim = repmat({''},size(dom));
for k=1:numel(dom)
    if ischar(dom{k}) || isstring(dom{k})
        names_lim{k} = char(dom{k});            dom{k} = 0;
    elseif isa(dom{k},'polynomial') && ~isempty(dom{k}.varname)
        names_lim(k) = pvar2varname(dom{k});    dom{k} = 0;
    end
end
[vars,lims] = parse_int_inputs(vars,dom,names_cp,dom_cp);
N = numel(vars);
% idcs_lim(k,j): the index in VARS of limit j of variable k, 0 for a number
idcs_lim = zeros(N,2);
if ~isempty(names_lim)
    [is_var,idcs_lim] = ismember(names_lim,vars);
    if any(~is_var & ~cellfun(@isempty,names_lim),'all')
        error("Limits of integration must be numbers or variables of integration, named as in 'vars'.")
    end
end

% KCop is a functional, no partial integral: each copy of the output once
% ('si' and 'si_1' are the same copy). cp(k): the copy that variable k is
[is_out,loc] = ismember(vars,names_cp);
cp = zeros(1,N);            cp(is_out) = idcs_cp(loc(is_out));
if ~isequal(sort(cp(is_out)),1:sum(deps_out))
    error("Partial integrals are not supported: variables of integration must include each copy of the output variables of the operator once.")
end

% Order of integration: each variable before the variables in its limits.
% There is none if a variable is in its own limits, directly or through
% the limits of others. The copies come last (their limits are copies)
is_dep = idcs_lim(:,1)==(1:N) | idcs_lim(:,2)==(1:N);      % l a limit of k
ord = zeros(1,0);
while numel(ord)<N
    is_left = true(1,N);        is_left(ord) = false;
    is_next = is_left & ~any(is_dep(is_left,:),1);
    if ~any(is_next)
        error("Limits of integration are circular: a variable appears in its own limits, directly or through the limits of other variables.")
    end
    ord = [ord,find(is_next)];  %#ok<AGROW> at most N variables
end
ord_W = ord(~is_out(ord));

% Intervals of the copies, as the intmat of Kop: for each copy of output
% variable i, the lower limit a or a copy of variable i, the upper limit b
% or a copy of variable i, as 0, d_i+1, or the number of that copy. A
% delta is at a copy only: an 'intop' has no evaluation at a or b
cs = [0,cumsum(deps_out)];
var_cp = repelem(1:numel(deps_out),deps_out);   % output variable of each copy
intmat = arrayfun(@(d) repmat([0,d+1],1,d),deps_out,'UniformOutput',false);
for k=find(is_out)
    i = var_cp(cp(k));
    row = NaN(1,2);
    for j=1:2
        l = idcs_lim(k,j);
        if l>0 && cp(l)>0 && var_cp(cp(l))==i
            row(j) = cp(l)-cs(i);
        elseif l==0 && lims(k,j)==dom_cp(loc(k),j)
            row(j) = (j-1)*(deps_out(i)+1);     % 0 for a, d_i+1 for b
        end
    end
    if any(isnan(row))
        error("Copies of the output variables must be integrated over their full domain, or between or at copies of the same variable (limits of '%s').",vars{k})
    end
    intmat{i}(2*(cp(k)-cs(i))+[-1,0]) = row;
end
% A copy between two copies: the iterated integral equals the integral over
% the region only if the lower limit is below the upper by the other limits
[~,~,is_ord] = intop.canonicalize_intmat(intmat,[],true);
if ~all(cellfun(@all,is_ord))
    error("Limits of the copies of the output variables are not ordered: a copy between two copies requires the lower to be below the upper by the other limits.")
end

% Limits: numbers, or the pvars of the variables
lims = num2cell(lims);
lims(idcs_lim>0) = cellfun(@(name) polynomial({name}),vars(idcs_lim(idcs_lim>0)),'UniformOutput',false);

end
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %



%% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %%
function [vars,lims,intmat,ord_W] = intmat2lims(intmat_in,vars_out,dom_out,deps_out)
% [VARS,LIMS,INTMAT,ORD_W] = INTMAT2LIMS(INTMAT_IN,VARS_OUT,DOM_OUT,
% DEPS_OUT) reads the region of integration of 'int' given as an intmat,
% with the outputs of dom2intmat: the intmat of Kop for the copies in the
% output of Cop, and the limits of the other copies, over which W is
% integrated first.
%
% INPUTS
% - intmat_in:  1 x M cell, for each output variable si a 1 x 2*Di array,
%               the region of the copies si_1,...,si_Di as a row of the
%               intmat of an 'intop'; the first deps_out(i) copies are
%               those of the output of Cop;
% - vars_out:   1 x M 'cellstr', the output variables of Cop;
% - dom_out:    M x 2 array, the domains of the output variables;
% - deps_out:   1 x M array, the number of copies of each output variable
%               in the output of Cop;
%
% OUTPUTS
% - vars:       1 x N 'cellstr', the names si_j of all copies, N=sum(Di);
% - lims:       N x 2 cell, the limits of the copies over which W is
%               integrated: a 'double' (ai or bi), or the 'polynomial'
%               pvar of a copy (empty for the copies of Cop);
% - intmat:     1 x M cell, the region of the copies of Cop, the intmat of
%               Kop;
% - ord_W:      1 x K array, the indices in VARS of the other copies, in
%               the order in which W is integrated over them: each before
%               the copies in its limits;

M = numel(vars_out);
if ~iscell(intmat_in) || numel(intmat_in)~=M
    error("The region must be a 1 x M cell, with a row for each of the M output variables of the operator.")
end
vars = cell(1,0);       lims = cell(0,2);
intmat = cell(1,M);     ord_W = zeros(1,0);
for i=1:M
    row = intmat_in{i};
    d = deps_out(i);        D = numel(row)/2;
    if size(row,1)>1 || D<d || D~=round(D) || any(row~=round(row) | row<0 | row>D+1)
        error("The region of the copies of '%s' must be a 1 x 2*D array of integers 0,...,D+1, with D>=%d.",vars_out{i},d)
    end
    L = row(1:2:end);       U = row(2:2:end);
    % The copies of Cop: from ai or a copy of Cop to bi or a copy of Cop,
    % as for VARS, DOM (an 'intop' has no evaluation at ai or bi)
    LU = [L(1:d),U(1:d)];
    if any(LU>d & LU<=D) || any(L(1:d)==D+1) || any(U(1:d)==0)
        error("Copies of the output variables must be integrated over their full domain, or between or at copies of the same variable (copies of '%s').",vars_out{i})
    end
    row_K = row(1:2*d);
    row_K(row_K==D+1) = d+1;
    intmat{i} = row_K;
    % Order: each copy before the copies in its limits, so first those in
    % the limits of no copy left. W is integrated over the other copies
    nv = numel(vars);
    left = 1:D;
    while ~isempty(left)
        is_next = ~ismember(left,[L(left),U(left)]);
        if ~any(is_next)
            error("Limits of integration are circular: a variable appears in its own limits, directly or through the limits of other variables.")
        end
        ord_W = [ord_W,nv+left(is_next & left>d)];    %#ok<AGROW> at most N copies
        left = left(~is_next);
    end
    % Limits of the other copies: ai, bi, or the pvar of a copy
    names = get_tensorvars(vars_out{i},D,'full');
    lim_vals = [{dom_out(i,1)},cell(1,D),{dom_out(i,2)}];
    for j=unique([L(d+1:D),U(d+1:D)])
        if j>=1 && j<=D
            lim_vals{j+1} = polynomial(names(j));
        end
    end
    lims_i = cell(D,2);
    lims_i(d+1:D,:) = [lim_vals(L(d+1:D)+1).',lim_vals(U(d+1:D)+1).'];
    vars = [vars,names];        lims = [lims;lims_i];     %#ok<AGROW> M variables
end
% A copy between two copies: the iterated integral equals the integral over
% the region only if the lower limit is below the upper by the other limits
[~,~,is_ord] = intop.canonicalize_intmat(intmat,[],true);
if ~all(cellfun(@all,is_ord))
    error("Limits of the copies of the output variables are not ordered: a copy between two copies requires the lower to be below the upper by the other limits.")
end

end
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %