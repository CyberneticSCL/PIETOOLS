function [prog,Kop] = piesos_intopvar(prog,matdim,vars,dom,deps,intmat,degs)
% [PROG,KOP] = PIESOS_INTOPVAR(PROG,MATDIM,VARS,DOM,DEPS,INTMAT,DEGS)
% adds a functional variable KOP to the PIESOS optimization program,
% acting on the function space
%   L2^{n}[[a1,b1]^q1 x ... x [aN,bN]^qN]
% where n = MATDIM(2), [ai,bi] = DOM(i,:), and qi = DEPS(i).
%
% INPUTS
% - prog:   'struct' specifying a PIESOS program (see `piesos_program`);
% - matdim: 1 x 2 array specifying the matrix dimensions of the functional;
% - vars:   1 x N cellstr specifying the names of the spatial variables in
%           terms of which the functional is defined. Can also be specified
%           as struct with fields 'in' and 'out', where field 'in' provides
%           the spatial variable names, and 'out' must be empty;
% - dom:    N x 2 array specifying for each of the spatial variables the
%           interval that it is defined on. Can also be specified as struct
%           with fields 'in' and 'out', with 'in' providing the spatial
%           domains, and 'out' must be empty. A single interval (1 x 2
%           array) is used for all variables;
% - deps:   1 x N array, specifying for each spatial variable the order to
%           which the input to the functional depends on this variable.
%           A scalar is used for all variables. A variable with deps(i)=0
%           is not included in the functional;
% - intmat: 1 x N cell (option 1) or char (option 2) specifying the
%           integrals to include in the functional:
%           OPTION 1: 1 x N cell, with element i an ntrms(i) x 2*deps(i)
%             array specifying the intervals of integration along
%             variable i to include in the functional variable, in the
%             'intop' format;
%           OPTION 2: a char indicating what kind of integrals to include
%             in the functional. An associated intmat in the 'intop'
%             format is generated accordingly. Must be one of
%               'full':   include all admissible integrals;
%               'nomult': exclude any delta functions in the kernel, so
%                         that the functional is defined by only proper
%                         integrals;
%               'lower':  include only the integral imposing the ordering
%                           si_1 >= si_2 >= ... >= si_deps(i)
%                         for each i in {1,...,N};
%               'upper':  include only the integral imposing the ordering
%                           si_1 <= si_2 <= ... <= si_deps(i)
%                         for each i in {1,...,N}.
%             Here, 'full' includes one integral for each possible
%             ordering of the deps(i) dummy variables of variable i, with
%             ties between variables taken as delta functions
%             (multipliers). 'nomult' includes the deps(i)! integrals
%             without ties, with the lower integral first and the upper
%             integral last.
%           The terms along the N variables form an ntrms(1) x ... x
%           ntrms(N) cell Kop.params, with a kernel for each combination
%           of terms;
% - degs:   (optional) struct with a single field, 'max' or 'full',
%           specifying the degrees of the monomials in the kernels. If
%           degs is specified as a cell or array, it is interpreted as
%           struct('max',degs). Defaults to struct('max',0), for constant
%           kernels.
%           FIELD 'max': the maximal degrees of the monomials, as
%             - a scalar, interpreted as the maximal cumulative degree of
%               the monomials in the dummy variables associated with each
%               of the variables;
%             - a 1 x N cell, with element i specifying the maximal
%               degrees in the dummy variables of variable i, for each of
%               the terms along variable i: a scalar, the maximal
%               cumulative degree, or a 1 x deps(i) array, the maximal
%               degree in each of the deps(i) dummy variables.
%           FIELD 'full': the full arrays of the degrees of the monomials,
%             each row of an array specifying the degrees of one monomial,
%             so that a single row specifies a single monomial. These can
%             be specified as
%             - ntrms(1) x ... x ntrms(N) cell array, with element
%               (i1,...,iN) specifying the monomials in the kernel
%               Each element is either an nmons x sum(deps) array, with a
%               column for each of the sum(deps) dummy variables of the
%               term (those of variable 1 first, in the order of the
%               variables as specified in "vars", and for each variable
%               in the order of its dummy variables), which allows
%               monomials in different variables to be coupled; or a
%               1 x N cell, with element i an nmons(i) x deps(i) array of
%               the monomials in the dummy variables of variable i, so
%               that the monomials in the kernel are all products of
%               these;
%             - 1 x N cell, with element i an nmons(i) x deps(i) array of
%               the monomials in the dummy variables of variable i, used
%               for each of the terms along variable i, as above.
%             In a 1 x N cell, the element for a variable with
%             deps(i)=0, which is discarded, must have no columns, e.g.
%             [].
%           The number of terms ntrms(i) along variable i is the number
%           of rows of intmat{i}; for a generated intmat these are the
%           integrals described above, in the order described above;
%
% OUTPUTS
% - prog:   'struct' specifying the same PIESOS program as the input, but
%           with the functional Kop added;
% - Kop:    'intop' object representing a functional acting on the
%           prescribed function space and defined in terms of the
%           specified integrals. The kernels are 'dpvar' objects, defined
%           by monomials of degrees determined by DEGS, in terms of the
%           dummy variables 'si_dum1',...,'si_dum<deps(i)>'. The
%           variables are sorted as in any 'intop' object;
%
% NOTES
% - A kernel in a term where a dummy variable si_j is set equal to
%   another one, i.e. intmat{i}(t,2*j-1)=intmat{i}(t,2*j), cannot
%   depend on si_j: monomials with nonzero degree in si_j are
%   discarded, and the kernel is zero if no monomials remain.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PIETOOLS - piesos_intopvar
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
% DJ, 10/15/2026: Initial coding;


% Check the matrix dimensions
if ~isnumeric(matdim) || numel(matdim)~=2 || any(abs(round(matdim))~=matdim)
    error("Matrix dimensions of the functional must be specified as 1 x 2 array of nonnegative integers.")
end
matdim = reshape(matdim,1,2);

% Check the input variables
if iscellstr(vars)
    varnames = reshape(vars,1,[]);
elseif isa(vars,'struct')
    if ~isfield(vars,'in')
        error("No input variables have been specified.")
    elseif ~iscellstr(vars.in)
        error("Input variables must be specified as 1 x N cellstr.")
    elseif isfield(vars,'out') && ~isempty(vars.out)
        error("'intop' objects with output variables are not supported.")
    end
    varnames = reshape(vars.in,1,[]);
else
    error("Variables must be specified as 1 x N cellstr, or as struct with field 'in'.")
end
N = numel(varnames);
if numel(unique(varnames))~=N
    error("Variable names must be distinct.")
end

% Check the domains
if isa(dom,'struct')
    if ~isfield(dom,'in')
        error("No domains have been specified for the input variables.")
    elseif isfield(dom,'out') && ~isempty(dom.out)
        error("No output variable domains can be specified for 'intop' objects.")
    end
    dom = dom.in;
end
if ~isnumeric(dom)
    error("Domains should be specified as N x 2 array for N spatial variables.")
elseif N==0
    domval = zeros(0,2);
elseif size(dom,2)~=2
    error("An upper and lower boundary of the domain must be speciified.")
elseif size(dom,1)==1
    % Assume same interval for all variables.
    domval = repmat(dom,N,1);
elseif size(dom,1)~=N
    error("Number of domains must match number of spatial variables.")
else
    domval = dom;
end
if any(domval(:,1)>=domval(:,2))
    error("Lower boundary of a domain must be smaller than the upper boundary.")
end

% Check the degree to which the input depends on each variable
if ~isnumeric(deps) || any(deps(:)<0) || any(round(deps(:))~=deps(:))
    error("Degree of dependence on each spatial variable must be specified as 1 x N array of nonnegative integers for N variables.")
elseif isscalar(deps)
    % Assume same degree in all variables
    deps = deps*ones(1,N);
elseif numel(deps)~=N
    error("Degree of dependence on each spatial variable must be specified as 1 x N array of nonnegative integers for N variables.")
end
deps = reshape(deps,1,N);

% Determine the intervals of integration
if isstring(intmat) && isscalar(intmat)
    % The user wants the function to generate an intmat
    intmat = char(intmat);
end
ntrms = ones(1,N);                  % number of terms along each variable
if isa(intmat,'cell')
    % Assume intervals are already specified
    if numel(intmat)~=N
        error("The input 'intmat' must be specified as 1 x N cell for N spatial variables.")
    end
    for i=1:N
        if ~isnumeric(intmat{i}) || size(intmat{i},2)~=2*deps(i)
            error("Element i of 'intmat' must contain 2*q columns for q the number of dummy variables associated with variable i.")
        end
        if deps(i)>0
            ntrms(i) = size(intmat{i},1);
            if ntrms(i)==0
                error("At least one interval of integration must be specified along each variable.")
            end
        end
    end
elseif ischar(intmat)
    % Generate intmat based on specified option
    if ~strcmp(intmat,'full') && ~strcmp(intmat,'nomult') && ~strcmp(intmat,'lower') && ~strcmp(intmat,'upper')
        error("Option for 'intmat' must be specified as one of 'full', 'nomult', 'lower', or 'upper'.")
    end
    intmat = build_intmat(intmat,deps);
    for i=1:N
        ntrms(i) = size(intmat{i},1);
    end
else
    error("The input 'intmat' must be specified as 1 x N cell or as char.")
end

% Check the format of the degrees: 
%  - a struct with field 'max' (maximal degrees: a scalar or a 1 x N cell)
%  - or 'full' (the monomials: an ntrms(1) x ... x ntrms(N) cell, 
%               or a 1 x N cell of the monomials of each variable).
if nargin<=6
    degs = struct('max',0);
elseif ~isstruct(degs)
    degs = struct('max',{degs});
end
if ~isscalar(degs)
    error("Degrees must be specified as a scalar struct..") 
elseif numel(fieldnames(degs))~=1 || ~(isfield(degs,'max') || isfield(degs,'full')) 
    error("Degrees must be specified as struct with a single field 'max' or 'full'.")
end
per_term = false;           % true if the monomials are specified for each term 
is_max = isfield(degs,'max');
if is_max
    degv = degs.max;
    if isnumeric(degv) && isscalar(degv)
        degs = repmat({degv},1,N);
    elseif iscell(degv) && numel(degv)==N
        degs = reshape(degv,1,N);
    else
        error("Maximal degrees must be specified as a scalar, or as 1 x N cell, with an element for each variable.") 
    end
else
    degv = degs.full;
    if ~iscell(degv)
        error("Full degrees must be specified as ntrms(1) x ... x ntrms(N) cell, with the monomials of each term, or as 1 x N cell, with the monomials of each variable.")
    end                                                                     
    % 1 x N cell of arrays with deps(i) columns for variable i
    is_var_cell = isequal(size(degv),[1,N]) && all(cellfun(@isnumeric,degv(:)));
    for i=1:N
        is_var_cell = is_var_cell && size(degv{i},2)==deps(i);
    end
    Nd = max(N,2);
    is_term_cell = ndims(degv)<=Nd && (isequal(size(degv,1:Nd),[ntrms,ones(1,Nd-N)]) || (N==1 && isvector(degv) && numel(degv)==ntrms(1)));
    if is_var_cell
        degs = reshape(degv,1,N);
    elseif is_term_cell
        per_term = true;
        if any(cellfun(@(a) ~(isnumeric(a) || iscell(a)),degv(:)))
            error("The monomials of a term must be specified as array, or as 1 x N cell.")
        elseif any(cellfun(@(a) iscell(a) && numel(a)~=N,degv(:))) 
            error("The monomials of a term specified as cell must be specified as 1 x N cell, with an element for each variable.")
        elseif any(cellfun(@(a) iscell(a) && ~all(cellfun(@(e,d) isnumeric(e) && size(e,2)==d,reshape(a,1,[]),num2cell(deps))),degv(:)))
            error("The monomials of each variable in a term must be specified as array with deps(i) columns for variable i, so without columns for a variable with deps(i)=0.")
        end
        degs = reshape(degv,[ntrms,1]);
    else
        error("Full degrees must be specified as ntrms(1) x ... x ntrms(N) cell, with an element for each term in the functional, or as 1 x N cell of arrays with deps(i) columns for variable i.")
    end
end

% Discard variables that the functional does not depend on,
% and sort the rest.
is_dep = deps>0;
varnames = reshape(varnames(is_dep),1,[]);    domval = domval(is_dep,:);
ntrms = ntrms(is_dep);          deps = deps(is_dep);
intmat = intmat(is_dep);
if per_term
    degs = cellfun(@(a) select_degs(a,is_dep),degs,'UniformOutput',false);
    degs = reshape(degs,[ntrms,1]);     % drop singleton terms of discarded variables
else
    degs = degs(is_dep);
end
N = numel(varnames);
deps_user = deps;
[varnames,sort_idx] = sort(varnames);
domval = domval(sort_idx,:);
ntrms = ntrms(sort_idx);        deps = deps(sort_idx);
intmat = intmat(sort_idx);
if per_term
    degs = cellfun(@(a) sort_degs(a,sort_idx,deps_user),degs,'UniformOutput',false);
    if N>=2
        degs = permute(degs,[sort_idx,N+1]);
    end
else
    degs = degs(sort_idx);
end

% Check that the intervals of integration are valid, with placeholder
% kernels, before declaring any decision variables, and take the names of
% the dummy variables from 'intop' itself
Kop = intop(repmat({0},[ntrms,1]),intmat,varnames,domval);
varnames_dum = Kop.varnames.in_dum;

% Declare the monomials in terms of which the kernels will be defined
Zcell = build_monoms(degs,per_term,is_max,intmat,deps,varnames_dum,ntrms);

% Finally, build the parameters: a matrix of polynomials in the monomials 
% Zcell, with a decision variable for each coefficient
params = cell([ntrms,1]);
for k = 1:numel(params)
    if isempty(Zcell{k})
        % No monomials left, so the kernel is zero
        params{k} = polynomial(spalloc(1,prod(matdim),0),zeros(1,0),cell(0,1),matdim);
    else
        [prog,params{k}] = sospolymatrixvar(prog,Zcell{k},matdim);
    end
end

% Declare the functional, checking the kernels as well (chkval=false)
Kop = intop(params,intmat,varnames,domval);


end
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %



%% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %%
function intmat = build_intmat(opt,deps)
% INTMAT = BUILD_INTMAT(OPT,DEPS) builds an 'intop' intmat structure based
% on the specified option, with no integrals for the variables with 
% DEPS(i)=0.
%
% NOTES:
% Each integral is defined by the ordering of the dummy variables
% si_1,...,si_d, as a vector r of ranks, with r(j) the rank of si_j in the
% increasing order of the values; variables with the same rank are equal.
% The variables of rank q are all set equal to the first of them, rep(q),
% as the kernel cannot depend on the others. Variable rep(q) is integrated
% over [a,si_rep(q+1)] if q<max(r), and over the full domain if q=max(r).
% Then the variables are ordered as r, and every ordering gives a different
% integral (intmat entries 0 and d+1 stand for the boundaries a and b of
% the domain).

N = numel(deps);
intmat = cell(1,N);
for i=1:N
    d = deps(i);
    if d==0
        intmat{i} = zeros(1,0);
        continue
    end
    switch opt
        case 'lower'
            R = d:-1:1;
        case 'upper'
            R = 1:d;
        otherwise
            R = weak_orders(d);
            if strcmp(opt,'nomult')
                R = R(max(R,[],2)==d,:);
            end
            % Fewest distinct values first, then decreasing ranks, so that
            % for d=2 the order is multiplier, lower, upper (as in monopmat)
            [~,idcs] = sortrows([max(R,[],2),-R]);
            R = R(idcs,:);
    end
    intmat_i = zeros(size(R,1),2*d);
    for t=1:size(R,1)
        r = R(t,:);
        nb = max(r);
        rep = zeros(1,nb);
        for b=1:nb
            rep(b) = find(r==b,1,'first');
        end
        for j=1:d
            b = r(j);
            if j~=rep(b)
                intmat_i(t,2*j-1:2*j) = [rep(b),rep(b)];
            elseif b==nb
                intmat_i(t,2*j-1:2*j) = [0,d+1];
            else
                intmat_i(t,2*j-1:2*j) = [0,rep(b+1)];
            end
        end
    end
    intmat{i} = intmat_i;
end

end
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %



%% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %%
function R = weak_orders(d)
% R = WEAK_ORDERS(D) returns all rank vectors (rows) of the weak orderings
% of D variables, so that R(t,j) is the rank of variable j in the t-th
% ordering, and R(t,:) contains all integers from 1 to max(R(t,:)).
% Built by choosing the variables with rank 1, and ordering the rest.

if d==0
    R = zeros(1,0);
    return
end
Rcell = cell(2^d-1,1);
for s=1:2^d-1
    is_low = logical(bitget(s,1:d));        % variables of rank 1
    nrest = d-nnz(is_low);
    if nrest==0
        Rcell{s} = ones(1,d);
    else
        Rrest = weak_orders(nrest);
        Rs = ones(size(Rrest,1),d);
        Rs(:,~is_low) = Rrest+1;
        Rcell{s} = Rs;
    end
end
R = vertcat(Rcell{:});

end
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %



%% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %%
function Zcell = build_monoms(degs,per_term,is_max,intmat,deps,varnames_dum,ntrms) 
% ZCELL = BUILD_MONOMS(DEGS,PER_TERM,IS_MAX,INTMAT,DEPS,VARNAMES_DUM,NTRMS)
% builds for each term in the functional the column vector of monomials in
% terms of which the associated kernel is defined, as a 'polynomial' object
% in the dummy variables, or [] if there are no such monomials.
% If DEGS is a 1 x N cell, with the maximal degrees (IS_MAX) or monomials 
% of each variable, or if it is a cell of terms with such a 1 x N cell for 
% each term, the monomials are the product of the monomials in the dummy 
% variables along each variable, with the first variable changing slowest. 
% If DEGS is a cell of terms with an array for a term, the array specifies 
% the monomials in all dummy variables of the term together.                               

N = numel(deps);
Zcell = cell([ntrms,1]);
sz = [ntrms,1,1];       % ind2sub needs at least 2 dimensions, also for N=0
idx = cell(1,numel(sz));
for k=1:numel(Zcell)
    [idx{:}] = ind2sub(sz,k);
    % Dummy variables equal to another one do not appear in the kernel
    is_sub = cell(1,N);
    names = cell(1,0);
    for i=1:N
        row = intmat{i}(idx{i},:);
        is_sub{i} = row(1:2:end)==row(2:2:end);
        names = [names, varnames_dum{i}(~is_sub{i})];
    end
    if per_term && isnumeric(degs{k})
        % Monomials in all dummy variables of the term together
        degmat = degs2degmat(degs{k},sum(deps),[is_sub{:}],false);
    else
        degmat = zeros(1,0);
        for i=1:N
            if ~per_term
                spec = degs{i};         % maximal degrees or monomials of variable i 
                max_degs = is_max;
            else
                spec = degs{k}{i};      % monomials of variable i in term k 
                max_degs = false;
            end
            Ei = degs2degmat(spec,deps(i),is_sub{i},max_degs);
            degmat = [repelem(degmat,size(Ei,1),1), repmat(Ei,size(degmat,1),1)];
        end
    end
    nZ = size(degmat,1);
    if nZ==0
        continue
    elseif isempty(names)
        Zcell{k} = polynomial(1);
    else
        Zcell{k} = polynomial(speye(nZ),degmat,names,[nZ,1]);
    end
end

end
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %



%% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %%
function a = select_degs(a,idx)
% A = SELECT_DEGS(A,IDX) restricts the degrees A of a term to the variables
% IDX, if specified as 1 x N cell. Arrays have a column for each dummy
% variable, and the variables that are discarded have none.

if iscell(a)                                                     
    a = a(idx);                                            
end                                                 
end
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %



%% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %%
function a = sort_degs(a,sort_idx,deps)
% A = SORT_DEGS(A,SORT_IDX,DEPS) reorders the monomials A of a term to 
% match the variables sorted as SORT_IDX, where DEPS is the degree of each 
% of the variables in the order specified. The columns of an array, one for 
% each of the dummy variables, are reordered in blocks of DEPS(i) columns.

if iscell(a)
    a = a(sort_idx);
else
    if size(a,2)~=sum(deps)
        error("The degrees of a term specified as array must have a column for each of the %d dummy variables.",sum(deps))
    end
    if numel(deps)>1
        cols = mat2cell(1:sum(deps),1,deps);
        a = a(:,[cols{sort_idx}]);
    end  
end

end
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %



%% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% 
function degmat = degs2degmat(degs,nvars,rm_vars,is_max)
% DEGMAT = DEGS2DEGMAT(DEGS,NVARS,RM_VARS,IS_MAX) computes a nmons x NVARS
% array DEGMAT of monomial degrees based on the degrees specified by DEGS.
% If IS_MAX is true, DEGS specifies maximal degrees: a scalar, the maximal
% cumulative degree, or a 1 x NVARS array, the maximal degree in each
% variable. Otherwise, each row of DEGS specifies the degrees of a monomial,
% with a column for each of the NVARS dummy variables, so that a single row
% specifies a single monomial. NVARS is the number of dummy variables of a
% variable, or of all variables of a term.
% Monomials which are not constant in variables RM_VARS are removed, and
% the columns associated with variables RM_VARS are removed from DEGMAT.

if ~isnumeric(degs) || isempty(degs) || any(degs(:)<0) || any(round(degs(:))~=degs(:))
    error("Degrees of the monomials must be specified as nonnegative integers.")
end
if is_max
    if isscalar(degs)
        % Maximal cumulative degree
        degmat = box_degs(degs*ones(1,nvars));
        degmat = degmat(sum(degmat,2)<=degs,:);
    elseif size(degs,1)==1 && size(degs,2)==nvars 
        % Maximal degree in each variable
        degmat = box_degs(degs);
    else
        error("Maximal degrees of the monomials in d dummy variables must be specified as scalar or 1 x d array. Use the field 'full' to specify monomials.")
    end
elseif size(degs,2)==nvars
    % Degrees of the monomials, in the rows 
    degmat = unique(degs,'rows');
else
    error("Degrees of the monomials in d dummy variables must be specified as nmons x d array, with a column for each dummy variable.")
end
degmat = sortrows(degmat(~any(degmat(:,rm_vars)>0,2),~rm_vars));

end



%% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %%
function degmat = box_degs(maxdegs)
% DEGMAT = BOX_DEGS(MAXDEGS) returns the degrees (rows) of all monomials in
% numel(MAXDEGS) variables, of degree at most MAXDEG(j) in variable j.

d = numel(maxdegs);
vecs = cell(1,d);
for j=1:d
    vecs{j} = 0:maxdegs(j);
end
G = cell(1,d);
[G{:}] = ndgrid(vecs{:});
degmat = zeros(numel(G{1}),d);
for j=1:d
    degmat(:,j) = G{j}(:);
end

end
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
