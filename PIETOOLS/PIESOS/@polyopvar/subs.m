function fnew = subs(f,var,val)
% FNEW = SUBS(F,VAR,VAL) takes a distributed polynomial F, and replaces
% the states, or the factors in the monomial, specified by VAR by the
% distributed polynomial or constant VAL.
%
% INPUTS
% - f:      'polyopvar' object representing a distributed polynomial,
%           taking the form f(x) = Kop*Z(x), e.g. for functional (intop)
%           Kop;
% - var:    2 OPTIONS
%           - 'cellstr' or 'polyopvar' object specifying the names of the
%             variables which to substitute for expression val;
%           - integer value or array with entries between 1 and
%             d=d1+...+dp, specifying which of the factors in the
%             monomial defining f to substitute. The factors are numbered
%             by their position in Z, so that the factors that are
%             multipliers from R (zeros in Z) are counted but cannot be
%             substituted;
% - val:    'polyopvar' object representing the polynomial function,
%               p(x)(s) = (Aop*x)(s) + (Bop*x^{o 2})(s) + ...,
%           with which to replace the specified states or factors, or a
%           constant (see NOTES), or a 'cellstr' of new names if var is a
%           'cellstr';
%
% OUTPUTS
% - fnew:   'polyopvar' object representing the distributed polynomial
%           defined by f, but with the factor "var" replaced by the
%           specified polynomial, e.g. if var=i<d1,
%           fnew(x):=Kop(x1^{o(i-1)} o p(x) o x1^{o(d1-i)} o ...
%                           ... o xp^{o dp});
%
% NOTES
% - f may have several rows, each a functional (intop) or a function
%   (tensopvar).
% - val may also be a constant: a 'double', 'polynomial' or 'dpvar'
%   (without decision variables) column vector, stacking the values of
%   all substituted states/factors, in order. If only some of the factors
%   in a monomial are replaced, e.g. x in x*y, they remain as known
%   factors (multipliers from R, zeros in Z, see 'polyopvar'), and the
%   result is cleaned (clean_poly).
% - If var specifies states, val has a row for each, in the same order,
%   and all factors in these states are replaced, e.g. subs(V,u,Top*x).
%   If val is a cellstr as well, the states are renamed. Several factors
%   in a monomial are then substituted one at a time (if val does not
%   contain the states replaced), whereas they are substituted all at
%   once if var specifies the factors by number.
% - The operator acting on a monomial may be nested, i.e. have factors
%   that are 'tensopvar' objects acting on several factors of the
%   monomial, see tensopvar/plus. Some or all factors can then be
%   substituted.
%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PIETOOLS - subs
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
% DJ, 10/28/2026: Substitute constants in functionals of degree d>1;

if ~isa(f,'polyopvar')
    error("Input must be of type 'polyopvar'.")
end
if ~f.chkval
    f = initialize(f);
end

if isnumeric(val) || isa(val,'polynomial') || isa(val,'dpvar')
    % Substitution by a constant (function) rather than a polynomial
    fnew = subs_const(f,var,val);
elseif iscellstr(var) || isa(var,'polyopvar')
    % Substitution of a variable for a polynomial
    fnew = subs_var(f,var,val);
elseif isnumeric(var)
    % Substitution of a monomial factor for a polynomial
    fnew = subs_factor(f,var,val);
else
    error("Variables or factors to substitute must be expressed as 'polyopvar' or numeric array, respectively.")
end

end
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %



%% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %%
function fnew = subs_var(f,var,val)
% FNEW = SUBS_VAR(F,VAR,VAL) replaces the states VAR of F (cellstr or
% polyopvar of individual states) by the rows of the polynomial VAL, or
% renames them if VAL is a cellstr too, fnew = f|_{var=val}

if iscellstr(var) && iscellstr(val)
    % Replace each of the variable names in 'var' by variable names 'val'
    [~,idcs1,idcs2] = intersect(f.varname,var(:));
    varname_new = f.varname(:);
    varname_new(idcs1) = val(idcs2);

    % Remove duplicate variables, canonicalizing to sorted order
    nvars_old = numel(f.varname);
    [varname_new,~,idcs2] = unique(varname_new);   % varname_unique(idcs2)=varname_new
    f.varname = varname_new;
    % Move the function space of each state to its new position. States renamed to
    % the same name must live on the same space
    dims_new = zeros(numel(varname_new),1);
    depmat_new = zeros(numel(varname_new),size(f.depmat.in,2));
    dims_new(idcs2) = f.dims.in;
    depmat_new(idcs2,:) = f.depmat.in;
    if ~isequal(dims_new(idcs2),f.dims.in(:)) || ~isequal(depmat_new(idcs2,:)~=0,f.depmat.in~=0)
        error("States that are renamed to the same name must be defined on the same function space.")
    end
    f.dims.in = dims_new;
    f.depmat.in = cast(depmat_new,'like',f.depmat.in);
    idcs2 = [0; idcs2(:)];
    [ridcs,cidcs,vals] = find(f.Z);
    f.Z = sparse(ridcs,cidcs,idcs2(vals+1),size(f.Z,1),size(f.Z,2));    % keep the zeros of Z
    if numel(varname_new)<nvars_old
        % Monomials that coincide after merging states must be added
        f = combine_terms(f);
    end
    fnew = f;
    return
elseif isa(var,'polyopvar')
    % Extract the names of the variables to substitute
    [isa_pvar,var_names] = ispvar(var);
    if ~isa_pvar
        error("Variables to substitute should be specified as array of individial states.")
    elseif numel(var_names)~=numel(var.varname)
        error("Variables to substitute should be unique.")
    end
    var = var_names;
end

% Names that are not states of f are ignored, together with their rows of val
is_state = ismember(var(:),f.varname);
if ~all(is_state)
    if isa(val,'polyopvar')
        if numel(val.C.dims.out)~=numel(var)
            error("Number of rows of polynomial expression should match number of variables.")
        end
        val = val_rows(val,find(is_state));
    end
    var = var(is_state);
    if isempty(var)
        fnew = f;
        return
    end
end

% Express the functional and new expression in terms of the same state
% variables
fvarname = f.varname;
[f,val] = common_states(f,val);

% Determine for each of the variables in "var" which variable in f they
% correspond to
[fvarname,fvar_idcs] = intersect(f.varname,fvarname);
[~,var_idcs,idcs2] = intersect(fvarname,var(:));
var_idcs = fvar_idcs(var_idcs);

% Update order of val to match order of var_idcs
if ~isa(val,'polyopvar')
    error("New values of variables must be specified as object of type 'cellstr' or 'polyopvar'.")
end
if numel(val.C.dims.out)~=numel(var)
    error("Number of rows of polynomial expression should match number of variables.")
end
val = val_rows(val,idcs2);

% The function space of val must match that of the states replaced. Compare
% by variable name: the output variables of val may be a subset of those of
% f, if the states depend on some of the variables only
[depmat_val,is_in] = depmat_in_vars(val,f.vars.in);
depmat_val = logical(depmat_val);
if ~all(is_in) || ~isequal(f.dims.in(var_idcs),val.dims.out) ||...
        ~isequal(logical(f.depmat.in(var_idcs,:)),depmat_val)
    error("Output space of expression must match the function space of the variables.")
end

% Loop over each of the monomials in Z
fnew = 0;
Z = f.Z;
C = f.C;
for j=1:size(Z,1)
    % Extract term j of the polynomial
    f_j = f;
    f_j.Z = Z(j,:);
    f_j.C.ops = C.ops(:,j);
    f_j.C.dims.in = C.dims.in(j);
    f_j.C.depmat.in = C.depmat.in(j,:);
    % Determine which factors in the monomial correspond to variables to
    % substitute, including repeated factors
    [is_fctr,val_idcs] = ismember(full(Z(j,:)),var_idcs);
    fctrs = find(is_fctr);
    val_idcs = val_idcs(is_fctr);
    if isempty(fctrs)
        fnew = fnew+f_j;
        continue
    end
    % Substitute the factors for the desired expressions: one at a time if
    % val does not contain the states replaced
    if numel(fctrs)>1 && ~any(ismember(full(val.Z(:)),var_idcs))
        f_j = subs_fctrs(f_j,var_idcs,val);
    else
        f_j = subs_factor(f_j,fctrs,val_rows(val,val_idcs));
    end
    fnew = fnew + f_j;
end

end
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %



%% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %%
function fnew = subs_fctrs(f,var_idcs,val)
% FNEW = SUBS_FCTRS(F,VAR_IDCS,VAL) substitutes all factors in the single
% monomial of F that are the states VAR_IDCS by the corresponding rows of
% VAL, one factor at a time, for VAL not containing the states VAR_IDCS.
%
% subs_factor merges equivalent integrals (combine_terms) after each
% factor. The factors are found by their state, as combine_terms sorts them.

[is_fctr,val_idcs] = ismember(full(f.Z(1,:)),var_idcs);
j = find(is_fctr,1,'first');
if isempty(j)
    fnew = f;
    return
end
fnew = subs_factor(f,j,val_rows(val,val_idcs(j)));
% Substitute the remaining factors in each of the monomials of the result
frest = 0;
for k=1:size(fnew.Z,1)
    if all(cellfun(@isempty,fnew.C.ops(:,k)))
        continue
    end
    fk = fnew;
    fk.Z = fnew.Z(k,:);
    fk.C.ops = fnew.C.ops(:,k);
    fk.C.dims.in = fnew.C.dims.in(k);
    fk.C.depmat.in = fnew.C.depmat.in(k,:);
    frest = frest + subs_fctrs(fk,var_idcs,val);
end
if isa(frest,'polyopvar')
    fnew = frest;
end
end
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %



%% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %%
function fnew = subs_factor(f,fctr,val)
% FNEW = SUBS_FACTOR(F,FCTR,VAL) replaces the factors FCTR (unique positive
% integers, by position in Z, so that factors that are multipliers from R
% (zeros) are counted but cannot be substituted) in the single monomial of F
% by the rows of the polynomial VAL, fnew = f|_{f.Z(fctr)=val}

if ~isa(fctr,'double') || any(round(abs(fctr))~=fctr) || any(fctr<1) || numel(unique(fctr))~=numel(fctr)
    error("Indices of factors to substitute must be specified as unique positive integers.")
end

% Make sure we are only substituting a single monomial
if size(f.Z,1)==0
    fnew = f;
    return
elseif size(f.Z,1)>1
    error("Substitution of individual factor in multiple monomials is not supported.")
end
% The monomial has nf factors: its degree, the number of nonzeros in Z (padded
% with zeros), and the multipliers from R at the other zeros of Z, if the
% operator has more factors
nf = num_factors(f);
zf = full(f.Z(1,:));        zf(end+1:nf) = 0;
if any(fctr>nf)
    error("Index of factor to substitute cannot exceed the number of factors in the monomial, i.e. its degree plus the multipliers from R.")
elseif any(zf(fctr)==0)
    error("Factors to substitute must be states, not multipliers from R (zeros in Z).")
end
% Constants are dealt with in subs_const
if ~isa(val,'polyopvar')
    error("Substitution of variables with non-'polyopvar' objects is currently not supported.")
end
if numel(fctr)~=numel(val.C.dims.out)
    error("Number of rows in the new expression must match number of factors to substitute.")
end

% Express the functional and new value in terms of the same variables
if ~isequal(f.varname,val.varname)
    [f,val] = common_states(f,val);
    zf = full(f.Z(1,:));        zf(end+1:nf) = 0;   % states of f may have been renumbered
end

% Perform substitution of the monomial defining f:
%   f.Z = (I1 o ... o Id)(x1 o ... o xd),   val = val.C*val.Z
%   --> f.Z|_{xi=val} = (I1 o ... o I{i-1} o val.C o ... o Id)
%                           (x1 o ... o x{i-1} o val.Z o ... o xd)
%                     = Rop * Z
% A constant (degree-0) val stays as a factor, a multiplier from R (zero in Z)
Rop = 1;
Zrows = {zeros(1,0)};       % new monomials, each as long as its number of factors
for j=1:nf
    % Check whether factor j is substituted
    is_fctr_j = fctr==j;
    if any(is_fctr_j)
        % Factor j is replaced by val(j), which must live on the space of factor j,
        % compared by variable name as the output variables of val_j may be a
        % subset of those of f
        dim_j = f.dims.in(zf(j));
        depmat_j = logical(f.depmat.in(zf(j),:));
        val_j = val_rows(val,is_fctr_j);
        [Dval_j,is_in_j] = depmat_in_vars(val_j,f.vars.in);
        if val_j.dims.out~=dim_j || ~is_in_j || any(Dval_j~=depmat_j)
            error("Input function space of factor to substitute must match output space of new expression.")
        end
        % Take tensor product of the coefficients
        Rop = otimes(Rop,val_j.C,2);
        % Take tensor product of the monomials, giving each monomial of val_j as many
        % entries as it has factors (zeros for multipliers from R)
        nval = size(val_j.Z,1);
        nfac_val = num_factors(val_j);
        Zval = cell(1,nval);
        for k=1:nval
            zk = full(val_j.Z(k,:));        zk(end+1:nfac_val(k)) = 0;
            Zval{k} = zk(1:nfac_val(k));
        end
        Zrows = reshape(cellfun(@(a,b) [a,b],repelem(Zrows,1,nval),repmat(Zval,1,numel(Zrows)),'UniformOutput',false),1,[]);
    elseif zf(j)==0
        % Factor j is a multiplier from R, retained as is: identity on R
        Rop = otimes(Rop,mat2tensopvar(1),2);
        Zrows = cellfun(@(a) [a,0],Zrows,'UniformOutput',false);
    else
        % Factor j is retained as is: identity on the space of the state f.Z(j)
        dim_j = f.dims.in(zf(j));
        depmat_j = logical(f.depmat.in(zf(j),:));
        vars_in_j = f.vars.in(depmat_j);
        dom_in_j = f.dom.in(depmat_j,:);
        vars = struct('in',{vars_in_j},'out',{vars_in_j});
        dom = struct('in',{dom_in_j},'out',{dom_in_j});
        Iop_j = mat2tensopvar(eye(dim_j),vars,dom);
        % Take tensor product with the coefficients and monomials
        Rop = otimes(Rop,Iop_j,2);
        Zrows = cellfun(@(a) [a,zf(j)],Zrows,'UniformOutput',false);
    end
end
% Pad the monomials with zeros after the last factor. The zeros before it are
% the positions of multipliers from R, and must stay in place
Z = zeros(numel(Zrows),max([0,cellfun(@numel,Zrows)]));
for k=1:numel(Zrows)
    Z(k,1:numel(Zrows{k})) = Zrows{k};
end

% Replace the old monomial with the new polynomial
%   f.C * f.Z|_{xi=val} = f.C * Rop * Z = Kop * Z
fnew = f;
fnew.Z = Z;
fnew.C = f.C*Rop;
fnew = combine_terms(fnew);

end
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %



function val = val_rows(val,idcs)
% Rows idcs of the vector of polynomials val
val.C.ops = val.C.ops(idcs,:);
val.C.dims.out = val.C.dims.out(idcs);
val.C.depmat.out = val.C.depmat.out(idcs,:);
val.dims.out = val.C.dims.out;
val.depmat.out = val.C.depmat.out;
end
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %



function [D,is_in] = depmat_in_vars(val,vars)
% [D,IS_IN] = DEPMAT_IN_VARS(VAL,VARS) returns the dependence of the rows of the
% polyopvar VAL on the variables VARS, by name: D(i,j) is the degree to which
% row i depends on VARS{j}. IS_IN(i) is false if row i depends on a variable
% not in VARS.
D = zeros(size(val.depmat.out,1),numel(vars));
[tf,loc] = ismember(val.vars.out,vars);
D(:,loc(tf)) = val.depmat.out(:,tf);
is_in = ~any(val.depmat.out(:,~tf)~=0,2);
end
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %



%% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %%
function fnew = subs_const(f,var,val)
% FNEW = SUBS_CONST(F,VAR,VAL) replaces the states or factors VAR of F by
% the constant (function) VAL, by applying the coefficient operators of F to
% VAL.
%
% INPUTS
% - var:    'cellstr' or 'polyopvar' array of individual states, or integer
%           array of factors in the monomial defining F (then F must have a
%           single monomial), by position in Z;
% - val:    column vector of type 'double', 'polynomial' (or 'dpvar' without
%           decision variables), stacking the values of the states/factors
%           in VAR, in order. A 'polynomial' may only be expressed in terms
%           of the spatial variables of the corresponding state;
%
% NOTES
% - A monomial in which all states are replaced becomes a constant (degree-0)
%   term, Kop*(val o ... o val). Functions are supported only if their output
%   is a pointwise (not tensor) product of the factors, as this relies on
%   apply_op of 'tensopmat'.
% - In a monomial in which only some states are replaced, the factors replaced
%   stay in the operator, as multipliers from R (zeros in Z, see `polyopvar`)
%   mapping the dummy input 1 to the constant function (in the variables the
%   state depends on), by subs_factor. Afterwards the result is cleaned
%   (clean_poly): consecutive multipliers from R are combined, monomials that
%   coincide are added, and states that no longer appear are removed;

% Constants cannot depend on decision variables
if isa(val,'dpvar')
    if ~isempty(val.dvarname)
        error("Substitution by a constant depending on decision variables is not supported.")
    end
    val = dpvar2poly(val);      % no decision variables: polynomial(dpvar) is not supported
end
if size(val,2)~=1
    error("Constants for substitution must be specified as column vector, stacking the values of all variables or factors to substitute.")
end
nZ = size(f.Z,1);
nfac = num_factors(f);          % number of factors of each monomial

% Determine for each monomial the values of its factors
%   Xcell{j}:   1 x nfac(j) cell with the values of the factors in monomial j (the dummy
%               input 1 for multipliers from R), or empty if not all states in monomial j
%               are substituted
%   Pcell{j}:   if only some of the states in monomial j are substituted, 1 x 2 cell with the
%               positions in Z(j,:) of these factors, and 1 x m cell with their values,
%               otherwise empty
Xcell = cell(nZ,1);
Pcell = cell(nZ,1);
if isa(var,'polyopvar')
    [isa_pvar,var] = ispvar(var);
    if ~isa_pvar
        error("Variables to substitute should be specified as array of individial states.")
    end
end
if iscellstr(var)
    % Substitute states by name
    [is_state,state_idcs] = ismember(var(:),f.varname);
    if ~all(is_state)
        error("Variables to substitute by a constant must be states of the polynomial.")
    elseif numel(unique(state_idcs))~=numel(state_idcs)
        error("Variables to substitute should be unique.")
    end
    dims = f.dims.in(state_idcs);
    if size(val,1)~=sum(dims)
        error("Number of rows of constant should match the sum of the dimensions of the variables to substitute.")
    end
    rows = [0;cumsum(dims(:))];
    vals_state = cell(numel(f.varname),1);
    for k=1:numel(state_idcs)
        vals_state{state_idcs(k)} = val(rows(k)+1:rows(k+1));
    end
    for j=1:nZ
        Zj = full(f.Z(j,f.Z(j,:)~=0));
        is_sub = reshape(~cellfun(@isempty,vals_state(Zj)),1,[]);
        if ~any(is_sub)
            continue
        elseif ~all(is_sub)
            zpos_j = find(f.Z(j,:));
            Pcell{j} = {zpos_j(is_sub), reshape(vals_state(Zj(is_sub)),1,[])};
            continue
        end
        Xcell{j} = num2cell(ones(1,nfac(j)));
        Xcell{j}(find(f.Z(j,:))) = vals_state(Zj);
    end
elseif isnumeric(var)
    % Substitute factors in the monomial by position
    if nZ~=1
        error("Substitution of individual factors is only supported for polynomials defined by a single monomial.")
    end
    zj = full(f.Z(1,:));        zj(end+1:nfac(1)) = 0;      % factors by position in Z
    Zj = zj(zj~=0);
    fctr = var(:)';
    if any(round(abs(fctr))~=fctr) || any(fctr<1) || any(fctr>numel(zj)) || numel(unique(fctr))~=numel(fctr)
        error("Indices of factors to substitute must be unique integers not exceeding the number of factors in the monomial, i.e. its degree plus the multipliers from R.")
    elseif any(zj(fctr)==0)
        error("Factors to substitute must be states, not multipliers from R (zeros in Z).")
    end
    dims = f.dims.in(zj(fctr));
    if size(val,1)~=sum(dims)
        error("Number of rows of constant should match the sum of the dimensions of the factors to substitute.")
    end
    rows = [0;cumsum(dims(:))];
    X = num2cell(ones(1,numel(zj)));    % dummy input 1 for multipliers from R
    for k=1:numel(fctr)
        X{fctr(k)} = val(rows(k)+1:rows(k+1));
    end
    if numel(fctr)==numel(Zj)
        Xcell{1} = X;
    else
        Pcell{1} = {fctr, X(fctr)};
    end
else
    error("Variables or factors to substitute must be expressed as 'cellstr' or 'polyopvar', or numeric array, respectively.")
end

% Nothing to do if no monomials are affected
sub_idcs = find(~cellfun(@isempty,Xcell))';
part_idcs = find(~cellfun(@isempty,Pcell))';
if isempty(sub_idcs) && isempty(part_idcs)
    fnew = f;
    return
end

% Monomials in which only some of the states are replaced: each factor replaced
% stays in the operator, as a multiplier from R mapping the dummy input 1 to the
% constant. One factor at a time, as the positions of the factors in Z do not change
fpart = 0;
for j=part_idcs
    f_j = f;
    f_j.Z = f.Z(j,:);
    f_j.C.ops = f.C.ops(:,j);
    f_j.C.dims.in = f.C.dims.in(j);
    f_j.C.depmat.in = f.C.depmat.in(j,:);
    zj = full(f.Z(j,:));
    for k=1:numel(Pcell{j}{1})
        pos = Pcell{j}{1}(k);       Yk = Pcell{j}{2}{k};
        % A constant function may only depend on the variables of the state
        deps = logical(f.depmat.in(zj(pos),:));
        if isa(Yk,'polynomial') && ~all(ismember(Yk.varname,f.vars.in(deps)))
            error("Polynomials for substitution may only be expressed in terms of the spatial variables of the variable they replace.")
        end
        % The constant function as degree-0 polyopvar, a multiplier from R to the
        % function space of the state (states left empty: determined by 'initialize')
        vars = struct('in',{cell(1,0)},'out',{f.vars.in(deps)});
        dom = struct('in',zeros(0,2),'out',f.dom.in(deps,:));
        g = polyopvar();
        g.varname = cell(0,1);
        g.Z = sparse(1,0);
        g.C = tensopmat(mat2tensopvar(Yk,vars,dom));
        f_j = subs_factor(f_j,pos,initialize(g));
    end
    fpart = fpart + f_j;
end

% Monomials in which all states are replaced
if isempty(sub_idcs)
    g0 = 0;
else
    if any(any(f.depmat.out>1))
        error("Substitution of a constant in a polynomial function whose output is a tensor product of its factors is not supported.")
    end
    % Apply the operators acting on the substituted monomials to the constants,
    %   Kop_j*(val_1 o ... o val_d),
    % for each row i of operators separately, as rows may depend on different
    % output variables
    nr = size(f.C.ops,1);
    Yrows = cell(nr,1);
    for j=sub_idcs
        Zj = full(f.Z(j,f.Z(j,:)~=0));
        % A functional of degree d>1 acts on the tensor product of the d values (apply_op of 'intop')
        zpos = find(f.Z(j,:));          % positions of the states in Z
        for k=1:numel(Zj)
            % A constant function may only depend on the variables of the state
            if isa(Xcell{j}{zpos(k)},'polynomial') && ~all(ismember(Xcell{j}{zpos(k)}.varname,f.vars.in(logical(f.depmat.in(Zj(k),:)))))
                error("Polynomials for substitution may only be expressed in terms of the spatial variables of the variable they replace.")
            end
        end
        for i=1:nr
            % Extract the operator in row i acting on monomial j
            Cij = f.C;
            Cij.ops = f.C.ops(i,j);
            Cij.dims.in = f.C.dims.in(j);
            Cij.depmat.in = f.C.depmat.in(j,:);
            Cij.dims.out = f.C.dims.out(i);
            Cij.depmat.out = f.C.depmat.out(i,:);
            if isempty(Yrows{i})
                Yrows{i} = apply_op(Cij,Xcell(j));
            else
                Yrows{i} = Yrows{i} + apply_op(Cij,Xcell(j));
            end
        end
    end

    % Declare the resulting constant (degree-0) term, mapping the dummy input 1
    % to Y, which is a function of the same output variables as f, row by row
    Gcell = cell(nr,1);
    for i=1:nr
        Y = Yrows{i};
        if isa(Y,'double')
            % Empty operators (zero) give a scalar double; expand to the size of row i
            Y = polynomial(Y.*ones(f.C.dims.out(i),1));
        end
        % apply_op of 'tensopvar' renames the output variables s to s_1 for operators with a
        % Kronecker (otype(1)=2) output; as the output depends at most linearly
        % on each variable (depmat.out<=1), rename them back
        for k=1:numel(f.vars.out)
            Y.varname(strcmp(Y.varname,[f.vars.out{k},'_1'])) = f.vars.out(k);
        end
        % Output variables on which row i depends
        if isempty(f.vars.out)
            is_out_i = false(1,0);
        else
            is_out_i = logical(f.depmat.out(i,:));
        end
        vars = struct('in',{cell(1,0)},'out',{f.vars.out(is_out_i)});
        dom = struct('in',zeros(0,2),'out',f.dom.out(is_out_i,:));
        Gcell{i} = mat2tensopvar(Y,vars,dom);
    end
    g0 = polyopvar();
    g0.Z = sparse(1,0);
    g0.C = tensopmat(Gcell);
    g0.varname = cell(0,1);         % no states: input spaces determined by 'initialize', as in constant2polyopvar
    g0 = initialize(g0);
end

% Add the remaining monomials, which are not affected by the substitution
is_kept = cellfun(@isempty,Xcell) & cellfun(@isempty,Pcell);
fnew = g0 + fpart;
if any(is_kept)
    frest = f;
    frest.Z = f.Z(is_kept,:);
    frest.C.ops = f.C.ops(:,is_kept);
    frest.C.dims.in = f.C.dims.in(is_kept);
    frest.C.depmat.in = f.C.depmat.in(is_kept,:);
    fnew = frest + fnew;
end

% After partial substitution, clean the result, which first combines the
% consecutive factors that are now multipliers from R (combine_factors)
if ~isempty(part_idcs)
    fnew = clean_poly(fnew);
end

end
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
