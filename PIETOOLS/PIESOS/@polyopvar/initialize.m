function f = initialize(f)
% F = INITIALIZE(F) checks whether all fields of the 'polyopvar' object F
% are properly specified (and do not conflict), and fills in any missing
% information.
%
% INPUTS
% - f:      m x 1 'polyopvar' object, with field f.C a 'tensopmat' object
%           and f.Z an nZ x d array of monomials, with nZ equal to the
%           number of columns of f.C.ops. Fields f.dims, f.vars, f.dom and
%           f.depmat may be left empty (in/out), in which case they are
%           inferred from f.C; otherwise they are checked against f.C.
%           The input variables f.vars.in may differ from those of f.C: a
%           variable of f.C that f lacks must not be acted on, and a
%           variable of f that f.C lacks must not appear in any monomial
%           (every monomial is acted on by f.C);
%
% OUTPUTS
% - f:      'polyopvar' object representing the same distributed
%           polynomial as the input, but with all fields filled in and
%           f.chkval = true, and with the nonzeros of each row of f.Z
%           moved to the leading columns (order kept), unless the operators
%           acting on the monomial have more factors than its degree: then
%           the zeros mark the factors that are multipliers from R, see
%           'polyopvar' and num_factors;
%
% NOTES
% - Only scalar or column-vector-valued 'polyopvar' objects are
%   supported: an error is thrown if f.matdim(2)>1.
%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PIETOOLS - initialize
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

% Decompose the distributed polynomial
C = f.C;
Z = f.Z;
xvarname = f.varname;
fdims = f.dims;
fvars = f.vars;
fdom = f.dom;
fdepmat = f.depmat;

% Make sure the coefficients are properly specified.
if ~isa(C,'tensopmat')
    error("Field 'C' must be specified as an object of type 'tensopmat' (was: %s).",class(C))
end
if ~issparse(Z)
    Z = sparse(Z);
end
if ~C.chkval
    try  C = initialize(C);
    catch
        error("The coefficient operator is not appropriately declared (try running initialize(C)).")
    end
end
Cdims = C.dims;
Cvars = C.vars;
Cdom = C.dom;
Cdepmat = C.depmat;

% Make sure variable names and monomials are properly specified
if ~iscellstr(xvarname) || numel(unique(xvarname))~=numel(xvarname)
    error("Variable names must be specifed as 1 x p 'cellstr' array with unique elements.")
end
p = numel(xvarname);
if size(Z,1)~=size(C.ops,2)
    error("The number of monomials must match the number of columns of operators in the coefficient operator.")
end
[~,~,Zvals] = find(Z);
if any(any(~ismember(Zvals,0:p)))
    error("Variable index in monomials must be specified as integer between 0 and p, for p state variables.")
end
% Left-justify nonzeros of each monomial, keeping factor order: [0,1,2,0],
% [1,0,2,0] -> [1,2,0,0]. Zeros are then padding (no operator factor), so C is
% unaffected. This is only done for monomials acted on by no more factors than
% their degree: otherwise the zeros are the positions of the factors that are
% multipliers from R, see `polyopvar`.
f.C = C;
[~,nops] = num_factors(f);          % number of factors acting on each monomial
deg = full(sum(Z~=0,2));            % degree of each monomial
[Zr,Zc,Zv] = find(Z);
Zr = Zr(:);     Zc = Zc(:);     Zv = Zv(:);     % (:) since find returns rows for row Z
mv = find(nops(Zr)<=deg(Zr));       % nonzeros to move
if ~isempty(mv)
    [Zr_mv,srt] = sort(Zr(mv));     % sort is stable: keeps factor order
    nZr = accumarray(Zr_mv,1,[size(Z,1),1]);
    row0 = cumsum(nZr) - nZr;
    Zc(mv(srt)) = (1:numel(Zr_mv))' - row0(Zr_mv);
end
Z = sparse(Zr,Zc,Zv,size(Z,1),size(Z,2));


% Make sure the other fields are properly specified
if ~isstruct(fdims) || ~isfield(fdims,'in') || ~isfield(fdims,'out')
    error("Field 'dims' must be specified as a struct with fields 'in' and 'out'.");
elseif ~isstruct(fvars) || ~isfield(fvars,'in') || ~isfield(fvars,'out')
    error("Field 'vars' must be specified as a struct with fields 'in' and 'out'.");
elseif ~isstruct(fdom) || ~isfield(fdom,'in') || ~isfield(fdom,'out')
    error("Field 'dom' must be specified as a struct with fields 'in' and 'out'.");
elseif ~isstruct(fdepmat) || ~isfield(fdepmat,'in') || ~isfield(fdepmat,'out')
    error("Field 'depmat' must be specified as a struct with fields 'in' and 'out'.");
end

% Make sure the output function spaces of f match those of C
if isempty(fdims.out)
    fdims.out = Cdims.out;
elseif ~isequal(fdims.out,Cdims.out)
    error("Output dimensions of the polynomial should match those of the coefficient operator.")
end
if isempty(fvars.out)
    fvars.out = Cvars.out;
elseif ~isequal(fvars.out,Cvars.out)
    error("Output spatial variables of the polynomial should match those of the coefficient operator.")
end
if isempty(fdom.out)
    fdom.out = Cdom.out;
elseif ~isequal(fdom.out,Cdom.out)
    error("Output domains of the spatial variables in the polynomial should match those of the coefficient operator.")
end
if isempty(fdepmat.out)
    fdepmat.out = Cdepmat.out;
elseif ~isequal(fdepmat.out,Cdepmat.out)
    error("Degrees of variables in output of polynomial should match those of the coefficient operator.")
end

% For the input function spaces, we have to be careful that p is defined in
% terms of monomials in x, rather than just x itself.
% First, check the input spatial variables
M = numel(fvars.in);
Nc = numel(Cvars.in);
if isempty(fvars.in)
    fvars.in = Cvars.in;
    M = numel(fvars.in);
    idcs2 = 1:M;
else
    if ~iscellstr(fvars.in)
        error("Input variables should be specified as 1 x M 'cellstr' array.")
    elseif numel(unique(fvars.in))~=numel(fvars.in)
        error("Input variable names in field 'vars.in' must be unique.")
    end
    % Variables of f need not coincide with those of C: a variable on only
    % one side is allowed if nothing depends on / acts on it.
    [~,~,idcs2] = unique([Cvars.in, fvars.in(:)'],'stable');
    idcs2 = idcs2(Nc+1:end);
    idcs2(idcs2>Nc) = 0;
end
inC = reshape(idcs2>0,1,[]);    % 1 x M: variable of f is a variable of C
if isempty(fdom.in)
    if ~all(inC)
        error("Domains of input variables that do not appear in the coefficient operator must be specified in field 'dom.in'.")
    end                                                                  
    fdom.in = Cdom.in(idcs2,:);
else
    if size(fdom.in,1)~=M || ~isequal(Cdom.in(idcs2(inC),:),fdom.in(inC,:))
        error("Input domains of the spatial variables in the polynomial should match those of the coefficient operator.")
    end
end
% Next check for monomials of degree 1 to determine the size and variable
% dependency
p = numel(xvarname);
dims_in = nan(p,1);
depmat_in = false(p,M);
depmat_known = false(p,1);
d = size(Z,2);
if d>=1
    for i=1:p
        % Degree 1 in state i: the zeros may be anywhere, e.g. [0,i] for c*x_i
        deg1_idx = find(sum(Z==i,2)==1 & deg==1,1,'first');
        % Guard against state i never appearing as a "pure" degree-1 
        % monomial in Z (e.g. f=x^2 with no separate linear x term).
        if ~isempty(deg1_idx)
            dims_in(i) = Cdims.in(deg1_idx);
            depmat_in(i,inC) = logical(Cdepmat.in(deg1_idx,idcs2(inC)));
            depmat_known(i) = true;
        else
            % Try to establish variable dependence from a pure power of
            % variable i.
            pure_pow_idx = find(all(Z==i | Z==0,2) & any(Z==i,2),1,'first');
            if ~isempty(pure_pow_idx)
                depmat_in(i,inC) = logical(Cdepmat.in(pure_pow_idx,idcs2(inC)));
                depmat_known(i) = true;
            end
        end
    end
end
% Assume scalar-valued state if not known
is_known = ~isnan(dims_in);
dims_in(~is_known) = 1;
if isempty(fdims.in)
    fdims.in = dims_in;
else
    if ~isequal(dims_in(is_known), fdims.in(is_known))
        error("Dependency array of state variables on spatial variables does not match that expected from the coefficient operator.")
    end
end
if isempty(fdepmat.in)
    fdepmat.in = depmat_in;
else
    % Only variables of C can be inferred; others are checked per monomial below
    if size(fdepmat.in,1)~=p || size(fdepmat.in,2)~=M || ...
            ~isequal(depmat_in(depmat_known,inC), fdepmat.in(depmat_known,inC))
        error("Dependency array of state variables on spatial variables does not match that expected from the coefficient operator.")
    end
end

% Finally, go through all monomials, to check if the set dimensions and
% depmats pose any conflicts
for i=1:size(Z,1)
    Zi = Z(i,Z(i,:)~=0);
    fdim_i = prod(fdims.in(Zi));
    Cdim_i = Cdims.in(i);
    if fdim_i~=Cdim_i
        error("Field 'dims.in' conflicts with field 'C.dims.in'.")
    end
    fdeps_i = sum(fdepmat.in(Zi,:),1);
    Cdeps_i = Cdepmat.in(i,:);
    % Every monomial is acted on by C (no empty columns), so it may not
    % depend on a variable C is not defined in.
    if any(fdeps_i(~inC))
        error("Monomial %d depends on spatial variables on which the coefficient operator is not defined.",i)
    end
    % Variables of C that f lacks must not be acted on; the rest must agree.
    if any(Cdeps_i(setdiff(1:Nc,idcs2(inC)))) || any(Cdeps_i(idcs2(inC))~=fdeps_i(inC))
        error("Field 'depmat.in' conflicts with 'C.depmat.in'.")
    end
end

% Monomials acted on by more factors than their degree: factor j acts on state
% Z(i,j), or is a multiplier from R if Z(i,j)=0, see `polyopvar`. Check that this
% matches the factors of the operator, as moving the zeros to the right, or
% placing them elsewhere, changes which factor acts on which state.
for i=find(nops(:)'>deg(:)')
    zc = find(Z(i,:));                  % positions of the states
    if any(zc>nops(i))
        error("Monomial %d has a state in position %d of Z, beyond the %d factors of the operators acting on it.",i,max(zc),nops(i))
    end
    is_state = false(1,nops(i));        is_state(zc) = true;
    for r=1:size(C.ops,1)
        if ~isa(C.ops{r,i},'tensopvar')
            continue
        end
        [nin,hasvar] = leaf_inputs(C.ops{r,i});
        if numel(nin)~=nops(i)
            error("Operators acting on monomial %d must have the same number of factors.",i)
        end
        % A zero in Z must be a multiplier from R, a state must have its size
        bad = ~is_state & (nin~=1 | hasvar);
        bad(zc) = nin(zc)~=reshape(fdims.in(full(Z(i,zc))),1,[]);
        if any(bad)
            error("Factor(s) %s of the operators acting on monomial %d do not match Z(%d,:): zeros must be multipliers from R (input dimension 1, no input variables), states must have the dimension of the state.",strjoin(string(find(bad)),","),i,i)
        end
    end
end

% Set the newly determined field values
f.C = C;
f.Z = Z;
f.varname = xvarname(:);
f.dims = fdims;
f.vars = fvars;
f.dom = fdom;
f.depmat = fdepmat;
f.chkval = true;

% Only column-vector/scalar polyopvar objects are supported for now
if f.matdim(2)>1                                                           
    error("polyopvar:matrix_valued_not_supported", ...                    
        ['Matrix-valued (multi-column) distributed polynomials are not ' ...
         'currently supported -- only scalar or column-vector-valued ' ...
         '''polyopvar'' objects (matdim(2)<=1) are. This object would have ' ...
         'matdim = [%d,%d].'],f.matdim(1),f.matdim(2))
end

end


function [nin,hasvar] = leaf_inputs(op)
% Column dimension (nin) and whether it has input variables (hasvar), for each
% factor in the first term of the (nested) 'tensopvar' op
nin = zeros(1,0);       hasvar = false(1,0);
fctrs = op.ops(1,:);
for j=1:numel(fctrs)
    opj = fctrs{j};
    if isa(opj,'tensopvar')
        [nin_j,hasvar_j] = leaf_inputs(opj);
    else
        nin_j = size(opj,2);            hasvar_j = ~isempty(opj.vars.in);
    end
    nin = [nin,nin_j];      hasvar = [hasvar,hasvar_j];
end
end
