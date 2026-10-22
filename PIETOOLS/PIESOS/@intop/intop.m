classdef (InferiorClasses={?polynomial,?dpvar,?nopvar,?ndopvar})intop
% INTOP is a class used to represent functionals Kop acting on
% infinite-dimensional states
%   v in L2[[a1,b1]^d1 x ... x [aN,bN]^dN],
% taking the form
%   Kop*v
%   = sum_{i1=1}^{q1}...sum_{iN=1}^{qN} int_{[a1,b1]^d1}...int_{[aN,bN]^dN}
%       II_{i1}(s1)*...*II_{iN}(sN)*K_{i1,...,iN}(s1,...,sN)*v(s1,...,sN)
%                                                           dsN ... ds1
% where the functions II_{j}(si) are indicator functions of the form
%       II_{j}(si) = chi_{[a_{Ai_{j,1}},b_{Bi_{j,1}}]}(si_1)*...
%                       chi_{[a_{Ai_{j,di}},b_{Bi_{j,di}}]}(si_di)
% where Ai and Bi are matrices specifying the domains of integration, and
% where chi is the characteristic function, so that
%                         { 1,          z in [x,y],    x~=y,
%       chi_{[x,y]}(z) := { 0,          z notin [x,y], x~=y,
%                         { delta(x,z), x==y,
% for delta the Dirac delta function.
%
% The kernels K_{i1,...,iN} are all expressed in terms of the same
% monomials Z.in in the dummy variables, and the same decision variables
% Zd, so that the parameters params{i1,...,iN} only hold sparse
% coefficients,
%   K_{i}(s) = kron(I_m,[1,d']) * params{i} * kron(I_n,z(s))
% where d is the vector of decision variables Zd, and z(s) the vector of
% monomials in the dummy variables s defined by Z.in.
%
% CLASS properties
% - params: q1 x q2 x ... x qN cell specifying the kernels defining the
%           functional. In particular, for each i1,...,iN, the element
%           params{i1,...,iN} is a sparse m*(nd+1) x n*nZ array defining
%           the m x n kernel K_{i1,...,iN}, where
%             nd = numel(Zd),   nZ = prod(cellfun(@(a) size(a,1),Z.in)).
%           The layout is that of the coefficients 'C' of a 'dpvar' object
%           with monomials Z: the rows are ordered with the matrix row
%           outer and the decision variable inner (one block of nd+1 rows
%           for each row of the matrix, constant term first), the columns
%           with the matrix column outer and the monomial inner;
% - Z:      struct with fields
%     in:   1 x N cell. Element i is a nmons(i) x deps(i) array of
%           nonnegative integers, each row specifying the degrees of one
%           monomial in the deps(i) dummy variables associated with
%           variable i, {vars.in{i}_dum1,...,vars.in{i}_dum<deps(i)>}. The
%           monomials z(s) in all dummy variables are the Kronecker product
%           of the monomials of the different variables, the first variable
%           changing slowest. The rows of Z.in{i} must be distinct;
%    out:   1 x 0 cell, as the functional has no output variables;
% - Zd:     nd x 1 cellstr specifying the names of the decision variables,
%           shared by all parameters. A parameter that does not depend on
%           all of them has zero rows for those it does not use. The names
%           must be distinct;
% - vars:   struct with fields
%     in:   1 x N cell specifying the names of the spatial variables used
%           to define the kernels. Note that only one variable name is
%           specified along each direction i in {1,...,N}, but each of
%           these variables may in fact itself be multivariate in the
%           kernel. If vars.in{i} = 'si', then the associated variable
%           names in the kernel will be
%               {'si_dum1',...,'si_dumdi'}
%           for di the "degree" of v in variable si;
%    out:   1 x M cell specifying the names of the spatial variables
%           defining the output of the operator. Should be empty (M=0) for
%           the object to be a functional;
% - dom:    struct with fields
%    in:    N x 2 array specifying for each of the spatial variables the
%           interval on which it is defined. Note that if si in [ai,bi],
%           then also si_j in [ai,bi] for each j in {1,...,di};
%   out:    M x 2 array specifying for each of the output variables the
%           interval on which it is defined (should be empty);
% - intmat: 1 x N cell defining the indicator functions used in each term
%           of the functional. In particular, for each k in {1,...,N},
%           intmat{k} is a qk x 2*dk array, with row i corresponding to
%           term i in the functional (along direction k), and columns
%           [2*(j-1)+1,2*j] specifying the interval over which variable
%           sk_j is integrated. Here, if
%           intmat{k}(i,2*(j-1)+1:2*j) = [p,r], then this implies that sk_j
%           is integrated over the domain [sk_p,sk_r], where we set
%           sk_0:=ak and sk_{dk+1}:=bk. If p=r, then the term acts as a
%           multiplier along variable sk_j, evaluating v at sk_j=sk_p.
%           The intervals must be ordered: if both limits of sk_j are
%           variables (p,r in {1,...,dk}, p~=r), sk_p<=sk_r must follow from
%           the limits of the other variables, e.g. because r has lower limit
%           p. Otherwise, chi_{[sk_p,sk_r]}(sk_j) would be nonzero for
%           sk_r<sk_p as well, depending on the order of integration. The
%           constructor rejects rows that are not ordered, e.g.
%           [0,4,1,3,0,4]. Use intop.canonicalize_intmat to split them into
%           ordered rows, [0,2,0,3,0,4], with the same kernel;
% - matdim: 1 x 2 array specifying the matrix dimensions [m,n] of the
%           functional (not stored: determined by params, Zd, and Z);
% - dims:   struct with fields 'in' and 'out', with dims.in = n and
%           dims.out = m the dimensions of the input and output of the
%           functional, as for 'tensopmat' objects (not stored:
%           determined by matdim);
% - deps:   struct with fields 'in' and 'out', with 'in' the 1 x N array of
%           the degrees dk of the variables sk in the functional (not
%           stored: determined by intmat), and 'out' empty;
%
% NOTES
% - vars.in and vars.out always store the spatial variables in the same
%   (sorted) order; this order corresponds positionally to dom.in/dom.out,
%   intmat, Z.in, and params.
% - The kernels may also be specified as a cell of 'double', 'polynomial',
%   or 'dpvar' objects in the constructor, which expresses them in terms of
%   the same monomials and decision variables. The method 'kernels' returns
%   the kernels in this form. The kernels may only depend on the dummy
%   variables of the functional, Kop.varnames.in_dum.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PIETOOLS - intop
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
% DJ, 10/22/2026: Require ordered intmat, see 'intop.canonicalize_intmat';


properties
    params = {sparse(0,0)};
    Z = struct('in',{cell(1,0)},'out',{cell(1,0)});
    Zd = cell(0,1);
    intmat = cell(1,0);
    vars = struct('in',{cell(1,0)},'out',{cell(1,0)});
    dom = struct('in',zeros(0,2),'out',zeros(0,2));
end
properties (Dependent)
    matdim;
    dims;
    deps;
end
properties (Hidden)
    chkval = false;  % indicates whether we have checked whether the object is properly declared
end
properties (Hidden,Dependent)
    int_order;
    varnames;
    pvars;
    depmat;
end


methods

function Kop = intop(params,intmat,vars,dom,Z,Zd,chkval)
% KOP = INTOP(PARAMS,INTMAT,VARS,DOM,Z,ZD,CHKVAL) constructs an 'intop'
% object from a cell PARAMS of sparse coefficient matrices in the monomials
% Z and decision variables ZD, as described above;
%
% KOP = INTOP(KERNELS,INTMAT,VARS,DOM,CHKVAL) constructs an 'intop' object
% from a cell KERNELS of kernels of type 'double', 'polynomial', or
% 'dpvar', expressing all kernels in terms of the union of the monomials
% and decision variables that they involve;
%
% KOP = INTOP(K) constructs a functional without integrals (N=0) from the
% constant K of type 'double', 'polynomial', or 'dpvar'.
%
% INPUTS
% - params: q1 x ... x qN cell, for qk the number of terms along direction
%           k. In the first syntax, element {i1,...,iN} is the sparse
%           m*(nd+1) x n*nZ coefficient matrix of the kernel K_{i1,...,iN}.
%           In the second syntax, it is the m x n kernel K_{i1,...,iN}
%           itself. In the third syntax, it is the single m x n constant K;
% - intmat: 1 x N cell, with element k a qk x 2*dk array of integers
%           specifying the indicator functions along direction k. If N=1,
%           may also be specified as the qk x 2*dk array;
% - vars:   1 x N cellstr or vector of N 'pvar' objects specifying the
%           spatial variables, or struct with fields 'in' and 'out', see
%           above;
% - dom:    N x 2 array of type 'double' specifying the interval of each
%           spatial variable (a 1 x 2 array is used for all variables), or
%           struct with fields 'in' and 'out', see above;
% - Z:      (first syntax) struct with fields 'in' and 'out', or the 1 x N
%           cell Z.in, see above;
% - Zd:     (optional, first syntax) nd x 1 cellstr specifying the names of
%           the decision variables. Defaults to cell(0,1);
% - chkval: (optional) logical. If true, the object is assumed to be
%           properly declared, and this is not checked. Defaults to false;
%
% OUTPUTS
% - Kop:    'intop' object;
%
% NOTES
% - With five inputs, the fifth is read as Z if it is a cell or a struct,
%   and as CHKVAL otherwise.
% - The spatial variables are sorted, and dom, intmat, params, and Z.in
%   are permuted accordingly.
% - Without inputs, INTOP returns an object with the default properties.
% - The rows of intmat must be ordered, see intop.canonicalize_intmat. This
%   is checked by isvalid; only set chkval to true if the intmat is known 
%   to be ordered, e.g. if all rows are those of PI operators.

if nargin==0
    return
elseif nargin==1
    % Declare a functional without integrals
    if isa(params,'double')
        params = {params};
        Kop = intop(params,cell(1,0),cell(1,0),zeros(0,2),true);
        return
    elseif isa(params,'polynomial') || isa(params,'dpvar')
        Kop = intop({params},cell(1,0),cell(1,0),zeros(0,2),true);
        return
    else
        error("For single input case, input must be of type 'double', 'polynomial', or 'dpvar'.")
    end
end
% Distinguish coefficients from kernels
is_coeff = nargin>=6 || (nargin==5 && (iscell(Z) || isstruct(Z)));
if is_coeff
    if nargin<=5
        Zd = cell(0,1);
    end
    if nargin<=6
        chkval = false;
    end
    if iscell(Z)
        Z = struct('in',{reshape(Z,1,[])},'out',{cell(1,0)});
    end
else
    if nargin<=4
        % By default, check the structure
        chkval = false;
    else
        chkval = Z;
    end
end
% Convert polynomial variables to varnames
if isa(vars,'polynomial')
    if ~ispvar(vars)
        error("Variables must be specified as 1 x N 'cellstr' array or vector of pvar objects.")
    end
    pvars = vars;
    vars = struct();
    vars.in = reshape(pvar2varname(pvars),1,[]);
    vars.out = cell(1,0);
elseif iscellstr(vars)
    invars = vars(:)';
    vars = struct();
    vars.in = invars;
    vars.out = cell(1,0);
end
% Convert matrix 'intmat' to cell if appropriate
if isa(intmat,'double') && isscalar(vars.in)
    % The functional acts only on one variable
    intmat = {intmat};
end
% Assume same domain for all variables if not specified
if isa(dom,'double')
    if size(dom,1)==1 && numel(vars.in)>1
        dom = repmat(dom,[numel(vars.in),1]);
    end
    domval = dom;
    dom = struct();
    dom.in = domval;
    dom.out = zeros(0,2);
end
% Store empty variables and intervals as 1x0 cells, as isvalid requires 1 x N
% (a 0x0 cell would not compare equal to the 1x0 cells of intop(c)); only
% the empty case is changed
if iscell(vars.in) && isempty(vars.in),     vars.in = cell(1,0);    end
if iscell(vars.out) && isempty(vars.out),   vars.out = cell(1,0);   end
if iscell(intmat) && isempty(intmat),       intmat = cell(1,0);     end
% Canonicalize the order of the spatial variables
N = numel(vars.in);
if N>=2
    [vars_in_sorted,sort_idx] = sort(vars.in);
    if ~isequal(sort_idx,1:N)
        % Permute domains, intmats, and parameters to match new order
        vars.in = vars_in_sorted;
        dom.in = dom.in(sort_idx,:);
        intmat = intmat(sort_idx);
        params = reshape(params,[size(params,1:N),1]);
        params = permute(params,[sort_idx,N+1]);
        if is_coeff
            % The monomials of the variables are permuted as well
            Zin = Z.in(sort_idx);
            [~,~,nZ] = Zdims(Z.in);
            zperm = Zpermmap(Z.in,sort_idx);
            for k=1:numel(params)
                params{k} = remap_cols(params{k},zperm,nZ);
            end
            Z.in = Zin;
        end
    end
end
% Express the kernels in terms of the same monomials and decision variables
if ~is_coeff
    [params,Z,Zd] = kernels2params(params,intmat,vars);
end
if isempty(Z.in),   Z.in = cell(1,0);   end
% Declare the fields
Kop.params = params;
Kop.Z = Z;
Kop.Zd = Zd(:);
Kop.intmat = intmat;
Kop.vars = vars;
Kop.dom = dom;

if ~chkval
    % Check if the functional is declared properly.
    [logval,errmsg] = isvalid(Kop);
    if ~logval
        error(errmsg)
    end
end
% Indicate that we already checked the function.
Kop.chkval = true;

end



%%
function matdim = get.matdim(Kop)
% MATDIM = GET.MATDIM(KOP) returns the matrix dimensions of the
% functional

% The dimensions follow from the size of the parameters, the number of
% decision variables, and the number of monomials
nd = numel(Kop.Zd);
nZ = prod(cellfun(@(a) size(a,1),Kop.Z.in));
sz = size(Kop.params{1});
matdim = [sz(1)/(nd+1),sz(2)/nZ];
end



%%
function dims = get.dims(Kop)
% DIMS = GET.DIMS(KOP) returns the dimensions of the input and output of
% the functional
%
% OUTPUTS
% - dims:   struct with fields 'in' and 'out', with dims.in = n and
%           dims.out = m, for [m,n] = Kop.matdim.

matdim = Kop.matdim;
dims = struct('in',matdim(2),'out',matdim(1));
end



%%
function depmat = get.depmat(Kop)
% DEPMAT = GET.DEPMAT(KOP) returns the degrees of the different
% variables as they appear in the functional
%
% OUTPUTS
% - depmat: struct with fields 'in' and 'out'. Each field is a 1 x N array,
%           specifying the degree of the input and output to the functional
%           in each of the N spatial variables.

% For each degree of a variable, an interval is specified by two
% columns in Kop.intmat
N = numel(Kop.intmat);
depmat2 = zeros(1,N);
for j=1:N
    depmat2(j) = size(Kop.intmat{j},2)./2;
end
% Set the depmat structure
depmat = struct();
depmat.out = zeros(1,0);    % functional maps to scalar
depmat.in = depmat2;

end



%%
function deps = get.deps(Kop)
% DEPS = GET.DEPS(KOP) returns the degrees of the different
% variables as they appear in the functional, as KOP.DEPMAT
%
% OUTPUTS
% - deps:   struct with fields 'in' and 'out'. Each field is a 1 x N array,
%           specifying the degree of the input and output to the functional
%           in each of the N spatial variables.

deps = Kop.depmat;

end



%%
function pvars = get.pvars(Kop)
% PVARS = GET.PVARS(KOP) returns the variables in terms of which the
% kernels defining Kop are defined.
%
% OUTPUTS
% - pvars:  struct with fields 'in' and 'out'. Each element is a 1 x N
%           cell, with each element a d1x1 'polynomial' array of 'pvar'
%           objects specifying the variables in terms of which Kop.params
%           are defined;

% Extract the variables and degrees
vars_tmp = Kop.vars;
N = numel(vars_tmp.in);
depmats = Kop.depmat;

% Initialize the structure, the functional has no output variables
pvars = struct();
pvars.in = cell(1,N);   pvars.out = cell(1,0);
% Declare the input variables
for i=1:N
    d_in = depmats.in(i);
    pvars_i_in = cell(d_in,1);
    for j=1:d_in
        pvars_i_in{j} = [vars_tmp.in{i},'_dum',num2str(j)];
    end
    pvars.in{i} = polynomial(pvars_i_in);
end

end



%%
function varnames = get.varnames(Kop)
% VARNAMES = GET.VARNAMES(KOP) returns the names of the variables in terms
% of which the kernels defining Kop are defined.
%
% OUTPUTS
% - varnames:   struct with fields 'in', 'in_dum', 'out', and 'out_dum'.
%               Each element is a 1 x N cell, with each element a 1xdi
%               cellstr array, specifying the names of the variables along
%               direction i in terms of which the functional is defined;

% Extract the variables and degrees
vars_tmp = Kop.vars;
N = numel(vars_tmp.in);
depmats = Kop.depmat;

% Initialize the structure, the functional has no output variables
varnames = struct();
varnames.in = repmat(cell(1,0),1,N);
varnames.in_dum = varnames.in;
for i=1:N
    % Declare the input variables
    d_in = depmats.in(i);
    varnames_i = cell(1,d_in);
    varnames_dum_i = cell(1,d_in);
    for j=1:d_in
        varnames_i{j} = [vars_tmp.in{i},'_',num2str(j)];
        varnames_dum_i{j} = [vars_tmp.in{i},'_dum',num2str(j)];
    end
    varnames.in{i} = varnames_i;
    varnames.in_dum{i} = varnames_dum_i;
end
varnames.out = cell(1,0);
varnames.out_dum = cell(1,0);

end



%%
function int_order = get.int_order(Kop)
% INT_ORDER = GET.INT_ORDER(KOP) returns the order in which integration
% must be performed for each of the variables in each of the terms in the
% functional KOP
%
% OUTPUTS
% - int_order:  1 x N cell, where element i is a qi x di array of integers.
%               Here, int_order{i}(j,k) = q implies that degree k of
%               variable i in term j must be integrated qth.
%
% EXAMPLE
% Consider N=1 and d=3, with
%   intmat{1} = [0,2,0,4,2,4; 0,4,0,1,0,2]
% specifying the integrals
%   int_{[a,b]^3} chi_{[a,s2]}(s1) chi_{[a,b]}(s2) chi_{[s2,b]} ... ds
%    + int_{[a,b]^3} chi_{[a,b]}(s1) chi_{[a,s1]}(s2) chi_{[a,s2]} ... ds
% representing
%   int_{a}^{b} int_{a}^{s2} int_{s2}^{b} ... ds2 ds1 ds3
%    + int_{a}^{b} int_{a}^{s1} int_{a}^{s2} ... ds3 ds2 ds1
% In the first term, variable s2 must be integrated last, as it is
% integrated over the full domain [a,b]. Variables s1 and s3 can be
% integrated in either order. In the second term, variable s3 must be
% integrated first, then s2, and then s1. We get
%   int_order{1} = [1,2,1; 3,2,1];

% Check that the functional is properly specified
if ~Kop.chkval
    [logval,errmsg] = isvalid(Kop);
    if ~logval
        error(errmsg)
    end
end

% Extract the intervals of integration along each spatial direction
intmats = Kop.intmat;
N = numel(intmats);
int_order = cell(1,N);
for i=1:N
    % Consider the intervals for the ith variable
    intmat_i = intmats{i};
    [nterms,di] = size(intmat_i);
    di = di/2;
    %
    var_order_i = zeros(nterms,di);
    for j=1:nterms
        % Split the intervals into lower and upper boundaries
        Aj = intmat_i(j,1:2:end);
        Bj = intmat_i(j,2:2:end);
        % Declare the order of the variables in term j
        lvls_j = zeros(1,di);
        full_int_vars = [0,di+1];   % variables of which a full integral is taken
        var_rem = 1:di;             % variables of which interval of integration is still variable
        lvl = 1;                    % order assigned to variables over which full integral is taken
        while ~isempty(var_rem)
            % Find all variables that are integrated over the full domain
            % at level "lvl"
            is_full_int = ismember(Aj,full_int_vars) & ismember(Bj,full_int_vars);
            new_vars = var_rem(is_full_int);
            lvls_j(new_vars) = lvl;
            % Add the fully integrated variables to the list
            full_int_vars = [full_int_vars,new_vars];
            var_rem = setdiff(var_rem,new_vars);
            Aj = Aj(~is_full_int);      Bj = Bj(~is_full_int);
            % Move to the next level of integrals
            lvl = lvl+1;
            if lvl>di+1
                error("Order of integrals appears to be ambiguous, something is wrong....")
            end
        end
        % The order of integration should be inner integral first
        var_order_i(j,:) = max(lvls_j)-lvls_j+1;
    end
    int_order{i} = var_order_i;
end

end

end

%%
% Static methods, defined in the files of the class folder
methods (Static)
    [intmat,params,is_ord] = canonicalize_intmat(intmat,params,check_only)
end

end
