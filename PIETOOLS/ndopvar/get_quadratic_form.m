function [C,Zd1,Zd2,Zdec] = get_quadratic_form(P,var1,var2,dvars,deg1,deg2)
% [C,ZD1,ZD2,ZDEC] = GET_QUADRATIC_FORM(P,VAR1,VAR2,DVARS,DEG1,DEG2)
% Compute coefficients C such that 
%   P(x,y;z)=(Im o Zd1(x))^T (Ik o [1;dvars])^T C (In o Zd2(y)),
% where x is var1 and y is var2, [m,n] are the dimensions of P,
% o denotes the Kronecker product, and
%   Zd1(x) = Zd1_{1}(x1) o Zd1_{2}(x2) o ... o Zd1_{p}(xp)
% for Zd1_{i}(xi) the vector of monomials of degree at most deg1(i) in xi;
% INPUTS
% - P:      m x n 'dpvar' class object in variables var1 and/or var2,
%           of degree at most d1 in var1 and d2 in var2 (not cumulative);
% - var1:   p1 x 1 array of type 'polynomial' ('pvar') or 'cellstr'
%           representing the variables x = (x1,...,xp1);
% - var2:   p2 x 1 array of type 'polynomial' ('pvar') or 'cellstr'
%           representing the variables y = (y1,...,yp2);
% - dvars:  q x 1 array of type 'dpvar' or 'cellstr' representing the names
%           of the decision variables
% - deg1:   p1 x 1 array of integers specifying the maximal degree of the 
%           monomials in each variable xi to use in the basis Zd1;
% - deg2:   p2 x 1 array of integers specifying the maximal degree of the 
%           monomials in each variable yi to use in the basis Zd2;
%
% OUTPUTS
% - C:      m*r1*(q+1) x n*r2 array of type 'sparse', 
%           for r1 = (deg1(1)+1)*...*(deg1(p1)+1),
%           specifying the coefficients representing P in quadratic form.
% - Zd1:    m*r1 x m array of type 'polynomial' representing the monomial
%           basis (Im o Zd1(x));
% - Zd2:    n*r2 x n array of type 'polynomial' representing the monomial
%           basis (In o Zd2(y));
% - Zdec:   m*r1*(q+1) x m*r1 'cellstr' specifying the vector of decision 
%           variables, (Ik o [1;dvars]).



%% % Check the inputs
% Check the input polynomial.
if isa(P,'double')
    P = dpvar(P);
elseif isa(P,'polynomial')
    P = poly2dpvar(P);  % <-- ideally, we would not perform this conversion
elseif ~isa(P,'dpvar')
    error("Input polynomial should be specified as object of type 'dpvar' or 'polynomial'")
end
% Extract degrees and coefficients
[m,n] = size(P);
Pvars = P.varname;
Pdvars = P.dvarname;
Pdegs = P.degmat;
coeffs = P.C;

% Check the primary variables.
if isa(var1,'char')
    var1 = {var1};
elseif isa(var1,'polynomial')
    if ~ispvar(var1)
        error("Variables should be specified as array of 'pvar' objects or 'cellstr'.")
    end
    var1_name = cell(numel(var1),1);
    for ii=1:numel(var1)
        var1_name(ii) = var1(ii).varname;
    end
    var1 = var1_name;
elseif ~iscellstr(var1) && ~isempty(var1)
    error("Variables should be specified as array of 'pvar' objects or 'cellstr'.")
end
% Check the secondary variables
if isa(var2,'char')
    var2 = {var2};
elseif isa(var2,'polynomial')
    if ~ispvar(var2)
        error("Variables should be specified as array of 'pvar' objects or 'cellstr'.")
    end
    var2_name = cell(numel(var2),1);
    for ii=1:numel(var2)
        var2_name(ii) = var2(ii).varname;
    end
    var2 = var2_name;
elseif ~iscellstr(var2) && ~isempty(var2)
    error("Variables should be specified as array of 'pvar' objects or 'cellstr'.")
end
var1 = var1(:);         var2 = var2(:);
nvars1 = numel(var1);   nvars2 = numel(var2);
% Check the decision variables
if isa(dvars,'char')
    dvars = {dvars};
elseif isa(dvars,'dpvar')
    dvars_name = cell(numel(dvars),1);
    for ii=1:numel(dvars)
        dvar_ii = combine(dvars(ii));
        if ~isscalar(dvar_ii.dvarname) || nnz(dvar_ii.C)>1 || nonzeros(dvar_ii.C)~=1
            error("Decision variables should be specified as a qx1 array of individual 'dpvar' objects.")
        end
        dvars_name(ii) = dvar_ii.dvarname;
    end
    dvars = dvars_name;
elseif ~iscellstr(dvars) && ~isempty(dvars)
    error("Decision variables should be specified as array of 'dpvar' objects or 'cellstr'.")
end
dvars = dvars(:);

% Check the degrees
if nargin<=4
    if isempty(Pdegs)
        deg1 = zeros(size(var1));
    else
        deg1 = max(max(Pdegs))*ones(size(var1));
    end
elseif isscalar(deg1)
    deg1 = deg1*ones(size(var1));
elseif numel(deg1)~=numel(var1)
    error("Number of maximal degrees in primary variables should match number of primary variables.")
end
deg1 = deg1(:);
d1 = prod(deg1+1);
if nargin<=5
    if isempty(Pdegs)
        deg2 = zeros(size(var2));
    else
        deg2 = max(deg1)*ones(size(var2));
    end
elseif isscalar(deg2)
    deg2 = deg2*ones(size(var2));
elseif numel(deg2)~=numel(var2)
    error("Number of maximal degrees in secondary variables should match number of secondary variables.")
end
deg2 = deg2(:);
d2 = prod(deg2+1);

% Check that P depends only on the specified variables
nvars = nvars1 + nvars2;
if numel(Pvars)>0 && (numel(Pvars)>nvars || any(~ismember(Pvars,[var1;var2])))
    error("Input polynomial can depend only on the specified variables.")
end
q = numel(dvars);       q_P = numel(Pdvars);
if q_P>0 && (q_P>q || any(~ismember(Pdvars,dvars))) 
    error("Input object can depend only on the specified decision variables.")
end




%% % Compute the coefficient matrix C

% First, reorder the variable names of P to match var1 and var2
nZ = size(Pdegs,1);
Pdegs_new = zeros(nZ,nvars);
vars = [var1;var2];
for ii=1:nvars
    var_idx = ismember(Pvars,vars(ii));
    if any(var_idx)
        Pdegs_new(:,ii) = Pdegs(:,var_idx);
    end
end
Pdegs1 = Pdegs_new(:,1:nvars1);
Pdegs2 = Pdegs_new(:,nvars1+1:nvars);
if any(max(Pdegs1,[],1)>deg1') || any(max(Pdegs2,[],1)>deg2')
    error("Degree of monomials exceeds specified maximal degree.")
end

% Also check in what order the decision variables defining P appear in
% dvars:
% Pdvars = dvars(new2old_idcs)
new2old_idcs = zeros(q_P,1);
for ii=1:q_P
    new2old_idcs(ii) = find(ismember(dvars,Pdvars(ii)),1,'first');
end
new2old_idcs = [1;new2old_idcs+1];      % account for constant term

% Now, all coefficients in C come from a coefficient in P.coeff.
% --> extract the nonzero coefficient from P.coeff, as well as the row and
%       column numbers in P.coeff in which they appear.
[r_idcs_old,c_idcs_old,Cvals] = find(coeffs);

% We need to determine which row and column each element of Cvals appears
% in the new C matrix. 

% Determine to which row in P each nonzero element in coeffs belongs
r_nums = ceil(r_idcs_old(:)/(q_P+1));
% Determine to which decision variable in P each nonzero element belongs
Pdvar_nums = r_idcs_old(:) - (r_nums-1)*(q_P+1);
% Determine to which column in P each nonzero element belongs
c_nums = ceil(c_idcs_old(:)/nZ);
% Determine to which monomial in P each nonzero element belongs
Z_nums = c_idcs_old(:) - (c_nums-1)*nZ;

% Set the number of rows and columns in a block in C associated to the same
% row and column number in P:
% rows r_plus(i)+1:r_plus(i+1) in C correspond to row r_nums(i) in P
r_shift1 = (r_nums-1).*d1*(q+1);
c_shift1 = (c_nums-1).*d2;

% Set the number of rows in a block in C associated to the same power of 
% variable ii:
% rows 1 through nrows(ii) all involve the same power of var1(ii);
if isempty(deg1)
    nrows = zeros(0,1);
else
    nrows = cumprod(flipud([deg1+1;q+1]));
    nrows = flipud([nrows(1:end-1)]);
end
r_shift2 = Pdegs1*nrows;
r_idcs = r_shift1 + r_shift2(Z_nums) + new2old_idcs(Pdvar_nums);

% Set the number of columns in a block in C associated to the same power of 
% variable ii:
% columns 1 through ncols(ii) all involve the same power of var2(ii);
if isempty(deg2)
    ncols = zeros(0,1);
else
    ncols = cumprod(flipud(deg2+1));
    ncols = flipud([1;ncols(1:end-1)]);
end
c_shift2 = Pdegs2*ncols;
c_idcs = c_shift1 + c_shift2(Z_nums) +1;


% Declare the matrix C
C = sparse(r_idcs,c_idcs,Cvals,m*d1*(q+1),n*d2);
if nargout==1
    return
end

% % % Build the monomial bases
% Inefficient implementation for now...
Zd1 = eye(m);
for ii=1:nvars1
    Zd_ii = monomials(pvar(var1{ii}),0:deg1(ii));
    Zd1 = kron(Zd1,Zd_ii(:));
end
Zd2 = eye(n);
for ii=1:nvars2
    Zd_ii = monomials(pvar(var2{ii}),0:deg2(ii));
    Zd2 = kron(Zd2,Zd_ii(:));
end
k = size(Zd1,1);
Zdec = dpvar(eye(k*(q+1)),zeros(1,0),{},dvars,[k,k*(q+1)])';

end