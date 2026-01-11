function [C,Zd1,Zd2] = get_quadratic_form(P,var1,var2,deg1,deg2)
% [C,ZD1,ZD2] = GET_QUADRATIC_FORM(P,VAR1,VAR2,DEG1,DEG2)
% Compute coefficients C such that P(x,y)=(Im o Zd1(x))^T C (In o Zd2(y)),
% where x is var1 and y is var2, [m,n] are the dimensions of P,
% o denotes the Kronecker product, and
%   Zd1(x) = Zd1_{1}(x1) o Zd1_{2}(x2) o ... o Zd1_{p}(xp)
% for Zd1_{i}(xi) the vector of monomials of degree at most deg1(i) in xi;
% INPUTS
% - R:      m x n 'polynomial' class object in variables var1 and/or var2,
%           of degree at most d1 in var1 and d2 in var2 (not cumulative)
% - var1:   p x 1 array of type 'polynomial' ('pvar') or 'cellstr'
%           representing the variables x = (x1,...,xp);
% - var2:   q x 1 array of type 'polynomial' ('pvar') or 'cellstr'
%           representing the variables y = (y1,...,yp);
% - deg1:   p x 1 array of integers specifying the maximal degree of the 
%           monomials in each variable xi to use in the basis Zd1;
% - deg2:   q x 1 array of integers specifying the maximal degree of the 
%           monomials in each variable yi to use in the basis Zd2;
%
% OUTPUTS
% - C:      m*d1 x n*d2 array of type 'sparse', 
%           for d1 = deg1(1)*...*deg1(p),
%           specifying the coefficients representing P in quadratic form.
% R(x,y)=(Im o Zd(x))^T C (In o Zd(y))



%% % Check the inputs
% Check the input polynomial.
if isa(P,'double')
    P = polynomial(P);
elseif ~isa(P,'polynomial')
    error("Input polynomial should be specified as object of type 'polynomial'")
end
% Extract degrees and coefficients
[m,n] = size(P);
Pvars = P.varname;
Pdegs = P.degmat;
coeffs = P.coeff;

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
elseif ~iscellstr(var2) && ~isempty(var2);
    error("Variables should be specified as array of 'pvar' objects or 'cellstr'.")
end
var1 = var1(:);         var2 = var2(:);
nvars1 = numel(var1);   nvars2 = numel(var2);

% Check the degrees
if nargin<=3
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
if nargin<=4
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
if numel(Pvars)>nvars || any(~ismember(Pvars,[var1;var2]))
    error("Input polynomial can depend only on the specified variables.")
end


%% % Compute the coefficient matrix C

% Re-order the variable names of P to match var1 and var2
Pdegs_new = zeros(size(Pdegs,1),nvars);
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


% Set the row and column number in the matrix-valued P associated with each
% column of P.coeff
r_nums = repmat((1:m),[1,n]);
c_nums = kron((1:n),ones(1,m));
% Set the number of rows and columns in a block in C associated to the same
% row and column number in P:
% rows r_plus(i)+1:r_plus(i+1) in C correspond to row r_nums(i) in P
r_plus = (r_nums-1).*d1;
c_plus = (c_nums-1).*d2;
% Set the number of rows and columns in a block in C associated to 
% the same power of variable ii:
% rows 1 through nrows(ii) all involve the same power of var1(ii);
% columns 1 through ncols(ii) all involve the same power of var2(ii);
if isempty(deg1)
    nrows = zeros(0,1);
else
    nrows = cumprod(flipud(deg1+1));
    nrows = flipud([1;nrows(1:end-1)]);
end
if isempty(deg2)
    ncols = zeros(0,1);
else
    ncols = cumprod(flipud(deg2+1));
    ncols = flipud([1;ncols(1:end-1)]);
end

% Extract the nonzero coefficients from P.coeff, and establish the
% corresponding row and column numbers in the matrix C
Cvals = [];
r_idcs = [];
c_idcs = [];
for ii=1:size(coeffs,1)
    % Get the coefficients associated with the iith monomial.
    Cvals_ii = coeffs(ii,:);
    is_nonzero = Cvals_ii~=0;
    Cvals_ii = Cvals_ii(is_nonzero);
    Cvals = [Cvals; Cvals_ii(:)];
    % Get the row numbers in C of each of the coefficients.
    Pdegs1_ii = Pdegs1(ii,:);
    r_idcs_ii = r_plus(is_nonzero) + Pdegs1_ii*nrows + 1;
    r_idcs = [r_idcs; r_idcs_ii(:)];
    % Get the column numbers in C of each of the coefficients.
    Pdegs2_ii = Pdegs2(ii,:);
    c_idcs_ii = c_plus(is_nonzero) + Pdegs2_ii*ncols + 1;
    c_idcs = [c_idcs; c_idcs_ii(:)];
end
% Declare the matrix C
C = sparse(r_idcs,c_idcs,Cvals,m*d1,n*d2);
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

end