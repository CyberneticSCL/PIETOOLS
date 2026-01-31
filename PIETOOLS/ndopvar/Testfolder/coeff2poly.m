function P = coeff2poly(C,dim,deg,vars)
% P = COEFF2POLY(C,DIM,DEG,VARS) returns a 'polynomial' object representing
% the function
%   P(s,t) = (Im o Z1(s))^T C (In o Z2(t))
% where [m,n] = DIM, s = VARS(:,1), t = VARS(:,2), and Z1 and Z2 are
% monomial basis of degree at most DEG in the variables s and t
%
% INPUTS
% - C:      m*d1 x n*d2 sparse array of coefficients, for 
%           d1=prod(degs(:,1) and d2 = prod(degs(:,2)+1);
% - dim:    1x2 array, [m,n], specifying the dimensions of P
% - deg:    N x 2 array, with element (i,j) specifying the maximal degree
%           of the monomials in variables vars(i,j);
% - vars:   N x 2 'polynomial' array specifying the variables (s1,...,sN) 
%           (column 1) and (t1,...,tN) (column 2);
%
% OUTPUTS
% - P:      m x n 'polynomial' object;
% 


% Determine the dimensions of the matrix-valued polynomial
m = dim(1);    n = dim(2);

if isempty(C)
    % Allow for empty coefficients
    P = polynomial(zeros(m,n));
    return
end

% Extract the polynomial variables
var1 = vars(:,1);
var2 = vars(:,2);
N = numel(var1);

% Get the degrees of the monomials
if isscalar(deg)
    deg = deg*ones(N,2);
elseif numel(deg)==N
    deg = [deg(:),deg(:)];
elseif ~all(size(deg)==[N,2])
    error("Monomial degrees should be specified as Nx2 array for N variables.")
end

% Build a full matrix of degrees in variables var1 and var2
degmat1 = zeros(1,0);       degmat2 = zeros(1,0);
varname1 = cell(N,1);       varname2 = cell(N,1);
for ii=1:N
    degmat1 = [kron(degmat1,ones(deg(ii,1)+1,1)), repmat((0:deg(ii,1))',[size(degmat1,1),1])];
    degmat2 = [kron(degmat2,ones(deg(ii,2)+1,1)), repmat((0:deg(ii,2))',[size(degmat2,1),1])];
    varname1(ii) = var1(ii).varname;
    varname2(ii) = var2(ii).varname;
end
% Build the monomial basis
%   Z1t = (Im o Z1(s))^T
nZ1 = size(degmat1,1);
Imat1 = sparse(repmat((1:nZ1)',[m,1]),m*(0:nZ1-1)'+(0:m-1)*m*nZ1+(1:m),1,nZ1,m*nZ1*m);
Z1t = polynomial(Imat1,degmat1,varname1,[m,nZ1*m]);

% Build the monomial basis
%   Z2 = (In o Z2(t))
nZ2 = size(degmat2,1);
Imat2 = sparse(repmat((1:nZ2)',[n,1]),(1:nZ2)'+(0:n-1)*nZ2,1,nZ2,n*nZ2*n);
Z2 = polynomial(Imat2,degmat2,varname2,[nZ2*n,n]);

% Declare the polynomial function
%   P(s,t) = (Im o Z1(s))^T C (In o Z2(t))
P = Z1t*C*Z2;

end