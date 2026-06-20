function P = rand_poly(dim,var,deg,nZ_max)
% P = RAND_POLY(DIM,VAR,DEG,NZ_MAX) randomly generates a matrix-valued
% polynomial of dimensions DIM, and of degree at most DEG in variables VAR.
% NZ_MAX is the maximal number of monomials to include in the polynomial
% (optional argument), independent of the number of variables or degrees.

% Deal with case of empty polynomial
if any(dim==0)
    P = polynomial(zeros(dim));
    return
end

% Check that the variables are properly specified
if isa(var,'double')
    % Generate N=var variable names
    N = var;
    var = cell(N,1);
    for i=1:N
        var{i} = ['s',num2str(i)];
    end
elseif isa(var,'char')
    var = {var};
elseif isa(var,'polynomial')
    if ~ispvar(var) && ~isempty(var)
        error("Variables must be specified as 'pvar' objects or 'cellstr'.")
    end
    var = pvar2varname(var);
end
var = var(:);
nvars = numel(var);

% Check that the degrees are properly specified.
if isscalar(deg)
    deg = deg*ones(1,nvars);
elseif numel(deg)~=nvars
    error("Number of maximal monomial degrees should match number of variables.")
end
    

m = dim(1);
n = dim(2);

% Randomly determine the number of monomials.
if nargin<=3
    nZ_max = prod(deg+1);
else
    nZ_max = min(prod(deg+1),nZ_max);
end
nZ = randi(nZ_max);
% Generate a random monomial basis.
degmat = zeros(nZ,nvars);
for jj=1:nvars
    degmat(:,jj) = randi(deg(1)+1,nZ,1)-1;
end
degmat = uniquerows_integerTable(degmat);
nZ = size(degmat,1);

% Generate random coefficients, on average two nonzero coefficients per
% element of the matrix-valued object
%is_zero = logical(randi(floor(2*(m*n+1)),[nZ,m*n])-1);
%C = (~is_zero).*(2*rand([nZ,m*n])-1);
dnsty = 2/nZ;
C = sprand(nZ,m*n,dnsty);

% Set the polynomial
P = polynomial(C,degmat,var,[m,n]);

end

