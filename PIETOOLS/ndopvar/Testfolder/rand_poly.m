function P = rand_poly(dim,var,deg)
% P = RAND_POLY(DIM,VAR,DEG) randomly generates a matrix-valued polynomial
% of dimensions DIM, and of degree at most DEG in variables VAR

% Check that the variables are properly specified
if isa(var,'char')
    var = {var};
elseif isa(var,'polynomial')
    if ~ispvar(var)
        error("Variables must be specified as 'pvar' objects or 'cellstr'.")
    end
    var_name = cell(numel(var),1);
    for ii=1:numel(var)
        var_name(ii) = var(ii).varname;
    end
    var = var_name;
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
nZ_max = prod(deg+1);
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
is_zero = logical(randi(floor(0.5*(m*n+1)),[nZ,m*n])-1);
C = (~is_zero).*(2*rand([nZ,m*n])-1);

% Set the polynomial
P = polynomial(C,degmat,var,[m,n]);

end

