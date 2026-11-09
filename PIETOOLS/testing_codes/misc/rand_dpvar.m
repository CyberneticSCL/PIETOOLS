function P = rand_dpvar(dim,dvar,var,deg,nZ_max,rho)
% P = RAND_POLY(DIM,DVAR,VAR,DEG,NZ_MAX) randomly generates a matrix-valued
% 'dpvar' of dimensions DIM, and of degree at most DEG in variables VAR and
% decision variables DVAR;
% NZ_MAX is the maximal number of monomials to include in the polynomial
% (optional argument), independent of the number of variables or degrees.
% rho is the density of the coefficient matric defining the dpvar object,
% defaults to 2*numel(dvar)/numel(C);

% Deal with case of empty polynomial
if any(dim==0)
    P = dpvar(zeros(dim));
    return
end

% Check that the decision variables are properly specified
if isa(dvar,'double')
    % Generate q=dvar decision variable names
    q = dvar;
    dvar = cell(q,1);
    for i=1:q
        dvar{i} = ['coeff',num2str(i)];
    end
elseif ~iscellstr(dvar)
    error("Decision variables must be specified as 'cellstr'.")
end
dvar = dvar(:);
ndvars = numel(dvar);

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
    
% Extract the dimensions
m = dim(1);
n = dim(2);

% Randomly determine the number of monomials.
if nargin<=4
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

% Generate random coefficients of desired sparsity
if nargin<=5
    %rho = 2*ndvars/((m*ndvars+1)*n*nZ);
    rho = 2/(nZ*ndvars);
end
C = sprand(m*(ndvars+1),n*nZ,rho);

% Set the polynomial
P = dpvar(C,degmat,var,dvar,[m,n]);

end

