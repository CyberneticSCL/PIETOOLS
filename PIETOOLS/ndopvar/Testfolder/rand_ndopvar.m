function Pop = rand_ndopvar(dim,deg,dom,var1,var2,dvarname)
% RAND_NDOPVAR Generates a random ndopvar object operator with monomials of 
% degree at most deg in each of the variables

Pop = opvar2d(); 
if numel(dim)~=2
    error("Only 3^N-PI operators are supported.")
end
m = dim(1);     n = dim(2);
N = size(dom,1);
q = numel(dvarname);
if isscalar(deg)
    deg = deg*ones(N,1);
end
d = prod(deg+1);

% Initialize an empty operator
Pop = ndopvar();
Pop.deg = deg;
Pop.dom = dom;
Pop.vars = [var1,var2];
Pop.dvarname = dvarname;

% Generate random coefficients
Pop.C = cell([3*ones(1,N),1]);
rdim = m*d*(q+1);
sz_C = size(Pop.C);
for ii=1:numel(Pop.C);
    % Determine the index of element ii along each dimension of the cell C
    idcs = cell(1,N);
    [idcs{:}] = ind2sub(sz_C,ii);
    idcs = cell2mat(idcs);
    % If element ii corresponds to an integral, we need to account for the
    % monomial basis in the associated dummy variable
    is_int = logical(idcs-1);
    cdim = n*prod(deg(is_int)+1);
    % Set sparse coefficients of dimension rdim x cdim
    rho = (q+10)/(rdim*cdim);
    Pop.C{ii} = sprand(rdim,cdim,rho);
end


end

