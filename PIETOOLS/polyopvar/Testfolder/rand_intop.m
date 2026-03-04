function Kop = rand_intop(dim,var,dom,deg)
% RAND_INTOP(DIM,VAR,DOM,DEG) generates a random intop object of dimensions
% DIM, in variables VAR on domain DOM, and with parameters of monomial
% degree at most deg

% Extract the dummy variables used for integration
if isa(var,'polynomial')
    varname = var.varname(:)';
else
    varname = var(:)';
    var = polynomial(varname)';
end
nvars = numel(varname);

% Declare the order of integration in each term in the functional
omat_full = 1;
for j=2:nvars
    omat_tmp = omat_full;
    omat_full = zeros(0,j);
    nr = size(omat_tmp,1);
    for k=1:j
        omat_full = [omat_full; [omat_tmp(:,1:k-1), j*ones(nr,1), omat_tmp(:,k:end)]];
    end
end
n_rtn = randi(size(omat_full,1));
rtn_idcs = randi(size(omat_full,1),[1,n_rtn]);
omat = omat_full(unique(rtn_idcs),:);
ntrms = size(omat,1);

% Declare the kernel in each integral
params = rand_poly([dim(1),dim(2)*ntrms],var,deg)+1;

% Declare the functional
Kop = intop(params,omat,varname,dom);

end