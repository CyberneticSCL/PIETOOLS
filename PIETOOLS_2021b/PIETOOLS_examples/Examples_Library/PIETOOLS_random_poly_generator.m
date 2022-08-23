function P = PIETOOLS_random_poly_generator(nr,nc,varname,dmax)
%% Random matrix-valued polynomial constructor
% Polynomial will consist of n rows and m columns
% Polynomial will exist on variables vars
% The polynomial will have maximal degree dmax

% Convert pvar to cellstr.
if isempty(varname)
    varname = cell(0,1);
elseif ispvar(varname)
    varname = varname.varname;
end

nvars = length(varname);
nm = nr*nc;

coeff = randi([0,1],dmax,nm).*randi([0,dmax+1],dmax,nm);
degmat = randi([0 dmax],dmax,nvars);
matdim = [nr nc];
P = polynomial(coeff,degmat,varname,matdim);

end