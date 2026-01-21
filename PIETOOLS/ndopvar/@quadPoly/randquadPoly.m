function F = randquadPoly(dim,nMons,var_s,var_t,maxdeg,density)
% F = quadPoly.randquadPoly([m,n],[ds,dt],var_s,var_t,density,degS,degT)

m = dim(1); n = dim(2);

ds = nMons(1);
dt = nMons(2);
degS = maxdeg(1);
degT = maxdeg(2);

Zs = randExps(ds, numel(var_s), degS);
Zt = randExps(dt, numel(var_t), degT);

C  = sprand(m*ds, n*dt, density);

ns = var_s;
nt = var_t;

F = quadPoly(C, Zs, Zt, [m n], ns, nt);
end

function Z = randExps(d, nvar, deg)
if isscalar(deg)
    Z = randi([0 deg], d, nvar);
else
    deg = deg(:).';
    Z = zeros(d,nvar);
    for j = 1:nvar
        Z(:,j) = randi([0 deg(j)], d, 1);
    end
end
end
