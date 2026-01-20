function F = randquadPoly(dim,nMons,vars,maxdeg,density)
% F = quadPoly.randquadPoly(m,n,ds,dt,nsVar,ntVar,density,degS,degT)
% degS,degT: scalar max degree or 1xnsVar / 1xntVar vectors (default 3)

m = dim(1); n = dim(2);

nsVar = vars(1); ntVar = vars(2);
ds = nMons(1);
dt = nMons(2);
degS = maxdeg(1);
degT = maxdeg(2);

Zs = randExps(ds, vars(1), degS);
Zt = randExps(dt, vars(2), degT);

C  = sprand(m*ds, n*dt, density);

ns = arrayfun(@(i) sprintf('s%d',i), 1:nsVar, 'UniformOutput', false);
nt = arrayfun(@(i) sprintf('t%d',i), 1:ntVar, 'UniformOutput', false);

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
