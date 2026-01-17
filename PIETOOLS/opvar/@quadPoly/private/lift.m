% @quadPoly/private/lift.m
function Clift = lift(C, m, n, dsOld, dtOld, dsNew, dtNew, mapS, mapT)
[Iold, Jold, Vold] = find(C);
if isempty(Vold)
    Clift = sparse(m*dsNew, n*dtNew);
    return;
end

br = floor((Iold-1)/dsOld) + 1;
ls = mod(Iold-1, dsOld) + 1;
bc = floor((Jold-1)/dtOld) + 1;
lt = mod(Jold-1, dtOld) + 1;

Inew = (br-1)*dsNew + mapS(ls);
Jnew = (bc-1)*dtNew + mapT(lt);

Clift = sparse(Inew, Jnew, Vold, m*dsNew, n*dtNew);
end
