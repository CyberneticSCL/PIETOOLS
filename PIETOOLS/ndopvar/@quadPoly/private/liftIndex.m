function Cup = liftIndex(C, m, n, mapS, mapT, dsU, dtU)
%LIFTINDEX Lift coefficient matrix by reindexing nnz entries using maps.
% Lifts C such that 
%       (I_m\otimes Zs_old')*C*(I_n\otimes Zt_old) 
%               = (I_m \otimes Zs_union')*Cup*(I_n \otimes Zt_union)
% mapS: dsOld×1 old->union map (values in 1..dsU)
% i.e., Zs_old = Zs_union(mapS)
% mapT: dtOld×1 old->union map (values in 1..dtU)
% i.e., Zt_old = Zt_union(mapT)

C = sparse(C);

dsOld = numel(mapS);
dtOld = numel(mapT);

[i, j, v] = find(C);

a  = mod(i-1, dsOld) + 1;
ib = floor((i-1) / dsOld);

b  = mod(j-1, dtOld) + 1;
jb = floor((j-1) / dtOld);

i2 = ib*dsU + mapS(a);
j2 = jb*dtU + mapT(b);

Cup = sparse(i2, j2, v, m*dsU, n*dtU);
end
