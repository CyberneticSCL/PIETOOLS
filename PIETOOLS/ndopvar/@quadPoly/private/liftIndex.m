function Cup = liftIndex(C, m, n, mapS, mapT, dsU, dtU)
%LIFTINDEXMAPS Lift coefficient matrix by reindexing nnz entries using maps.
%
% mapS: dsOld×1 old->union map (values in 1..dsU)
% mapT: dtOld×1 old->union map (values in 1..dtU)

C = sparse(C);

dsOld = numel(mapS);
dtOld = numel(mapT);

[i, j, v] = find(C);

a  = mod(i-1, dsOld) + 1;
ib = floor((i-1) / dsOld);

b  = mod(j-1, dtOld) + 1;
jb = floor((j-1) / dtOld);

i2 = ib * dsU + mapS(a);
j2 = jb * dtU + mapT(b);

Cup = sparse(i2, j2, v, m*dsU, n*dtU);
end
