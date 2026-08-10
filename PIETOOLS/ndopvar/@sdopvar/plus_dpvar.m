function [CC,ZdC] = plus_dpvar(CA,CB,ZdA,ZdB)
% Given two decision parameter structures of 
% CA = CA.A + CA.Bt*ZdA
% CB = CB.A + CB.Bt*ZdB,
% this function returns the objects representing
% the following
% CC = CC.C + CC.Bt*ZdC
nA = numel(ZdA);

[ZdC,~,ic] = unique([ZdA(:); ZdB(:)]);
mapA = ic(1:nA);
mapB = ic(nA+1:end);


% Build coefficient matrix over common parameter set
Btout = sparse(size(CA.Bt,1), numel(ZdC));
Btout(:,mapA) = CA.Bt;
Btout(:,mapB) = Btout(:,mapB) + CB.Bt;

Aout = CA.A + CB.A;
CC = struct('A',Aout, 'Bt', Btout);
end