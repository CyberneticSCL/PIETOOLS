function [CC,ZdC] = plus_decparam(CA,CB,ZdA,ZdB)
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
Bout = sparse(numel(ZdC),size(CA.B,2));
Bout(mapA,:) = CA.B;
Bout(mapB,:) = Bout(mapB,:) + CB.B;

Aout = CA.A + CB.A;
CC = struct('A',Aout, 'B', Bout);
end