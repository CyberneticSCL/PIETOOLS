function PIE = rescalePIE(PIE,I)

if nargin==1
    I = [-1,1];
end

PIE.T = transl(PIE.T,I);
PIE.Tw = transl(PIE.Tw,I);
PIE.Tu = transl(PIE.Tu,I);
PIE.A = transl(PIE.A,I);
PIE.B1 = transl(PIE.B1,I);
PIE.B2 = transl(PIE.B2,I);
PIE.C1 = transl(PIE.C1,I);
PIE.C2 = transl(PIE.C2,I);
PIE.D11 = transl(PIE.D11,I);
PIE.D12 = transl(PIE.D12,I);
PIE.D21 = transl(PIE.D21,I);
PIE.D22 = transl(PIE.D22,I);

if isfield(PIE,'L')
PIE.L = transl(PIE.L,I);
end
if isfield(PIE,'K')
PIE.K = transl(PIE.K,I);
end
end