function [A,B] = sync_basis(A,B)
% [A,B] = SYNC_BASIS(A,B) expresses the 'intop' objects A and B in terms of
% the same monomials Z.in and the same decision variables Zd, taking the
% unions of those of A and B, so that the parameters can be added.
%
% INPUTS
% - A:      'intop' object;
% - B:      'intop' object with the same number of spatial variables as A,
%           and the same degrees of these variables;
%
% OUTPUTS
% - A:      'intop' object representing the same functional as the input
%           A, expressed in terms of the monomials Z.in and the decision
%           variables Zd that A and B share;
% - B:      'intop' object representing the same functional as the input
%           B, with B.Z.in = A.Z.in and B.Zd = A.Zd;
%
% NOTES
% - If the monomials/decision variables already agree, nothing is done.
% - The unions are ordered by first appearance, first in A and then in B.
%

% DJ, 10/15/2026: Initial coding;

% Decision variables
if ~isequal(A.Zd,B.Zd)
    nA = numel(A.Zd);
    [Zd,~,ic] = unique([A.Zd;B.Zd],'stable');
    posA = ic(1:nA);        posB = ic(nA+1:end);
    nT = numel(Zd);
    mA = A.matdim(1);       mB = B.matdim(1);
    for k=1:numel(A.params)
        A.params{k} = remap_dvars(A.params{k},mA,posA,nT);
    end
    for k=1:numel(B.params)
        B.params{k} = remap_dvars(B.params{k},mB,posB,nT);
    end
    A.Zd = Zd;      B.Zd = Zd;
end

% Monomials
if ~isequal(A.Z.in,B.Z.in)
    N = numel(A.Z.in);
    Zu = cell(1,N);
    for i=1:N
        Zu{i} = unique([A.Z.in{i};B.Z.in{i}],'rows','stable');
    end
    [~,~,nZ] = Zdims(Zu);
    zmapA = Zmap(A.Z.in,Zu);        zmapB = Zmap(B.Z.in,Zu);
    for k=1:numel(A.params)
        A.params{k} = remap_cols(A.params{k},zmapA,nZ);
    end
    for k=1:numel(B.params)
        B.params{k} = remap_cols(B.params{k},zmapB,nZ);
    end
    A.Z.in = Zu;        B.Z.in = Zu;
end

end
