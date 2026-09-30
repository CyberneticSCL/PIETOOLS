function [A,B,Zd] = CombineDecisionBasis(A,B)
% Combine and align decision variables for two sdopvar objects.
if ~isa(A,'sdopvar') || ~isa(B,'sdopvar')
    error('Inputs must be sdopvar objects.');
end
ZdA = A.Zd.';
ZdB = B.Zd.';
% Zd contains decision variables d \in A and 
% decision variables d\in B | d \nin A, in  the order that they appear in ZdB
% When the two lists already agree, which is the normal case once the
% operands have been put on a common basis, skip the setdiff: it costs a
% string comparison per decision variable and dominates for large nd.
if isequal(ZdA,ZdB)                                                         % MMP, 08/29/2026
    Zd = ZdA(:);                                                            % MMP, 08/29/2026
    A = ChangeDecVar(A,Zd);     B = ChangeDecVar(B,Zd);                     % MMP, 08/29/2026
    return                                                                  % MMP, 08/29/2026
end                                                                         % MMP, 08/29/2026
ZdB_new = setdiff(ZdB,ZdA,'stable');
Zd = [ZdA(:); ZdB_new(:)];
% change both decision variables in A and B to combined variables in Zd
A = ChangeDecVar(A,Zd);
B = ChangeDecVar(B,Zd);
end
