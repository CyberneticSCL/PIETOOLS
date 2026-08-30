function [A,B,Zd] = CombineDecisionBasis(A,B)
% Combine and align decision variables for two sdopvar objects.
if ~isa(A,'sdopvar') || ~isa(B,'sdopvar')
    error('Inputs must be sdopvar objects.');
end
ZdA = A.Zd.';
ZdB = B.Zd.';
% Zd contains decision variables d \in A and 
% decision variables d\in B | d \nin A, in  the order that they appear in ZdB
ZdB_new = setdiff(ZdB,ZdA,'stable');
Zd = [ZdA(:); ZdB_new(:)];
% change both decision variables in A and B to combined variables in Zd
A = ChangeDecVar(A,Zd);
B = ChangeDecVar(B,Zd);
end
