function [A,B,Zd] = CombineDecisionBasis(A,B)
% Combine and align decision variables for two sdopvar objects.
if ~isa(A,'sdopvar') || ~isa(B,'sdopvar')
    error('Inputs must be sdopvar objects.');
end
ZdA = cellstr(string(A.Zd(:)));
ZdB = cellstr(string(B.Zd(:)));
% Zd contains decision variables d \in A and 
% decision variables d\in B | d \nin A, in  the order that they appear in ZdB
Zd = [ZdA; setdiff(ZdB,ZdA,'stable')];
% change both decision variables in A and B to combined variables in Zd
A = ChangeDecVar(A,Zd);
B = ChangeDecVar(B,Zd);
end
