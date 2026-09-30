function Cout = lrmultiply_batch(L, R, C)
% Given L = {L1, L2, L3,...} and R = {R1, R2, R3,...}
% this function performs vec(L*C*R)

Lprod = cellfun(@prod, L);
Rprod = cellfun(@prod, R);

Cout = lrmultiply(Lprod,C,Rprod);
end