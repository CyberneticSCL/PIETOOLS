% left shift monomials compute transformation of
% (I_p x Z1) G1 (I_q x Z2) G2 = (I x Z3) G3
% where Z1, Z2 are monomials
% G1,G2 are matrices

clc;clear;
rng(1)

pvar s1 s2 s3 s4 s5 s6 s7;
svar = [s1 s2 s3 s4  s5 s6 s7];

d = 4;    % dimension of monomial vector m1^d
m1 = 6; % size of Z1 m1^d
m2 = 5; % size of Z2 m2^d
p = 1;
q = 1;
deg = 7;

pG1 = p*m1^d;
pG2 = q*m2^d;
rG1 = q;
rG2 = 1;

var1 = sort(randsample(length(svar), d));
varA = svar(var1);
var2 = sort(randsample(length(svar), d));
varB = svar(var2);
Z1 = cell(1, d);
Z2 = cell(1, d);
for idx = 1:d
    Z1{idx} = sort(randsample(deg, m1))-1;
    Z2{idx} = sort(randsample(deg, m2))-1;
end
% Z1_full = monomials(svar(1:d), 0:deg).degmat;
% Z1_idx = sort(randsample(length(Z1_full), m1));
% Z2_idx = sort(randsample(length(Z1_full), m2));
% Z1 = Z1_full(Z1_idx, :);
% Z2 = Z1_full(Z2_idx, :);
G1 = sprand(pG1, rG1, 0.1);
G2 = sprand(pG2, rG2, 0.1);
% create monomial vector 
tic
[varsC2,ZC2,CC2] = leftShiftMonomials_old(varA.varname,Z1,{G1},varB.varname,Z2,{G2});
t2 = toc

tic
[varsC1,ZC1,CC1] = leftShiftMonomials_SS(varA.varname,Z1,{G1},varB.varname,Z2,{G2});
t1 = toc

isequal(varsC1, varsC2)
isequal(ZC1, ZC2)
max(abs(CC1{1} - CC2{1}), [], 'all')