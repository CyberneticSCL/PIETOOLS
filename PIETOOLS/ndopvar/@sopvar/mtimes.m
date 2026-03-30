function C = mtimes(A,B)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% C = mtimes(A,B) composes two sopvar operators C=AB where the output
% domain of B must be the input domain of A
% 
% Date: 1/16/26
% Version: 1.0
% 
% INPUT
% A: sopvar class object L_2^q[S4,S5] to L_2^r[S6,S5]
% B: sopvar class object L_2^p[S1,S3] to L_2^q[S2,S3]
%
% OUTPUT
% C = AB:  L_2^p[S7,S8] to L_2^r[S9,S8]
% 
% Where S7= S1, S3 cap S4                   %S1=S1, S3b= S3 cap S4 
%       S8= S3 cap S5                       %S3a= S3 cap S5 
%       S9= S6, S2 cap S5                   %S4=S6 ,S2a= S2 cap S5 
% NOTES:
% For support, contact M. Peet, Arizona State University at mpeet@asu.edu
% or S. Shivakumar at sshivak8@asu.edu



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PIETools - mtimes
%
% Copyright (C)2026  M. Peet, S. Shivakumar
%
% This program is free software; you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation; either version 2 of the License, or
% (at your option) any later version.
%
% This program is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
%
% You should have received a copy of the GNU General Public License
% along with this program; if not, write to the Free Software
% Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% If you modify this code, document all changes carefully and include date
% authorship, and a brief description of modifications
%
% Initial coding MMP, SS  - 1_16_2026
%
% Error handling: Checks to ensure A and B are compatible
if A.dims(2)~=B.dims(1)
    error('number of output components of B is different than number of input components of A');
end
if any(~strcmp(A.vars.in,B.vars.out))
    error('names of B output variables differ from names of A input variables');
end
if any(any(A.dom.in~=B.dom.out))
    error('output domain of B differs from input domain of A');
end

%B [s1,s3a,s3b] -> [s2a,s2b,s3a,s3b]
%A [s2a,s2b,s3a,s3b] -> [s4,s2a,s3a]
vs1 = B.vars_S1;
vs3a = intersect(A.vars_S3,B.vars_S3);
vs3b = setdiff(B.vars_S3,vs3a);
vs4 = A.vars_S2;
vs2a = setdiff(A.vars_S3,vs3a);
vs2b = setdiff(A.vars_S1,vs3b);

Cvars = struct('in',B.vars.in,'out',A.vars.out);
Cdom = struct('in',B.dom.in,'out',A.dom.out);
Cdims = [A.dims(1),B.dims(2)];

n3 = numel(vs3a);

if n3 == 0  % params is 1x1 cell
    Cparams = {zeros(A.dims(1), B.dims(2))};
elseif n3 == 1  % params is 3x1 cell
    Cparams = repmat({zeros(A.dims(1), B.dims(2))}, 3, 1);   % <-- key change
else % params is 3x3x...x3 cell
    Cparams = repmat({zeros(A.dims(1), B.dims(2))}, 3*ones(1,n3));
end

% Outline
% 1. Extract A_\alpha parameters from A and B_\beta parameters from B 
% 2. Extract Monomials and Coefficients from each set of parameters: C^A_\alpha and C^B_\beta
% 3a. Compute the D_{\alpha_a,\beta_b} 
% 3b. Compute the G_\alpha,\beta terms and express as Z C_{\alpha,\beta}^G Z 
% 4. Redistribute C^A_\alpha G_\alpha,\beta C^B_\beta as Z C^C_{\alpha,\beta} Z
% 5. Map \sum_{\alpha,\beta} C^C_{\alpha,\beta} to \sum_\gamma C^C_\gamma
% End: construct quadpolys from Z1,Z2, C^C_\gamma

% Now, just need to compute the parameters!
% I would suggest looping over the gammas and creating a set of terms to
% add for each gamma, since for some alpha,betas they get added to two of
% the gammas
nalpha = numel(B.vars_S3);
nbeta = numel(A.vars_S3);
alphaIdx = zeros(numel(B.params),nalpha);
betaIdx = zeros(numel(A.params),nbeta);

if nalpha==0
    Bsize = [1,1];
elseif nalpha==1
    Bsize = [3,1];
else
    Bsize = repmat(3,1,nalpha);
end
if nbeta==0
    Asize = [1,1];
elseif nbeta==1
    Asize = [3,1];
else
    Asize = repmat(3,1,nbeta);
end
for i=1:numel(B.params)
    tmp = cell(1,nalpha);
    [tmp{:}] = ind2sub(Bsize(:),i);
    alphaIdx(i,:) = [tmp{:}];
end
for i=1:numel(A.params)
    tmp = cell(1,nbeta);
    [tmp{:}] = ind2sub(Asize(:),i);
    betaIdx(i,:) = [tmp{:}];
end

% split alpha and beta as [alpha_a, alpha_b] and [beta_a, beta_b]
[idx2a,~] = ismember(A.vars_S3,v2a);
beta_a = betaIdx(:,idx2a);
beta_b = betaIdx(:,~idx2a);
[idx3a,~] = ismember(B.vars_S3,v3a);
alpha_a = alphaIdx(:,idx3a);
alpha_b = alphaIdx(:,~idx3a);

[idx2b,~] = ismember(B.vars_S2,v2b);
dom2b = B.dom_2(idx2b,:);

% do int(ZL*ZR', s2b) = (I⊗Zhat_s2a')*G_s3a(s3a)*(I⊗Zhat_s3b)
[Zhat_s2a,G_s3a,Zhat_s3b] = ...
        int_2b(A.ZR,B.ZL,strrep(A.nt,'t','s'),B.ns, vs2a,vs2b,vs3a,vs3b, dom2b);

[Cs2a, KZhat_betaa_s2a] = int_monomial(Zhat_s2a,beta_a);
[Cs3b, KZhat_alphab_s3b] = int_monomial(Zhat_s3b,alpha_b);

[CZL, CZR, Cparams] = termCompose(A,B,alphaIdx,betaIdx);

C = sopvar(Cparams,Cvars,CZR,CZL,Cdom,Cdims);
end


function [Z2aout, G3aout, Z3bout] = int_2b(ZL, ZR, ZLvar, ZRvar, s2a,s2b,s3a,s3b,lims)
% This performs the factorization
% int(ZL*ZR',s2b,0,1) = (Im\otimes Z2a') G(s3a) (In\otimes Z3b)
% where m = length(kron(ZL)) and n = length(kron(ZR))

% note that
% int_0^1 (Z2a⊗Z2b⊗Z3a⊗Z3b)(Zp2a⊗Zp2b⊗Zp3a⊗Zp3b)' d2b
% = int_0^1 (Z2a*Zp2a')⊗(Z2b*Zp2b')⊗(Z3a*Zp3a')⊗(Z3b*Zp3b') d2b
% = (Z2a*Zp2a')⊗C2b⊗(Z3a*Zp3a')⊗(Z3b*Zp3b')
% for some constant matrix C2b

% first separate ZL and ZR into 2a,2b,3a,3b
ZL2a = cell(1,length(s2a));
ZR2a = cell(1,length(s2a));
ZL2b = cell(1,length(s2b));
ZR2b = cell(1,length(s2b));
ZL3a = cell(1,length(s3a));
ZR3a = cell(1,length(s3a));
ZL3b = cell(1,length(s3b));
ZR3b = cell(1,length(s3b));

for i=1:length(s2a)
    [idx,loc] = ismember(s2a(i),ZLvar);
    if ~idx, ZL2a{i} = 0; else, ZL2a{i} = ZL{loc}; end
    [idx,loc] = ismember(s2a(i),ZRvar);
    if ~idx, ZR2a{i} = 0; else, ZR2a{i} = ZR{loc};
    end
end
for i=1:length(s2b)
    [idx,loc] = ismember(s2b(i),ZLvar);
    if ~idx, ZL2b{i} = 0; else, ZL2b{i} = ZL{loc}; end
    [idx,loc] = ismember(s2b(i),ZRvar);
    if ~idx, ZR2b{i} = 0; else, ZR2b{i} = ZR{loc};
    end
end
for i=1:length(s3a)
    [idx,loc] = ismember(s3a(i),ZLvar);
    if ~idx, ZL3a{i} = 0; else, ZL3a{i} = ZL{loc}; end
    [idx,loc] = ismember(s3a(i),ZRvar);
    if ~idx, ZR3a{i} = 0; else, ZR3a{i} = ZR{loc};
    end
end
for i=1:length(s3b)
    [idx,loc] = ismember(s3b(i),ZLvar);
    if ~idx, ZL3b{i} = 0; else, ZL3b{i} = ZL{loc}; end
    [idx,loc] = ismember(s3b(i),ZRvar);
    if ~idx, ZR3b{i} = 0; else, ZR3b{i} = ZR{loc};
    end
end

% perform ZL*ZR' for each variable group 2a,2b,3a,3b
for i=1:length(s2a)
    Z2a{i} = ZL2a{i} + ZR2a{i}';
end
for i=1:length(s2b)
    Z2b{i} = ZL2b{i} + ZR2b{i}';
end
for i=1:length(s3a)
    Z3a{i} = ZL3a{i} + ZR3a{i}';
end
for i=1:length(s3b)
    Z3b{i} = ZL3b{i} + ZR3b{i}';
end

% integrate Z2b part over [a,b]
% int_a^b Z2b d2b = 
%       kron_i (b(i)^(Z2b{i}+1)-a(i)^(Z2b{i}+1))/(Z2b{i}+1)
C2b = 1; a = lims(:,1); b = lims(:,2);
for i=1:length(s2b)
    E = Z2b{i};
    Ci = (b(i).^(E+1)-a(i).^(E+1))./(E+1);
    C2b = kron(C2b,Ci);
end

% condense Z2a matrix into (I2a \otimes Z2ap')*K2a
Z2ap = cell(1,length(s2a));
K2a = 1;
for i=1:length(s2a)
    E = Z2a{i};
    u = unique(E(:),'stable');
    Z2ap{i} = u;
    Ki = sparse(length(u), numel(E));
    for j=1:numel(E)
        Ki(find(u==E(j),1), j) = 1;
    end
    K2a = kron(K2a,Ki);
end

% condense Z3b matrix into K3b*(I3b \otimes Z3bp)
Z3bp = cell(1,length(s3b));
K3b = 1;
for i=1:length(s3b)
    E = Z3b{i};
    u = unique(E(:),'stable');
    Z3bp{i} = u;
    Ki = sparse(numel(E), length(u));
    for j=1:numel(E)
        Ki(j,find(u==E(j),1)) = 1;
    end
    K3b = kron(K3b,Ki);
end

% condense Z3a matrix into (I3a0 ⊗ Z3anew')*K3a
Z3anew = cell(1,length(s3a));
K3a = 1;
for i=1:length(s3a)
    E = Z3a{i};
    u = unique(E(:),'stable');
    Z3anew{i} = u;
    Ki = sparse(length(u), numel(E));
    for j=1:numel(E)
        Ki(find(u==E(j),1), j) = 1;
    end
    K3a = kron(K3a,Ki);
end

% now we have 
% (Z2a*Zp2a')⊗C2b⊗(Z3a*Zp3a')⊗(Z3b*Zp3b')
%  = ((I2a \otimes Z2anew')*K2a)
%      ⊗C2b
%        ⊗(I3a0 ⊗ Z3anew')*K3a
%          ⊗(K3b*(I3b \otimes Z3bnew))

% let us focus on finding
% K2a⊗C2b⊗((I3a0 ⊗ Z3anew')*K3a)⊗K3b
tempA = kron(K2a,C2b);
tempB = K3b;
r3a  = prod(cellfun(@length,ZL3a));
c3a  = prod(cellfun(@length,ZR3a));
nz3a = prod(cellfun(@length,Z3anew));

mA = size(A,1);  nA = size(A,2);
mB = size(B,1);  nB = size(B,2);

% perfect shuffle: maps [q,nB] ordering to [nB,q] ordering
q = nz3a;
p = nB;
idx = reshape(1:q*p,[q,p]);
idx = reshape(idx.',1,[]);
P3a = sparse(1:q*p, idx, 1, q*p, q*p);

I3a = speye(mA*r3a*mB);

C3a = kron(kron(kron(A,speye(r3a)),speye(mB)),speye(nz3a)) ...
    * kron(speye(nA*r3a), P3a) ...
    * kron(kron(speye(nA),K3a),speye(nB));


% now we have
% (Z2a*Zp2a')⊗C2b⊗(Z3a*Zp3a')⊗(Z3b*Zp3b')
%  = ((I2a⊗Z2anew')⊗I⊗I⊗I)
%      *(I3a⊗Z3anew')*C3a
%        *(I⊗I⊗I⊗Z3bnew)
% we need (Im⊗ Z2anew')*(I3a⊗Z3anew')*C3a*(In⊗Z3bnew)

% first, we find permuation P2a such that
% ((I2a⊗Z2anew')⊗I⊗I⊗I) = (Im⊗ Z2anew')*P2a
% (I⊗I⊗I⊗Z3bnew) does not need any permutation
% sizes
n2a = prod(cellfun(@length,Z2anew));   % size of Z2anew block
n3a = prod(cellfun(@length,Z3anew));   % size of Z3anew block
nI3a = size(C3a,1) / n3a;              % row identity size in (I3a ⊗ Z3anew')*C3a
nrest = nI3a;                          % same block that sits to the right of Z2anew

idx = reshape(1:n2a*nrest,[n2a,nrest]);
idx = reshape(idx.',1,[]);   % perfect shuffle on [iz2a, irest]

Pblk = sparse(1:n2a*nrest, idx, 1, n2a*nrest, n2a*nrest);

nI2a = prod(cellfun(@length,ZL2a));    % size of I2a block
P2a = kron(speye(nI2a), Pblk);

% We want:
% P2a*(I3a ⊗ Z3anew') = (I3a_new ⊗ Z3anew')*Q2a
%
% Since P2a only permutes identity-type indices and does not touch Z3anew,
% Q2a is just the induced permutation on the identity-row block of C3a,
% lifted across the untouched Z3anew dimension.

% first extract the permutation vector from P2a
p_outer = zeros(size(P2a,1),1);
for j = 1:size(P2a,2)
    p_outer(j) = find(P2a(:,j),1);
end

% reshape identity rows against the untouched Z3anew block
tmp = reshape(1:nI3a*n3a, [nI3a, n3a]);
tmp = tmp(p_outer, :);
q = tmp(:);
Q2a = sparse(1:length(q), q, 1, length(q), length(q));

C3a = Q2a*C3a;

%finally, we have 
% (Im⊗ Z2anew')*(I3anew⊗Z3anew')*C3a*(In⊗Z3bnew)
Z2aout = Z2anew;
Z3bout = Z3bnew;
G3aout = struct('C', C3a, 'Z', Z3anew);
end