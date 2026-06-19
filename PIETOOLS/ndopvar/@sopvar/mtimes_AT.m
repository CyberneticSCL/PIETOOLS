function C = mtimes(A,B)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% C = mtimes(A,B) composes two sopvar operators C=A*B where the output
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

% Compute the domain of the composition
Cdom = struct('in',B.dom.in,'out',A.dom.out);

% Compute the dimensions of the composition
Cdims = [A.dims(1),B.dims(2)];

n3 = numel(vs3a); % size of passthrough variable

% initialize the parameter cell of appropriate size
if n3 == 0  % params is 1x1 cell
    Cparams = {zeros(A.dims(1), B.dims(2))};
elseif n3 == 1  % params is 3x1 cell
    Cparams = repmat({zeros(A.dims(1), B.dims(2))}, 3, 1);   
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

nalpha = numel(B.vars_S3); % pass through variables in B, s3a+s3b
nbeta = numel(A.vars_S3); % pass through variables in A, s2a+s3a
alphaIdx = cell(1,nalpha); % map parameter to operator time index
betaIdx = cell(1,nbeta);

Csize = size(Cparams);
Bsize = size(B.params);
Asize = size(A.params);

% convert linear index to multi-index for B parameters, e.g, if B is 3x3
% map 1,2...9 to [1,1],[2,1],...[3,3]
[alphaIdx{:}] = ind2sub(Bsize,1:numel(B.params));
if nalpha==0
    alphaIdx = [];
elseif nalpha==1
    alphaIdx = cell2mat(alphaIdx).';
else
    alphaIdx = cell2mat(alphaIdx(:)).';
end
[betaIdx{:}] = ind2sub(Asize,1:numel(A.params));
if nbeta==0
    betaIdx = [];
elseif nbeta==1
    betaIdx = cell2mat(betaIdx).';
else
    betaIdx = cell2mat(betaIdx(:)).';
end

% split alpha and beta as [alpha_a, alpha_b] and [beta_a, beta_b]
[idx2a,~] = ismember(A.vars_S3,vs2a);
beta_a = unique(betaIdx(:,idx2a),'rows','stable');
beta_b = unique(betaIdx(:,~idx2a),'rows','stable');
[idx3a,~] = ismember(B.vars_S3,vs3a);
alpha_a = unique(alphaIdx(:,idx3a),'rows','stable');
alpha_b = unique(alphaIdx(:,~idx3a),'rows','stable');

% REMEMBER beta_a, beta_b, alpha_a, alpha_b no longer have unique rows
% above. Need to revisit functions below to ensure this is handled
% appropriately!!!

[idx2b,~] = ismember(B.vars_S2,vs2b);
dom2b = B.dom_2(idx2b,:);
dom2a = A.dom_3(idx2a,:);
dom3a = B.dom_3(idx3a,:);
dom3b = B.dom_3(~idx3a,:);

% let len(A.ZR) = p, len(B.ZL) = q.
% do int(A.ZR*B.ZL', s2b) =
% (I_p⊗PL')(I_p⊗Zhat_s2a')*G_s3a(s3a)*(I_q⊗Zhat_s3b)(I_q⊗PR)
% where G_s3a is size len(Zhat_s2a)*p x len(Zhat_s3b)*q = g1 x g2
[Zhat_s2a,G_s3a,Zhat_s3b, PL, PR] = ...
        int_2b(A.ZR,B.ZL,A.vars.in,vs2a,vs2b,vs3a,vs3b,dom2b);
Gs3adim = [size(G_s3a.C,1)/prod(cellfun(@numel,G_s3a.Z)),size(G_s3a.C,2)];

[Cs2a_betaa, KZhat_s2a] = int_monomial(Zhat_s2a,beta_a,dom2a);   % note that KZhat_s2a has been expanded to ensure independence of beta_a
[Cs3b_alphab, KZhat_s3b] = int_monomial(Zhat_s3b,alpha_b,dom3b);  % same for KZhat_s3b
% Cs2a_betaa, Cs3b_alphab are cell with different coefficients for
% different beta_a/alpha_b stored in same order as beta_a and alpha_b variables



% now, we have
% int(A.ZR*B.ZL', s2a,s2b,s3a,s3b) 
% = (I_p⊗PL')(I_p⊗KZhat_s2a(s2a)'*Cs2a_betaa')*int(G_s3a(s3a),s3a)*(I_q⊗Cs3b_alphab*KZhat_s3b(s3b))(I_q⊗PR)
% = (I_p⊗PL')(I_p⊗KZhat_s2a(s2a)')*(I_p⊗Cs2a_betaa')*int(G_s3a(s3a),s3a)*(I_q⊗Cs3b_alphab)*(I_q⊗KZhat_s3b(s3b))(I_q⊗PR)
% next we need to perform int(G_s3a(s3a),s3a) = (I_g1⊗barZ_3aL')*C_(gam,beta,alp)*(I_g2⊗barZ_3aR)
[C_gam_alp_beta,barZ_3aL,barZ_3aR] = int_semisep(G_s3a,beta_b,alpha_a,dom3a, Csize); 




% looking at the original composition, we have  K^C_gamma = 
% sum_{betaa,betab,alphaa,alphab) (I⊗A.ZL')CA(betaa,betab)*int((I⊗A.ZR)*(I⊗B.ZL'),s2a,s2b,s3a,s3b)*CB(alphaa,alphab)(I⊗B.ZR)
% = sum (I⊗A.ZL')CA(betaa,betab)*(I⊗int(A.ZR*B.ZL',s2a,s2b,s3a,s3b))*CB(alphaa,alphab)(I⊗B.ZR)
% = sum (I⊗A.ZL')CA(betaa,betab)
%          *(I⊗((I_p⊗PL')(I⊗KZhat_s2a(s2a)')*(I⊗Cs2a_betaa')*(I⊗barZ_3aL')*C_(gam,betab,alphaa)*(I⊗barZ_3aR)*(I⊗Cs3b_alphab)*(I⊗KZhat_s3b(s3b))(I_q⊗PR)))
%             *CB(alphaa,alphab)(I⊗B.ZR)

% split (I⊗ABCDE) as (I⊗A)(I⊗B)(I⊗C)....Then
% sum (I⊗A.ZL')CA(betaa,betab)
%       *(I⊗((I⊗KZhat_s2a(s2a)')*(I⊗Cs2a_betaa')*(I⊗barZ_3aL')*C_(gam,betab,alphaa)*(I⊗barZ_3aR)*(I⊗Cs3b_alphab)*(I⊗KZhat_s3b(s3b))))
%          *CB(alphaa,alphab)(I⊗B.ZR)
% = sum (I⊗A.ZL')CA(betaa,betab)*(I⊗PL')*(I⊗KZhat_s2a(s2a)')*(I⊗Cs2a_betaa')
%         *(I⊗barZ_3aL')*(I⊗C_(gam,betab,alphaa))*(I⊗barZ_3aR)
%          *(I⊗Cs3b_alphab)*(I⊗KZhat_s3b(s3b))))*(I⊗PR)*CB(alphaa,alphab)(I⊗B.ZR)
% since none of the monomials depend on betas and alphas, there must be a
% template permutation operation that repeats for all indices to compute
% resulting composition in the quadPoly form (I⊗C.ZL') CC[gamma] (I⊗C.ZR) 


% first we build CA and CB indexed as CA{betaa}{betab} and
% CB{alphaa}{alphab}
colsbetaa = find(idx2a); colsbetab = find(~idx2a);
if isempty(colsbetaa)
CA = A.params';
else
CA = reshape(permute(A.params, [colsbetaa,colsbetab]), 3^numel(colsbetaa), 3^numel(colsbetab)); % this is a cell 
end
colsalphaa = find(idx3a); colsalphab = find(~idx3a);
if isempty(colsalphab)
    CB = B.params;
else
CB = reshape(permute(B.params, [colsalphaa,colsalphab]), 3^numel(colsalphaa), 3^numel(colsalphab));
end

% sum_betaa (I⊗A.ZL')CA(betaa,betab)*((I⊗PL'))(I⊗KZhat_s2a(s2a)')*(I⊗Cs2a_betaa')
% = sum_betaa (I⊗ZLtemp(varLtemp)')CLtemp_betaab
PLbig = kron(eye(size(CA{1},2)/size(PL,2)),PL');

CA_PL = cell(size(CA));
for i = 1:numel(CA)
    CA_PL{i} = CA{i}*PLbig;
end

[varLtemp,ZLtemp,CLtemp_betaab]=leftShiftMonomials_SS( ...
    A.vars.out,A.ZL,CA_PL,vs2a,KZhat_s2a, ...
    cellfun(@(x) kron(eye(Gs3adim(1)),x), Cs2a_betaa, UniformOutput=false));

% [varLtemp,ZLtemp,CLtemp_betaab]=leftShiftMonomials_SS( ...
%     A.vars.out,A.ZL,CA*kron(eye(size(CA{1},2)/size(PL,2)),PL'),vs2a,KZhat_s2a, ...
%     cellfun(@(x) kron(eye(Gs3adim(1)),x), Cs2a_betaa, UniformOutput=false));
CLtemp = cell(1,size(CLtemp_betaab,2));
for j=1:size(CLtemp_betaab,2)
    CLtemp{j} = CLtemp_betaab{1,j};
    for i=2:size(CLtemp_betaab,1)
        CLtemp{j} = CLtemp{j} + CLtemp_betaab{i,j};   % adding along beta_a
    end
end


% same for right side
% sum_alphab (I⊗Cs3b_alphab)*(I⊗KZhat_s3b(s3b))))*(I⊗PR)*CB(alphaa,alphab)(I⊗B.ZR)
% = sum_alphab (I⊗CRtemp_alphaab)*(I⊗ZRtemp(varRtemp))
CBt = CB.';% IT DOES NOT TRANSPOSE MATRIX ROW/COLUMNS
CBt = cellfun(@transpose, CBt, 'UniformOutput', false);

PRbig = kron(eye(size(CBt{1},1)/size(PR,2)),PR');

CBt_PR = cell(size(CBt));
for i = 1:numel(CBt)
    CBt_PR{i} = CBt{i}*PRbig;
end


cell_transpose_1 = cellfun(@(x) kron(eye(Gs3adim(2)),x), Cs3b_alphab, UniformOutput=false).';
cell_transpose_1 = cellfun(@transpose, cell_transpose_1, 'UniformOutput', false);
[varRtemp,ZRtemp,CRtemp_alphaba]=leftShiftMonomials_SS( ...
    B.vars.in,B.ZR,CBt_PR,vs3b,KZhat_s3b, cell_transpose_1);
% [varRtemp,ZRtemp,CRtemp_alphaba]=leftShiftMonomials_SS( ...
%     B.vars.in,B.ZR,CB.'*kron(eye(size(CB{1},1)/size(PR,2)),PR'),vs3b,KZhat_s3b, ...
%     cellfun(@(x) kron(eye(Gs3adim(2)),x), Cs3b_alphab, UniformOutput=false).');
% [varRtemp,ZRtemp,CRtemp_alphaba]=leftShiftMonomials_SS( ...
%     B.vars.in,B.ZR,CB.',vs3b,KZhat_s3b, kron(eye(B.dim(1)),Cs3b_alphab)');  % verify that CB.' not just transposes cell rows/columns but also matrix rows/columns
% note above that since we used leftShiftMonomials function by transposing
% everything the output is not CRtemp(alpha_a,alpha_b) but CRtemp(alpha_b,alpha_a)
% so to sum along alpha_b we will again sum along first index and not 2nd.
CRtemp = cell(size(CRtemp_alphaba,2),1);
for j=1:size(CRtemp_alphaba,2)
    CRtemp{j} = CRtemp_alphaba{1,j};
    for i=2:size(CRtemp_alphaba,1)
        CRtemp{j} = CRtemp{j}+CRtemp_alphaba{i,j};  % adding along alpha_b
    end
end
% now CLtemp is function of beta_b, and CRtemp is function of alpha_a

% now we have
% sum (I⊗A.ZL')CA(betaa,betab)*(I⊗KZhat_s2a(s2a)')*(I⊗Cs2a_betaa')
%         *(I⊗barZ_3aL')*(I⊗C_(gam,betab,alphaa))*(I⊗barZ_3aR)
%          *(I⊗Cs3b_alphab)*(I⊗KZhat_s3b(s3b))))*CB(alphaa,alphab)(I⊗B.ZR)
% = sum (I⊗ZLtemp(varLtemp)')*CLtemp(betab)*
%         *(I⊗barZ_3aL')*(I⊗C_(gam,betab,alphaa))*(I⊗barZ_3aR)
%          *CRtemp(alphaa)*(I⊗ZRtemp(varRtemp))

% now shift barZ_3aL too
% (I⊗ZLtemp(varLtemp)')*CLtemp(betab)*(I⊗barZ_3aL')
[Cvarsout,CZL,CLtemp_betab]=leftShiftMonomials_SS( ...
    varLtemp,ZLtemp,CLtemp,vs3a,barZ_3aL,{eye(Gs3adim(1)*prod(cellfun(@numel,barZ_3aL)))});

% likewise for right side
[Cvarsin,CZR,CRtemp_alphaa]=leftShiftMonomials_SS( ...
    varRtemp,ZRtemp,CRtemp.',vs3a,barZ_3aR,{eye(Gs3adim(2)*prod(cellfun(@numel,barZ_3aR)))});


for i=1:length(Cparams)
    Cparams{i} = 0;
    for j=1:size(beta_b,1)
        for k=1:size(alpha_a,1)
            Cparams{i} = Cparams{i} + CLtemp_betab{j}*C_gam_alp_beta{i,j,k}*CRtemp_alphaa{k}';
        end
    end
end

Cvars = struct('in', {Cvarsin}, 'out', {Cvarsout});
C = sopvar(Cparams,Cvars,CZR,CZL,Cdom,Cdims);
end

function [C_gam_alp_beta,ZL,ZR] = int_semisep(G,idxbeta,idxalpha,lims,Csize)
% suppose G(s3a) is size p x q; then G.CG size is p*len(G.ZG) x q
% this performs
% int(G(s3a), s3a, lims) = (I_p⊗int(G.Z'))*G.CG
%   = sum_{g,a,b} (I⊗ZL(s3a)') C_{g,a,b} (I ⊗ ZR(s3a'))
% G is structure storing ZG and CG where G(s3a)=(I_p⊗ZG(s3a)')*CG

% ZG(s3a) = ZG(s3a_1)\otimes ... ZG(s3a_i)...
% int(ZG(s3a),indicator(i)) = 
%   ZG(s3a_1)\otimes ...\otimes Ci*(ZLunique(s3a_i)\otimes ZLunique(t3a_i))...

ns3a = size(idxbeta,2);
if ns3a==0
    C_gam_alp_beta = {[1]};
    ZL = [];
    ZR = [];
    return
end
a = lims(:,1); b = lims(:,2);

CG = G.C; ZG = G.Z;
p = size(CG,1)/prod(cellfun(@numel,ZG));
q = size(CG,2);

ZL = ZG; ZLunique = ZL; 
ZLuniqueIdx = ZLunique;
for i=1:numel(ZG)
    mons = [0;ZG{i};ZG{i}+1];
    ZL{i} = mons;   % this has monomials [si^0; si^ZG{i}; si^{ZG{i}+1}]; note there may be repetition
    [ZLunique{i},~,ZLuniqueIdx{i}] = unique(mons);
end

C_gam_alp_beta = repmat({0},3^ns3a,size(idxbeta,1),size(idxalpha,1));

gamIdx = cell(1,ns3a);
[gamIdx{:}] = ind2sub(Csize,1:3^ns3a);
if ns3a==0
    gamIdx = [];
elseif ns3a==1
    gamIdx = cell2mat(gamIdx).';
else
    gamIdx = cell2mat(gamIdx(:)).';
end

for i=1:size(C_gam_alp_beta,2)
    idxB = idxbeta(i,:)';
    for j=1:size(C_gam_alp_beta,3)
        idxA = idxalpha(j,:)';
        for k=1:size(C_gam_alp_beta,1)
            gam = gamIdx(k,:);
            C=1;          
            for l=1:ns3a
            E = ZG{l}; nE = length(E);  % E is exponent in l-th variable
            Eint = ZL{l}; nEint = length(Eint); % Eint is exponent in l-th variable after integration+multiplier
            Ci = zeros(nE,nEint^2);   % basically int(ZG{i}*I_{gam,alp,bet}) = Ci*(ZL{2*i+1}(si)\otimes ZL{2*i+2}(ti)) 
            key = 100*gam(l)+10*idxB(l)+idxA(l);
            switch key
                case 111 % G(s_3a_i)
                    Ci(1:nE,nEint+1:nEint:nE*nEint+1) = eye(nE); 
                case 221  % G(t_3a_i);
                    Ci(1:nE,2:nE+1) = eye(nE);
                case 331  % G(t_3a_i);
                    Ci(1:nE,2:nE+1) = eye(nE);
                case 212  % G(s_3a_i)
                    Ci(1:nE,nEint+1:nEint:nE*nEint+1) = eye(nE); 
                case 222% int_{t_3a_i}^{s_3a_i} eta_3a_i^E deta_3a_i = s_3a_i^(E+1)/(E+1) - t_3a_i^(E+1)/(E+1);                                    
                    Ci(1:nE,nE+2:2*nE+1) = -diag(1./(E+1));  % -t^(E+1)/(E+1)
                    Ci(1:nE,(nE+1)*nEint+1:nEint:(2*nE)*nEint+1) = diag(1./(E+1)); % +s^(E+1)/(E+1)
                % case 222  % int_{s_3a_i}^{t_3a,i} eta_3a_i^E deta_3a_i = t_3a_i^{E+1}/(E+1) -  s_3a_i^(E+1)/(E+1); % THIS may need to change
                %     Ci(1:nE,nE+2:2*nE+1) = diag(1./(E+1));
                %     Ci(1:nE,(nE+1)*nEint+1:nEint:(2*nE)*nEint+1) = -diag(1./(E+1));
                case 232  % int_{s_3a_i}^{b(i)} eta_3a_i^E deta_3a_i = b(i)^{E+1}/(E+1) -  s_3a_i^(E+1)/(E+1);
                    Ci(1:nE,1) = b(l).^(E+1)./(E+1);
                    Ci(1:nE,(nE+1)*nEint+1:nEint:(2*nE)*nEint+1) = -diag(1./(E+1));
                case 332  % int_{t_3a_i}^{b(i)} eta_3a_i^E deta_3a_i = b(i)^{E+1}/(E+1) -  t_3a_i^(E+1)/(E+1);
                    Ci(1:nE,1) = b(l).^(E+1)./(E+1);
                    Ci(1:nE,nE+2:2*nE+1) = -diag(1./(E+1));
                case 313  % G(s_3a_i)
                    Ci(1:nE,nEint+1:nEint:nE*nEint+1) = eye(nE); 
                case 223  % int_{a(i)}^{t_3a_i} eta_3a_i^E deta_3a_i = t_3a_i^(E+1)/(E+1) - a(i)^{E+1}/(E+1);
                    Ci(1:nE,1) = -a(l).^(E+1)./(E+1);
                    Ci(1:nE,nE+2:2*nE+1) = diag(1./(E+1));
                case 323  % int_{a(i)}^{s_3a_i} eta_3a_i^E deta_3a_i = s_3a_i^(E+1)/(E+1) - a(i)^{E+1}/(E+1);
                    Ci(1:nE,1) = -a(l).^(E+1)./(E+1);
                    Ci(1:nE,(nE+1)*nEint+1:nEint:(2*nE)*nEint+1) = diag(1./(E+1));
                case 333  % int_{s_3a_i}^{t_3a_i} eta_3a_i^E deta_3a_i= t_3a_i^(E+1)/(E+1) - s_3a_i^(E+1)/(E+1);
                    Ci(1:nE,(nE+1)*nEint+1:nEint:(2*nE)*nEint+1) = -diag(1./(E+1)); % -s^(E+1)/(E+1)
                    Ci(1:nE,nE+2:2*nE+1) = diag(1./(E+1));  % +t^(E+1)/(E+1)
                % case 333  % int_{t_3a_i}^{s_3a,i} eta_3a_i^E deta_3a_i = s_3a_i^{E+1}/(E+1) -  t_3a_i^(E+1)/(E+1);  % THIS may need to change
                %     Ci(1:nE,(nE+1)*nEint+1:nEint:(2*nE)*nEint+1) = diag(1./(E+1));
                %     Ci(1:nE,nE+2:2*nE+1) = -diag(1./(E+1));
                otherwise % cases 211,311,121,321,131,231,112,312,122,322,132,113,213,123,133,233  Ci=0
                    lenZG = prod(cellfun(@numel, ZLunique));
                    C = zeros(nE,lenZG^2);
                    break;  % if one of the Ci is zero, kron is zero so exit loop
            end
            Rs = zeros(nEint,length(ZLunique{l}));  % Rs*ZLunique{l} = ZL
            Rs((1:nEint).'+(ZLuniqueIdx{l}-1)*nEint) = 1;
            R = kron(Rs,Rs); % (Rs\otimes Rs) (ZLunique(si)\otimes ZLunique(ti)) = (ZL(si)\otimes ZL(ti))
            Ci = Ci*R; % Ci*(Rs\otimes Rs)(ZLunique(si)\otimes ZLunique(ti)) = Ci(ZL(si)\otimes ZL(ti))
            C = kron(C,Ci);
            end
            % convert (I_p⊗Zint(s,s_dum)'*C')*CG =(I_p⊗Zint(s,s_dum)')(I_p⊗C')*CG= (I_p⊗ZintNew(s)')Cnew(I_q⊗ZintNew(s_dum))
            % where Zint(s,s_dum) = (Zint(s1)\otimes Zint(s1_dum))\otimes... \otimes (Zint(si)\otimes Zint(si_dum))
            Cnew = rearrangeCoef(kron(eye(p),C'),CG,p,q, ZLunique);
            C_gam_alp_beta{k,i,j} = C_gam_alp_beta{k,i,j}+Cnew;
        end
    end
end
ZL = ZLunique; ZR = ZL;
end

function Cout = rearrangeCoef(Cmul,C,p,q,Zint)
% convert (I_p⊗Zintst(s,s_dum)'*C1')*C =(I_p⊗Zintst(s,s_dum)')(I_p⊗C1')*C=
% (I_p⊗Zint(s)')Cnew(I_q⊗Zint(s_dum))  
% Cmul = (I_p⊗C1') and q = size(C,2)

n = length(Zint);  % number of s3a
ns = cellfun(@numel,Zint);  % size of monomials in s3a
Cout = full(Cmul*C); 
[nrows, ~] = size(Cout);
nz = nrows/p;  % nz must be ns^2  since we have Z(s3a)\otimes Z(t3a) mixed
K = reshape(1:nz,repelem(fliplr(ns),2));  % first generate lengths for each variable
% K corresponds to columns (len(tn), len(sn), ...., len(t1), len(s1))
K = permute(K,[1:2:2*n,2:2:2*n]);  % separate s and t variables
% now, K corresponds to columns (len(tn), ...,len(t1), len(sn),...len(s1))
Cout = reshape(Cout,nz,p,q);
Cout = Cout(K(:),:,:);  % now Cout columns are ordered as in K columns
Cout = reshape(Cout,prod(ns),prod(ns),p,q);  % now pull out ti and si into two directions Cout(ti,si,p,q)
Cout = permute(Cout,[2,3,1,4]);  % change order to Cout(si,p,ti,q)
Cout = reshape(Cout,p*prod(ns),q*prod(ns)); % club si*p and ti*q
end


function [Cidx, Zint] = int_monomial(Z,idxAll,lims)
% this performs 
% no integration if idx(i) = 1
% int(si^Z{i}, si, a(i), si) if idx(i) = 2
% int(si^Z{i}, si, si, b(i)) if idx(i) = 3
% returns C*Zint = int(Z,idx,lims)
if isempty(Z)
    Cidx{1} = 1;
    Zint = {[]};
    return;
end

Cidx = cell(1,size(idxAll,1));  % storing C_betaa for all beta_a (or C_alphab for all alpha_b)

Zint = cell(size(Z));
a = lims(:,1); b = lims(:,2);

for k = 1:length(Cidx)
idx = idxAll(k,:)';   % extract a specific beta_a or alpha_b
C = 1;                % this is C corresponding to that beta_a/alpha_b int(Z,idx) = C*Zint
for i=1:length(idx)
    E = Z{i};
    nE = length(E);
    switch idx(i)
        case 1
            % no integration along si
            Ci = eye(nE);
            Zint{i} = Z{i};
        case 2
            % int_a^si si^E dsi = si^(E+1)/(E+1) - a^(E+1)/(E+1);
            Ci = zeros(nE,nE+1);
            Ci(:,1) = -a(i).^(E+1)./(E+1);  % constant of integration
            Ci(:,2:end) = diag(1./(E+1));   % si^(E+1)/(E+1) terms
            Zint{i} = [0; Z{i}+1];
        case 3
            % int_si^b si^E dsi = b^(E+1)/(E+1)-si^(E+1)/(E+1);
            Ci = zeros(nE,nE+1);
            Ci(:,1) = b(i).^(E+1)./(E+1);  % constant of integration
            Ci(:,2:end) = -diag(1./(E+1));   % -si^(E+1)/(E+1) terms
            Zint{i} = [0; Z{i}+1];
    end
    C = kron(C,Ci);
end
Cidx{k} = C;
end
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [Z2aout, G3aout, Z3bout, PL, PR] = int_2b(ZL, ZR, Zvar, s2a,s2b,s3a,s3b,lims)
% This performs the factorization
% int(ZL*ZR',s2b,0,1) = PL'*(Im\otimes Z2a') G(s3a) (In\otimes Z3b)*PR
% where m = length(kron(ZL)) and n = length(kron(ZR))
%
% note that PL and PR are permutation matrices which account for the fact
% that Zvar may not be ordered as [s2a,s2b,s3a,s3b]

% note that
% int_0^1 (ZL2a⊗ZL2b⊗ZL3a⊗ZL3b)(ZR2a⊗ZR2b⊗ZR3a⊗ZR3b)' d2b
% = int_0^1 (ZL2a*ZR2a')⊗(ZL2b*ZR2b')⊗(ZL3a*ZR3a')⊗(ZL3b*ZR3b') d2b
% = (ZL2a*ZR2a')⊗C2b⊗(ZL3a*ZR3a')⊗(ZL3b*ZR3b')
% for some constant matrix C2b

% Zvar must appear in both for composition to be well-defined, but may not
% be ordered as [s2a,s2b,s3a,s3b]

Zvar_ordered = [s2a,s2b,s3a,s3b];

% first compute the permutation matrices so that
%       Z_ordered = P*Z_original
% hence
%       Z_original = P'*Z_ordered
dimsL = cellfun(@numel,ZL);
dimsR = cellfun(@numel,ZR);

PL = permute_kron_basis(Zvar,Zvar_ordered,dimsL);
PR = permute_kron_basis(Zvar,Zvar_ordered,dimsR);

% first separate ZL and ZR into Z2a,Z2b,Z3a,Z3b
[~,loc] = ismember(s2a,Zvar);
ZL2a = ZL(loc);
ZR2a = ZR(loc);
[~,loc] = ismember(s2b,Zvar);
ZL2b = ZL(loc);
ZR2b = ZR(loc);
[~,loc] = ismember(s3a,Zvar);
ZL3a = ZL(loc);
ZR3a = ZR(loc);
[~,loc] = ismember(s3b,Zvar);
ZL3b = ZL(loc);
ZR3b = ZR(loc);

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

% condense Z2a matrix into unique monomials
Z2anew = cell(1,length(s2a));
idx2a = cell(1,length(s2a));
for i=1:length(s2a)
    E = Z2a{i};
    [Z2anew{i},~,idx] = unique(E(:));
    idx2a{i} = reshape(idx,size(E));
end

% integrate Z2b part over [a,b]
% int_a^b Z2b d2b = 
%       kron_i (b(i)^(Z2b{i}+1)-a(i)^(Z2b{i}+1))/(Z2b{i}+1)
C2b = cell(1,length(s2b));
a = lims(:,1); b = lims(:,2);
for i=1:length(s2b)
    E = Z2b{i};
    C2b{i} = (b(i).^(E+1)-a(i).^(E+1))./(E+1);
end

% condense Z3a matrix into unique monomials
Z3anew = cell(1,length(s3a));
idx3a = cell(1,length(s3a));
for i=1:length(s3a)
    E = Z3a{i};
    [Z3anew{i},~,idx] = unique(E(:));
    idx3a{i} = reshape(idx,size(E));
end

% condense Z3b matrix into unique monomials
Z3bnew = cell(1,length(s3b));
idx3b = cell(1,length(s3b));
for i=1:length(s3b)
    E = Z3b{i};
    [Z3bnew{i},~,idx] = unique(E(:));
    idx3b{i} = reshape(idx,size(E));
end

% ordered monomial dimensions
dimL2a = cellfun(@numel,ZL2a);
dimL2b = cellfun(@numel,ZL2b);
dimL3a = cellfun(@numel,ZL3a);
dimL3b = cellfun(@numel,ZL3b);

dimR2a = cellfun(@numel,ZR2a);
dimR2b = cellfun(@numel,ZR2b);
dimR3a = cellfun(@numel,ZR3a);
dimR3b = cellfun(@numel,ZR3b);

dimsLord = [dimL2a,dimL2b,dimL3a,dimL3b];
dimsRord = [dimR2a,dimR2b,dimR3a,dimR3b];

m = prod_empty(dimsLord);
n = prod_empty(dimsRord);

% sizes of the condensed monomial bases
nz2a = prod_empty(cellfun(@numel,Z2anew));
nz3a = prod_empty(cellfun(@numel,Z3anew));
nz3b = prod_empty(cellfun(@numel,Z3bnew));

% generate ordered Kronecker multi-indices
IL = kron_indices(dimsLord);
IR = kron_indices(dimsRord);

% index offsets in the ordered variables
n2a = length(s2a);
n2b = length(s2b);
n3a = length(s3a);
n3b = length(s3b);

L2a = 1:n2a;
L2b = n2a + (1:n2b);
L3a = n2a+n2b + (1:n3a);
L3b = n2a+n2b+n3a + (1:n3b);

R2a = 1:n2a;
R2b = n2a + (1:n2b);
R3a = n2a+n2b + (1:n3a);
R3b = n2a+n2b+n3a + (1:n3b);

% build C3a directly so that
%
% int(ZL*ZR',s2b)
%   = PL'*(Im⊗Z2anew')*(Im*nz2a⊗Z3anew')*C3a*(In⊗Z3bnew)*PR
%
% equivalently, G3aout.C = C3a and G3aout.Z = Z3anew
II = zeros(m*n,1);
JJ = zeros(m*n,1);
VV = zeros(m*n,1);

cnt = 0;
for i=1:m
    for j=1:n

        % coefficient from integrating the s2b variables
        c = 1;
        for k=1:n2b
            c = c*C2b{k}(IL(i,L2b(k)),IR(j,R2b(k)));
        end

        % index in Z2anew
        idx = zeros(1,n2a);
        for k=1:n2a
            idx(k) = idx2a{k}(IL(i,L2a(k)),IR(j,R2a(k)));
        end
        k2a = kron_sub2ind(cellfun(@numel,Z2anew),idx);

        % index in Z3anew
        idx = zeros(1,n3a);
        for k=1:n3a
            idx(k) = idx3a{k}(IL(i,L3a(k)),IR(j,R3a(k)));
        end
        k3a = kron_sub2ind(cellfun(@numel,Z3anew),idx);

        % index in Z3bnew
        idx = zeros(1,n3b);
        for k=1:n3b
            idx(k) = idx3b{k}(IL(i,L3b(k)),IR(j,R3b(k)));
        end
        k3b = kron_sub2ind(cellfun(@numel,Z3bnew),idx);

        % row and column in C3a
        row = k3a + ((k2a + (i-1)*nz2a)-1)*nz3a;
        col = k3b + (j-1)*nz3b;

        cnt = cnt + 1;
        II(cnt) = row;
        JJ(cnt) = col;
        VV(cnt) = c;
    end
end

C3a = sparse(II(1:cnt),JJ(1:cnt),VV(1:cnt),m*nz2a*nz3a,n*nz3b);

% finally, we have 
% int(ZL*ZR',s2b)
%   = PL'*(Im⊗Z2anew')*G3aout(s3a)*(In⊗Z3bnew)*PR
Z2aout = Z2anew;
Z3bout = Z3bnew;
G3aout = struct('C', C3a, 'Z', {Z3anew});
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function P = permute_kron_basis(Zvar,Zvar_ordered,dims)
% This computes the permutation matrix P such that
%       Z_ordered = P*Z_original
% where Z_original is kron ordered according to Zvar and Z_ordered is kron
% ordered according to Zvar_ordered.

if isempty(dims)
    P = 1;
    return
end

[~,idx] = ismember(Zvar_ordered,Zvar);

if any(idx==0)
    error('variables used in int_2b are not contained in Zvar');
end

N = prod_empty(dims);
dims_ordered = dims(idx);

if length(dims)==1
    p = (1:N).';
else
    % Original Kronecker multi-indices, in kron ordering
    r = cell(1,length(dims));
    [r{:}] = ind2sub(fliplr(dims),1:N);
    I = flipud(vertcat(r{:})).';   % rows are [i1 i2 ... in]

    % Same monomial, but indexed in the permuted Kronecker basis
    Ip = I(:,idx);

    % Convert permuted multi-indices back to linear kron indices
    rp = num2cell(fliplr(Ip),1);
    p = sub2ind(fliplr(dims_ordered),rp{:}).';
end

% Permutation matrix satisfying P*Z_original = Z_ordered
P = sparse(p,1:N,1,N,N);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function I = kron_indices(dims)
% This returns the multi-indices associated with kron ordering.

if isempty(dims)
    I = zeros(1,0);
    return
end

N = prod_empty(dims);

if length(dims)==1
    I = (1:N).';
else
    r = cell(1,length(dims));
    [r{:}] = ind2sub(fliplr(dims),1:N);
    I = flipud(vertcat(r{:})).';
end
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function idx = kron_sub2ind(dims,sub)
% This converts a multi-index into a linear index in kron ordering.

if isempty(dims)
    idx = 1;
    return
end

if length(dims)==1
    idx = sub;
else
    sub = num2cell(fliplr(sub));
    idx = sub2ind(fliplr(dims),sub{:});
end
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function p = prod_empty(x)
% This treats empty products as one.

if isempty(x)
    p = 1;
else
    p = prod(x);
end
end