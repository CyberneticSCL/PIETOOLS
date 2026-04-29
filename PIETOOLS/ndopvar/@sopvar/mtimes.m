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


Cdom = struct('in',B.dom.in,'out',A.dom.out);
Cdims = [A.dims(1),B.dims(2)];

n3 = numel(vs3a);

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
% convert linear index to multi-index for B parameters
for i=1:numel(B.params)
    tmp = cell(1,nalpha);
    [tmp{:}] = ind2sub(Bsize(:),i);
    alphaIdx(i,:) = [tmp{:}];
end
% convert linear index to multi-index for A parameters
for i=1:numel(A.params)
    tmp = cell(1,nbeta);
    [tmp{:}] = ind2sub(Asize(:),i);
    betaIdx(i,:) = [tmp{:}];
end

% split alpha and beta as [alpha_a, alpha_b] and [beta_a, beta_b]
[idx2a,~] = ismember(A.vars_S3,vs2a);
beta_a = betaIdx(:,idx2a);
beta_b = betaIdx(:,~idx2a);
[idx3a,~] = ismember(B.vars_S3,vs3a);
alpha_a = alphaIdx(:,idx3a);
alpha_b = alphaIdx(:,~idx3a);

% REMEMBER beta_a, beta_b, alpha_a, alpha_b no longer have unique rows
% above. Need to revisit functions below to ensure this is handled
% appropriately!!!

[idx2b,~] = ismember(B.vars_S2,vs2b);
dom2b = B.dom_2(idx2b,:);
dom2a = A.dom_3(idx2a,:);
dom3a = B.dom_3(idx3a,:);
dom3b = B.dom_3(~idx3a,:);

% do int(A.ZR*B.ZL', s2b) = (I⊗Zhat_s2a')*G_s3a(s3a)*(I⊗Zhat_s3b)
[Zhat_s2a,G_s3a,Zhat_s3b] = ...
        int_2b(A.ZR,B.ZL,A.vars.in,B.vars.out,vs2a,vs2b,vs3a,vs3b,dom2b);

[Cs2a_betaa, KZhat_s2a] = int_monomial(Zhat_s2a,beta_a,dom2a);   % note that KZhat_s2a has been expanded to ensure independence of beta_a
[Cs3b_alphab, KZhat_s3b] = int_monomial(Zhat_s3b,alpha_b,dom3b);  % same for KZhat_s3b
% Cs2a_betaa, Cs3b_alphab are cell with different coefficients for
% different beta_a/alpha_b stored in same order as beta_a and alpha_b variables



% now, we have
% int(A.ZR*B.ZL', s2a,s2b,s3a,s3b) 
% = (I⊗KZhat_s2a(s2a)'*Cs2a_betaa')*int(G_s3a(s3a),s3a)*(I⊗Cs3b_alphab*KZhat_s3b(s3b))
% = (I⊗KZhat_s2a(s2a)')*(I⊗Cs2a_betaa')*int(G_s3a(s3a),s3a)*(I⊗Cs3b_alphab)*(I⊗KZhat_s3b(s3b))
% next we need to perform int(G_s3a(s3a),s3a) = (I⊗barZ_3aL')*C_(gam,beta,alp)*(I⊗barZ_3aR)
[C_gam_alp_beta,barZ_3aL,barZ_3aR] = int_semisep(G_s3a,beta_b,alpha_a,dom3a); 


% looking at the original composition, we have  K^C_gamma = 
% sum_{betaa,betab,alphaa,alphab) (I⊗A.ZL')CA(betaa,betab)*int((I⊗A.ZR)*(I⊗B.ZL'),s2a,s2b,s3a,s3b)*CB(alphaa,alphab)(I⊗B.ZR)
% = sum (I⊗A.ZL')CA(betaa,betab)*(I⊗int(A.ZR*B.ZL',s2a,s2b,s3a,s3b))*CB(alphaa,alphab)(I⊗B.ZR)
% = sum (I⊗A.ZL')CA(betaa,betab)
%          *(I⊗((I⊗KZhat_s2a(s2a)')*(I⊗Cs2a_betaa')*(I⊗barZ_3aL')*C_(gam,betab,alphaa)*(I⊗barZ_3aR)*(I⊗Cs3b_alphab)*(I⊗KZhat_s3b(s3b))))
%             *CB(alphaa,alphab)(I⊗B.ZR)

% split (I⊗ABCDE) as (I⊗A)(I⊗B)(I⊗C)....Then
% sum (I⊗A.ZL')CA(betaa,betab)
%       *(I⊗((I⊗KZhat_s2a(s2a)')*(I⊗Cs2a_betaa')*(I⊗barZ_3aL')*C_(gam,betab,alphaa)*(I⊗barZ_3aR)*(I⊗Cs3b_alphab)*(I⊗KZhat_s3b(s3b))))
%          *CB(alphaa,alphab)(I⊗B.ZR)
% = sum (I⊗A.ZL')CA(betaa,betab)*(I⊗KZhat_s2a(s2a)')*(I⊗Cs2a_betaa')
%         *(I⊗barZ_3aL')*(I⊗C_(gam,betab,alphaa))*(I⊗barZ_3aR)
%          *(I⊗Cs3b_alphab)*(I⊗KZhat_s3b(s3b))))*CB(alphaa,alphab)(I⊗B.ZR)
% since none of the monomials depend on betas and alphas, there must be a
% template permutation operation that repeats for all indices to compute
% resulting composition in the quadPoly form (I⊗C.ZL') CC[gamma] (I⊗C.ZR) 


% first we build CA and CB indexed as CA{betaa}{betab} and
% CB{alphaa}{alphab}
colsbetaa = find(idx2a); colsbetab = find(~idx2a);
if isempty(colsbetaa)
CA = A.params;
else
CA = reshape(permute(A.params, [colsbetaa,colsbetab]), 3^numel(colsbetaa), 3^numel(colsbetab)); % this is a cell 
end
colsalphaa = find(idx3a); colsalphab = find(~idx3a);
if isempty(colsalphab)
    CB = B.params;
else
CB = reshape(permute(B.params, [colsalphaa,colsalphab]), 3^numel(colsalphaa), 3^numel(colsalphab));
end

% sum_betaa (I⊗A.ZL')CA(betaa,betab)*(I⊗KZhat_s2a(s2a)')*(I⊗Cs2a_betaa')
% = sum_betaa (I⊗ZLtemp(varLtemp)')CLtemp_betaab
[varLtemp,ZLtemp,CLtemp_betaab]=leftShiftMonomials_SS( ...
    A.vars.out,A.ZL,CA,vs2a,KZhat_s2a, ...
    cellfun(@(x) kron(eye(A.dims(2)),x), Cs2a_betaa, UniformOutput=false));
CLtemp = cell(1,size(CLtemp_betaab,2));
for j=size(CLtemp_betaab,2)
    CLtemp{j} = CLtemp_betaab{1,j};
    for i=2:size(CLtemp_betaab,1)
        CLtemp{j} = CLtemp{j} + CLtemp_betaab{i,j};   % adding along beta_a
    end
end


% same for right side
% sum_alphab (I⊗Cs3b_alphab)*(I⊗KZhat_s3b(s3b))))*CB(alphaa,alphab)(I⊗B.ZR)
% = sum_alphab (I⊗CRtemp_alphaab)*(I⊗ZRtemp(varRtemp))
[varRtemp,ZRtemp,CRtemp_alphaba]=leftShiftMonomials_SS( ...
    B.vars.in,B.ZR,CB.',vs3b,KZhat_s3b, ...
    cellfun(@(x) kron(eye(B.dims(1)),x), Cs3b_alphab, UniformOutput=false).');
% [varRtemp,ZRtemp,CRtemp_alphaba]=leftShiftMonomials_SS( ...
%     B.vars.in,B.ZR,CB.',vs3b,KZhat_s3b, kron(eye(B.dim(1)),Cs3b_alphab)');  % verify that CB.' not just transposes cell rows/columns but also matrix rows/columns
% note above that since we used leftShiftMonomials function by transposing
% everything the output is not CRtemp(alpha_a,alpha_b) but CRtemp(alpha_b,alpha_a)
% so to sum along alpha_b we will again sum along first index and not 2nd.
CRtemp = cell(1,size(CRtemp_alphaba,2));
for j=size(CRtemp_alphaba,2)
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
    varLtemp,ZLtemp,CLtemp,vs3a,barZ_3aL,eye(A.dim(2)*prod(cellfun(@numel,barZ_3aL))));

% likewise for right side
[Cvarsin,CZR,CRtemp_alphaa]=leftShiftMonomials_SS( ...
    varRtemp,ZRtemp,CRtemp.',vs3a,barZ_3aR,eye(B.dim(1)*prod(cellfun(@numel,barZ_3aR))));


for i=1:length(Cparams)
    Cparams{i} = 0;
    for j=1:size(beta_b,1)
        for k=1:size(alpha_a,1)
            Cparams{i} = Cparams{i} + CLtemp_betab{j}*C_gam_alp_beta{i,j,k}*CRtemp_alphaa{k}';
        end
    end
end

Cvars = struct('in', Cvarsin, 'out', Cvarsout);
C = sopvar(Cparams,Cvars,CZR,CZL,Cdom,Cdims);
end

function [C_gam_alp_beta,ZL,ZR] = int_semisep(G,idxbeta,idxalpha,lims)
% this performs
% int(G(s3a), s3a, lims) = sum_{g,a,b} (I⊗ZL(s3a)') C_{g,a,b} (I ⊗ ZR(s3a'))
% G is structure storing ZG and CG where G(s3a)=(I_q⊗ZG(s3a)')*CG

CG = G.C; ZG = G.Z;
q = size(CG,1)/prod(cellfun(@numel,ZG));

Zint = cellfun(@(x) [0; x+1], ZG, UniformOutput=false);
ZL = Zint; ZR = Zint;

n = size(idxbeta,2);
a = lims(:,1); b = lims(:,2);
C_gam_alp_beta = cell(3^n,size(idxbeta,1),size(idxalpha,1));
if n==0 
    Csize=[1,1]; 
elseif n==1
    Csize=[1,3]; 
else 
    Csize = repmat(3,1,n); 
end

for i=1:size(C_gam_alp_beta,2)
    idxB = idxbeta(i,:)';
    for j=1:size(C_gam_alp_beta,3)
        idxA = idxalpha(j,:)';
        gam = mapBetaAlpha2Gamma(idxB,idxA);
        C=1;
        shiftMonomial = zeros(size(idxB));
        for k=1:length(idxA)
            E = ZG{k}; nE = length(E);
            Ci = zeros(2*nE+1,2*nE+1);
            key = 100*gam(k)+10*idxB(k)+idxA(k);
            breakflag =0;
            switch key
                case 111 % G(s_3a_i)
                    Ci(1:nE,1:nE) = eye(nE);
                case 221  % G(t_3a_i);
                    Ci(1:nE,[1,nE+2:2*nE]) = eye(nE);
                case 331  % G(t_3a_i);
                    Ci(1:nE,[1,nE+2:2*nE]) = eye(nE);
                case 212  % G(s_3a_i)
                    Ci(1:nE,1:nE) = eye(nE);
                case 222  % int_{s_3a_i}^{t_3a,i} eta_3a_i^E deta_3a_i = t_3a_i^{E+1}/(E+1) -  s_3a_i^(E+1)/(E+1); % THIS may need to change
                    Ci(1:nE,nE+2:2*nE+1) = diag(1./(E+1));
                    Ci(1:nE,2:nE+1) = -diag(1./(E+1));
                case 232  % int_{s_3a_i}^{b(i)} eta_3a_i^E deta_3a_i = b(i)^{E+1}/(E+1) -  s_3a_i^(E+1)/(E+1);
                    Ci(1:nE,1) = b(k).^(E+1)./(E+1);
                    Ci(1:nE,2:nE+1) = -diag(1./(E+1));
                case 332  % int_{t_3a_i}^{b(i)} eta_3a_i^E deta_3a_i = b(i)^{E+1}/(E+1) -  t_3a_i^(E+1)/(E+1);
                    Ci(1:nE,1) = b(k).^(E+1)./(E+1);
                    Ci(1:nE,nE+2:2*nE+1) = -diag(1./(E+1));
                case 313  % G(s_3a_i)
                    Ci(1:nE,1:nE) = eye(nE);
                case 223  % int_{a(i)}^{t_3a_i} eta_3a_i^E deta_3a_i = t_3a_i^(E+1)/(E+1) - a(i)^{E+1}/(E+1);
                    Ci(1:nE,1) = -a(k).^(E+1)./(E+1);
                    Ci(1:nE,nE+2:2*nE+1) = diag(1./(E+1));
                case 323  % int_{a(i)}^{s_3a_i} eta_3a_i^E deta_3a_i = s_3a_i^(E+1)/(E+1) - a(i)^{E+1}/(E+1);
                    Ci(1:nE,1) = -a(k).^(E+1)./(E+1);
                    Ci(1:nE,2:nE+1) = diag(1./(E+1));
                case 333  % int_{t_3a_i}^{s_3a,i} eta_3a_i^E deta_3a_i = s_3a_i^{E+1}/(E+1) -  t_3a_i^(E+1)/(E+1);  % THIS may need to change
                    Ci(1:nE,2:nE+1) = diag(1./(E+1));
                    Ci(1:nE,nE+2:2*nE+1) = -diag(1./(E+1));
                otherwise % cases 211,311,121,321,131,231,112,312,122,322,132,113,213,123,133,233  Ci=0
                    breakflag = 1;
            end
            if breakflag
                lenZG = prod(cellfun(@(x) length(x)+1, ZG));
                C = zeros(2*lenZG-1,2*lenZG-1);
                break;  % if one of the Ci is zero, kron is zero so exit loop
            else
                C = kron(C,Ci);
            end
        end
        Cnew = rearrangeCoef(kron(eye(q),C'),blkdiag(0,CG,CG),q);  % convert (I_q⊗Zint(s,s_dum)'*C')*CG =(I_q⊗Zint(s,s_dum)')(I_q⊗C')*CG= (I_q⊗Zint(s)')Cnew(I⊗Zint(s_dum))
    end
    
    gam = num2cell(gam);
    for gidx = 1:numel(gam)
        gamlin = sub2ind(Csize,gam{gidx,:});
        C_gam_alp_beta{gamlin,i,j} = Cnew;
    end
end
end

function Cout = rearrangeCoef(Cmul,C,q)
% convert (I_q⊗Zint(s,s_dum)'*C1')*CG =(I_q⊗Zint(s,s_dum)')(I_q⊗C1')*CG=
% (I_q⊗Zint(s)')Cnew(I⊗Zint(s_dum))  
% Cmul = (I_q⊗C1') and C = CG
Cout = Cmul*C;
[nrows, N] = size(Cout);
nz = nrows/q;
nA = (nz-1)/2;
Cnew = zeros(q*(nA+1), N*(nA+1));
for iq = 1:q
    rows_old = (iq-1)*nz + (1:nz);
    rows_new = (iq-1)*(nA+1) + (1:(nA+1));

    for j = 1:N
        cols_new = (j-1)*(nA+1) + (1:(nA+1));

        cvec = Cout(rows_old, j);

        c0 = cvec(1);
        cs = cvec(2:nA+1);
        cd = cvec(nA+2:end);

        % Matrix such that:
        % [1;ZA(s)]' * M * [1;ZA(s_dum)]
        % = c0 + cs'*ZA(s) + cd'*ZA(s_dum)
        M = [c0,  cd.';
             cs,  zeros(nA,nA)];

        Cnew(rows_new, cols_new) = M;
    end
end

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

Zint = cellfun(@(x) [0; x+1], Z, UniformOutput=false);
a = lims(:,1); b = lims(:,2);

for k = 1:length(Cidx)
idx = idxAll(k,:)';   % extract a specify beta_a or alpha_b
C = 1;                % this is C corresponding to that beta_a/alpha_b
for i=1:length(idx)
    E = Z{i};
    nE = length(E);
    Ci = zeros(nE,nE+1);
    switch idx(i)
        case 1
            % no integration along si
            Ci(:,1:end-1) = eye(nE);
        case 2
            % int_a^si si^E dsi = si^(E+1)/(E+1) - a^(E+1)/(E+1);
            Ci(:,1) = -a(i).^(E+1)./(E+1);  % constant of integration
            Ci(:,2:end) = diag(1./(E+1));   % si^(E+1)/(E+1) terms
        case 3
            % int_si^b si^E dsi = b^(E+1)/(E+1)-si^(E+1)/(E+1);
            Ci(:,1) = b(i).^(E+1)./(E+1);  % constant of integration
            Ci(:,2:end) = -diag(1./(E+1));   % -si^(E+1)/(E+1) terms
    end
    C = kron(C,Ci);
end
Cidx{i} = C;
end
end

<<<<<<< HEAD
function [Z2aout, G3aout, Z3bout] = int_2b(ZL, ZR, ZLvar, ZRvar,s2a,s2b,s3a,s3b,lims)
=======

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [Z2aout, G3aout, Z3bout] = int_2b(ZL, ZR, ZLvar, ZRvar, s2a,s2b,s3a,s3b,lims)
>>>>>>> 091796003159514d05a24a16bb99c5dccbac2651
% This performs the factorization
% int(ZL*ZR',s2b,0,1) = (Im\otimes Z2a') G(s3a) (In\otimes Z3b)
% where m = length(kron(ZL)) and n = length(kron(ZR))

% note that
% int_0^1 (ZL2a⊗ZL2b⊗ZL3a⊗ZL3b)(ZR2a⊗ZR2b⊗ZR3a⊗ZR3b)' d2b
% = int_0^1 (ZL2a*ZR2a')⊗(ZL2b*ZR2b')⊗(ZL3a*ZR3a')⊗(ZL3b*ZR3b') d2b
% = (ZL2a*ZR2a')⊗C2b⊗(ZL3a*ZR3a')⊗(ZL3b*ZR3b')
% for some constant matrix C2b

% first separate ZL and ZR into Z2a,Z2b,Z3a,Z3b
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

% condense Z2a matrix into (I2a \otimes Z2anew')*K2a
Z2anew = cell(1,length(s2a));
K2a = 1;
for i=1:length(s2a)
    E = Z2a{i};
    [u,~,idx] = unique(E(:));
    Z2anew{i} = u;
    Ki = sparse(idx,1:numel(E), 1, numel(u),numel(E));  % len(u)xlen(E) sparse matrix 
    K2a = kron(K2a,Ki);
end

% condense Z3b matrix into K3b*(I3b \otimes Z3bnew)
Z3bnew = cell(1,length(s3b));
K3b = 1;
for i=1:length(s3b)
    E = Z3b{i};
    [u,~,idx] = unique(E(:));
    Z3bnew{i} = u;
    Ki = sparse(1:numel(E),idx,1,numel(E),numel(u));
    K3b = kron(K3b,Ki);
end

% condense Z3a matrix into (I3a ⊗ Z3anew')*K3a
Z3anew = cell(1,length(s3a));
K3a = 1;
for i=1:length(s3a)
    E = Z3a{i};
    [u,~,idx] = unique(E(:));
    Z3anew{i} = u;
    Ki = sparse(idx,1:numel(E),1,length(u), numel(E));
    K3a = kron(K3a,Ki);
end

% now we have 
% (Z2a*Zp2a')⊗C2b⊗(Z3a*Zp3a')⊗(Z3b*Zp3b')
%  = ((I2a ⊗ Z2anew')*K2a)
%      ⊗C2b
%        ⊗(I3a ⊗ Z3anew')*K3a
%          ⊗(K3b*(I3b ⊗ Z3bnew))

% let us focus on finding
% K2a⊗C2b⊗((I3a ⊗ Z3anew')*K3a)⊗K3b
A = kron(K2a,C2b);
B = K3b;
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
n2a = prod(cellfun(@length,Z2ap));   % size of Z2anew block
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
Z2aout = Z2ap;
Z3bout = Z3bp;
G3aout = struct('C', C3a, 'Z', {Z3anew});
end



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function G = mapBetaAlpha2Gamma(alpha, beta)
% mapAlphaBetaToGamma  Map multi-indices alpha,beta (values in {1,2,3})
% to all possible gamma multi-indices under the specified per-component rule.
%
% Inputs:
%   alpha, beta : 1xd or d x 1 integer vectors with entries in {1,2,3}
%
% Output:
%   G : N x d integer matrix. Each row is one gamma multi-index.
%       N = product over i of number of possible gamma_i values.

    % Ensure row vectors
    if iscell(alpha)
        alpha = cell2mat(alpha);
    end
    if iscell(beta)
        beta = cell2mat(beta);
    end

    alpha = alpha(:).';
    beta  = beta(:).';
    d = numel(alpha);

    if numel(beta) ~= d
        error('alpha and beta must have the same length.');
    end

    % Build per-dimension option lists for gamma_i
    opts = cell(1, d);
    for i = 1:d
        opts{i} = gammaOptions(alpha(i), beta(i));
    end

    % Cartesian product across dimensions to form all gamma multi-indices
    grids = cell(1, d);
    [grids{:}] = ndgrid(opts{:});   % each grids{i} is an array of size |opts1|x...x|optsd|
    
    N = numel(grids{1});
    G = zeros(N, d);
    for i = 1:d
        G(:, i) = grids{i}(:);
    end
end

function gi = gammaOptions(a, b)
% gammaOptions  Return the vector of possible gamma values for one component.

    % Mapping table as cell array: row=a, col=b
    T = { ...
        1,   2,    3; ...
        2,   2,   [2 3]; ...
        3,  [2 3], 2 ...
    };

    gi = T{a, b};
end