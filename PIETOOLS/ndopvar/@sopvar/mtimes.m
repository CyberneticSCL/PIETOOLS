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

%B [s1,s3a,s3b] -> [s2a,s2b,s3a,s3b] S1=s1, S3 = [s3a,s3b], S2 = [s2a,s2b]
%A [s2a,s2b,s3a,s3b] -> [s4,s2a,s3a] S1=[s2b,s3b], S3 = [s2a,s3a], S2 = s4
vs1 = B.vars_S1; vs4 = A.vars_S2;
vs3a = intersect(A.vars_S3,B.vars_S3);
vs3b = setdiff(B.vars_S3,vs3a);
vs2a = setdiff(A.vars_S3,vs3a);
vs2b = setdiff(A.vars_S1,vs3b);

% Compute the vars, domain, and dimensions of the composition, C = A*B
Cdims = [A.dims(1),B.dims(2)];

n3a = numel(vs3a); % size of passthrough variable
if n3a == 0
    Csize = [1, 1];
elseif n3a == 1
    Csize = [3, 1];
else
    Csize = 3*ones(1,n3a);
end
% n3a=0 => Cparams (1x1), n3a=1 => Cparams (3x1), else Cparams (3x3x...x3) n3a times
Cparams = repmat({zeros(Cdims)}, Csize);
% Initialization of C parameters done. Monomials C.ZL/ZR will be done at the end.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

nalpha = numel(B.vars_S3); % pass through variables in B, s3a+s3b
nbeta = numel(A.vars_S3); % pass through variables in A, s2a+s3a

% convert linear index to multi-index for B parameters, e.g, if B is 3x3
% map 1,2...9 to [1,1],[2,1],...[3,3]
if nalpha == 0
    alphaIdx = zeros(1,0);
else
    alphaIdx = fliplr(dec2base(0:numel(B.params)-1,3,nalpha)-'0')+1;
end
if nbeta == 0
    betaIdx = zeros(1,0);
else
    betaIdx = fliplr(dec2base(0:numel(A.params)-1,3,nbeta)-'0')+1;
end

% split alpha and beta as [alpha_a, alpha_b] and [beta_a, beta_b]
% alpha_a corresponds to s3a vars, alpha_b to s3b
% beta_a corresponds to s2a vars, beta_b to s3a
% rebuild CA and CB indexed as CA{beta_a}{beta_b} and
% CB{alpha_a}{alpha_b}
% Desired sorted group orders:
%   beta_a  follows vs2a
%   beta_b  follows vs3a
%   alpha_a follows vs3a
%   alpha_b follows vs3b

[tf,colsbetaa] = ismember(vs2a,A.vars_S3);
if any(~tf)
    error('mtimes: vs2a not found in A.vars_S3.');
end

[tf,colsbetab] = ismember(vs3a,A.vars_S3);
if any(~tf)
    error('mtimes: vs3a not found in A.vars_S3.');
end

[tf,colsalphaa] = ismember(vs3a,B.vars_S3);
if any(~tf)
    error('mtimes: vs3a not found in B.vars_S3.');
end

[tf,colsalphab] = ismember(vs3b,B.vars_S3);
if any(~tf)
    error('mtimes: vs3b not found in B.vars_S3.');
end

beta_a  = unique(betaIdx(:,colsbetaa),'rows','stable');
beta_b  = unique(betaIdx(:,colsbetab),'rows','stable');

alpha_a = unique(alphaIdx(:,colsalphaa),'rows','stable');
alpha_b = unique(alphaIdx(:,colsalphab),'rows','stable');

if nbeta == 0
    Aparams = reshape(A.params,1,1);
    CA = Aparams;

elseif nbeta == 1
    Aparams = reshape(A.params,3,1);

    if numel(colsbetaa) == 1
        % single beta variable belongs to beta_a
        CA = reshape(Aparams,3,1);
    else
        % single beta variable belongs to beta_b
        CA = reshape(Aparams,1,3);
    end

else
    Aparams = reshape(A.params,3*ones(1,nbeta));

    permA = [colsbetaa,colsbetab];

    CA = permute(Aparams,permA);
    CA = reshape(CA,3^numel(colsbetaa),3^numel(colsbetab));
end
if nalpha == 0
    Bparams = reshape(B.params,1,1);
    CB = Bparams;

elseif nalpha == 1
    Bparams = reshape(B.params,3,1);

    if numel(colsalphaa) == 1
        % single alpha variable belongs to alpha_a
        CB = reshape(Bparams,3,1);
    else
        % single alpha variable belongs to alpha_b
        CB = reshape(Bparams,1,3);
    end

else
    Bparams = reshape(B.params,3*ones(1,nalpha));

    permB = [colsalphaa,colsalphab];

    % For true N-D case, permB has nalpha entries.
    CB = permute(Bparams,permB);
    CB = reshape(CB,3^numel(colsalphaa),3^numel(colsalphab));
end
if ~isequal(size(CA),[size(beta_a,1),size(beta_b,1)])
    error('mtimes: CA indexing failed. Expected CA(beta_a,beta_b).');
end

if ~isequal(size(CB),[size(alpha_a,1),size(alpha_b,1)])
    error('mtimes: CB indexing failed. Expected CB(alpha_a,alpha_b).');
end
% if numel(A.vars_S3)<=1
%     CA = (A.params)'; % nothing to reshape, no pass-through variables
% else
%     CA = reshape(permute(A.params, [colsbetaa,colsbetab]), ...
%         3^numel(colsbetaa), 3^numel(colsbetab)); % this is a cell 
% end
% if numel(B.vars_S3)<=1
%     CB = B.params; % nothing to reshape, no pass-through variables
% else
%     CB = reshape(permute(B.params, [colsalphaa,colsalphab]), ...
%         3^numel(colsalphaa), 3^numel(colsalphab));
% end

[tf,idx] = ismember(vs2b,B.vars_S2);
if any(~tf)
    error('mtimes: vs2b not found in B.vars_S2.');
end
dom2b = B.dom_2(idx,:);

[tf,idx] = ismember(vs2a,A.vars_S3);
if any(~tf)
    error('mtimes: vs2a not found in A.vars_S3.');
end
dom2a = A.dom_3(idx,:);

[tf,idx] = ismember(vs3a,B.vars_S3);
if any(~tf)
    error('mtimes: vs3a not found in B.vars_S3.');
end
dom3a = B.dom_3(idx,:);

[tf,idx] = ismember(vs3b,B.vars_S3);
if any(~tf)
    error('mtimes: vs3b not found in B.vars_S3.');
end
dom3b = B.dom_3(idx,:);

% len(A.ZR) = p, len(B.ZL) = q
p = prod(cellfun(@numel,A.ZR)); q = prod(cellfun(@numel,B.ZL));
% do int(A.ZR*B.ZL', s2b) =
% PL(I_p⊗Zhat_s2a')*G_s3a(s3a)*(I_q⊗Zhat_s3b)PR
% where G_s3a is size len(Zhat_s2a)*p x len(Zhat_s3b)*q = g1 x g2
[Zhat_s2a,G_s3a,Zhat_s3b, PL, PR] = ...
        int_2b(A.ZR,B.ZL,A.vars.in,vs2a,vs2b,vs3a,vs3b,dom2b);
Gs3adim = [p*prod(cellfun(@numel,Zhat_s2a)),q*prod(cellfun(@numel,Zhat_s3b))];

% alpha_b is attached to a B-pass-through variable which A fully integrates
% out.  Swapping the order of integration reverses lower and upper one-sided
% limits: int_z int_a^z = int_t^b and int_z int_z^b = int_a^t.
alpha_b_int = alpha_b;
alpha_b_int(alpha_b==2) = 3;
alpha_b_int(alpha_b==3) = 2;
[Cs2a_betaa, KZhat_s2a] = int_monomial(Zhat_s2a,beta_a,dom2a);   % note that KZhat_s2a has been expanded to ensure independence of beta_a
[Cs3b_alphab, KZhat_s3b] = int_monomial(Zhat_s3b,alpha_b_int,dom3b);  % same for KZhat_s3b
% Cs2a_betaa, Cs3b_alphab are cell with different coefficients for
% different beta_a/alpha_b stored in same order as beta_a and alpha_b variables



% now, we have
% int(A.ZR*B.ZL', s2a,s2b,s3a,s3b) 
% = (PL)(I_p⊗KZhat_s2a(s2a)'*Cs2a_betaa')*int(G_s3a(s3a),s3a)*(I_q⊗Cs3b_alphab*KZhat_s3b(s3b))(PR)
% = PL(I_p⊗KZhat_s2a(s2a)')*(I_p⊗Cs2a_betaa')*int(G_s3a(s3a),s3a)*(I_q⊗Cs3b_alphab)*(I_q⊗KZhat_s3b(s3b))PR
% next we need to perform int(G_s3a(s3a),s3a) = (I_g1⊗barZ_3aL')*C_(gam,beta,alp)*(I_g2⊗barZ_3aR)
[C_gam_alp_beta,barZ_3aL,barZ_3aR] = int_semisep(G_s3a,beta_b,alpha_a,dom3a, Csize); 




% looking at the original composition, we have  K^C_gamma = 
% sum_{betaa,betab,alphaa,alphab) (I⊗A.ZL')CA(betaa,betab)*int((I⊗A.ZR)*(I⊗B.ZL'),s2a,s2b,s3a,s3b)*CB(alphaa,alphab)(I⊗B.ZR)
% = sum (I⊗A.ZL')CA(betaa,betab)*(I⊗int(A.ZR*B.ZL',s2a,s2b,s3a,s3b))*CB(alphaa,alphab)(I⊗B.ZR)
% = sum (I⊗A.ZL')CA(betaa,betab)
%          *(I⊗(PL(I⊗KZhat_s2a(s2a)')*(I⊗Cs2a_betaa')*(I⊗barZ_3aL')*C_(gam,betab,alphaa)*(I⊗barZ_3aR)*(I⊗Cs3b_alphab)*(I⊗KZhat_s3b(s3b))PR))
%             *CB(alphaa,alphab)(I⊗B.ZR)

% split (I⊗ABCDE) as (I⊗A)(I⊗B)(I⊗C)....Then
% sum (I⊗A.ZL')CA(betaa,betab)
%       *(I⊗(PL(I⊗KZhat_s2a(s2a)')*(I⊗Cs2a_betaa')*(I⊗barZ_3aL')*C_(gam,betab,alphaa)*(I⊗barZ_3aR)*(I⊗Cs3b_alphab)*(I⊗KZhat_s3b(s3b))PR))
%          *CB(alphaa,alphab)(I⊗B.ZR)
% = sum (I⊗A.ZL')CA(betaa,betab)*(I⊗PL)*(I⊗KZhat_s2a(s2a)')*(I⊗Cs2a_betaa')
%         *(I⊗barZ_3aL')*(I⊗C_(gam,betab,alphaa))*(I⊗barZ_3aR)
%          *(I⊗Cs3b_alphab)*(I⊗KZhat_s3b(s3b))))*(I⊗PR)*CB(alphaa,alphab)(I⊗B.ZR)
% since none of the monomials depend on betas and alphas, there must be a
% template permutation operation that repeats for all indices to compute
% resulting composition in the quadPoly form (I⊗C.ZL') CC[gamma] (I⊗C.ZR) 




% sum_betaa (I⊗A.ZL')CA(betaa,betab)*((I⊗PL'))(I⊗KZhat_s2a(s2a)')*(I⊗Cs2a_betaa')
% = sum_betaa (I⊗ZLtemp(varLtemp)')CLtemp_betaab
% -------------------------------------------------------------------------
% Left side:
%
% sum_beta_a
%
%   (I \otimes A.ZL') CA(beta_a,beta_b)
%       * (I \otimes PL)
%       * (I \otimes KZhat_s2a')
%       * (I \otimes Cs2a_beta_a')
%
% is rewritten as
%
%   (I \otimes ZLtemp') CLtemp{beta_b}.
%
% Since int_monomial returns
%
%   I_beta Z = Cs2a_beta * KZhat_s2a,
%
% the coefficient used after KZhat_s2a' is Cs2a_beta'.
% -------------------------------------------------------------------------

nBetaA = size(beta_a,1);
nBetaB = size(beta_b,1);

CA = reshape(CA,nBetaA,nBetaB);

if numel(Cs2a_betaa) ~= nBetaA
    error('mtimes: Cs2a_betaa must have one cell for each beta_a.');
end

PLbig = kron(eye(size(CA{1},2)/size(PL,2)),PL);

CA_PL = cell(nBetaA,nBetaB);

for ia = 1:nBetaA
    for ib = 1:nBetaB
        CA_PL{ia,ib} = CA{ia,ib}*PLbig;
    end
end

qL = size(CA_PL{1},2);

CB_s2a = cell(nBetaA,1);

for ia = 1:nBetaA
    CB_s2a{ia} = kron(eye(qL),Cs2a_betaa{ia}.');
end

[varLtemp,ZLtemp,CLtemp_betaab] = leftShiftMonomials_SS( ...
    A.vars.out,A.ZL,CA_PL, ...
    vs2a,KZhat_s2a,CB_s2a);

CLtemp = cell(1,nBetaB);

for ib = 1:nBetaB
    CLtemp{ib} = CLtemp_betaab{1,ib};

    for ia = 2:nBetaA
        CLtemp{ib} = CLtemp{ib} + CLtemp_betaab{ia,ib};
    end
end
% -------------------------------------------------------------------------
% Right side:
%
% We need to absorb
%
%   (I \otimes Cs3b_alpha_b)
%   (I \otimes KZhat_s3b)
%   (I \otimes PR)
%   CB(alpha_a,alpha_b)
%   (I \otimes B.ZR)
%
% into the final right monomial basis.
%
% Reuse leftShiftMonomials_SS by transposing the matrix product:
%
%   [ ... CB(alpha_a,alpha_b) (I \otimes B.ZR) ]'
%
% This requires transposing both the cell layout and each matrix inside CB.
% -------------------------------------------------------------------------

nAlphaA = size(alpha_a,1);
nAlphaB = size(alpha_b,1);

% Make absolutely sure CB is indexed as CB(alpha_a,alpha_b)
CB = reshape(CB,nAlphaA,nAlphaB);

% Build CBt(alpha_b,alpha_a), with each matrix transposed.
CBt = cell(nAlphaB,nAlphaA);

for ia = 1:nAlphaA
    for ib = 1:nAlphaB
        CBt{ib,ia} = CB{ia,ib}.';
    end
end

% Apply PR' on the transposed side.
PRbig = kron(eye(size(CBt{1},2)/size(PR,1)),PR.');

CBt_PR = cell(size(CBt));

for i = 1:numel(CBt)
    CBt_PR{i} = CBt{i}*PRbig;
end

% Build one CB_s3b cell per alpha_b.
qR = size(CBt_PR{1},2);

CB_s3b = cell(nAlphaB,1);

for ib = 1:nAlphaB
    CB_s3b{ib} = kron(eye(qR),Cs3b_alphab{ib}.');
end

% Now CA rows = alpha_b, and CB_s3b has one cell per alpha_b.
[varRtemp,ZRtemp,CRtemp_alphaba] = leftShiftMonomials_SS( ...
    B.vars.in,B.ZR,CBt_PR, ...
    vs3b,KZhat_s3b,CB_s3b);

% CRtemp_alphaba is indexed as (alpha_b,alpha_a).
% Sum over alpha_b.
CRtemp = cell(nAlphaA,1);

for ia = 1:nAlphaA
    CRtemp{ia} = CRtemp_alphaba{1,ia};

    for ib = 2:nAlphaB
        CRtemp{ia} = CRtemp{ia} + CRtemp_alphaba{ib,ia};
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

% -------------------------------------------------------------------------
% Absorb barZ_3aL into the left side.
%
% Current left expression:
%
%   (I \otimes ZLtemp') CLtemp{beta_b}
%
% Need:
%
%   (I \otimes ZLtemp') CLtemp{beta_b}
%       * (I \otimes barZ_3aL')
%
% This is represented by multiplying with an identity coefficient after
% the barZ_3aL monomial factor.
% -------------------------------------------------------------------------

qLbar = size(CLtemp{1},2);
nBarL = prod(cellfun(@numel,barZ_3aL));
if isempty(barZ_3aL)
    nBarL = 1;
end

CB_barL = {speye(qLbar*nBarL)};

[Cvarsout,CZL,CLtemp_betab] = leftShiftMonomials_SS( ...
    varLtemp,ZLtemp,CLtemp, ...
    vs3a,barZ_3aL,CB_barL);


% -------------------------------------------------------------------------
% Absorb barZ_3aR into the right side.
%
% CRtemp{alpha_a} is still the coefficient for the transposed right factor.
%
% Current transposed right expression:
%
%   (I \otimes ZRtemp') CRtemp{alpha_a}
%
% Need:
%
%   (I \otimes ZRtemp') CRtemp{alpha_a}
%       * (I \otimes barZ_3aR')
%
% After this, we transpose the final coefficient in Cparams construction.
% -------------------------------------------------------------------------

CRtemp_row = CRtemp.';   % 1-by-nAlphaA

qRbar = size(CRtemp_row{1},2);
nBarR = prod(cellfun(@numel,barZ_3aR));
if isempty(barZ_3aR)
    nBarR = 1;
end

CB_barR = {speye(qRbar*nBarR)};

[Cvarsin,CZR,CRtemp_alphaa] = leftShiftMonomials_SS( ...
    varRtemp,ZRtemp,CRtemp_row, ...
    vs3a,barZ_3aR,CB_barR);


% -------------------------------------------------------------------------
% Assemble Cparams.
%
% CLtemp_betab{j}        : left coefficient for beta_b(j)
% C_gam_alp_beta{i,j,k} : middle coefficient for gamma(i), beta_b(j), alpha_a(k)
% CRtemp_alphaa{k}'     : right coefficient for alpha_a(k), transposed back
% -------------------------------------------------------------------------
nMid = A.dims(2);

for i = 1:numel(Cparams)

    Cparams{i} = sparse(size(CLtemp_betab{1},1),size(CRtemp_alphaa{1},1));

    for j = 1:size(beta_b,1)
        for k = 1:size(alpha_a,1)

            L = CLtemp_betab{j};
            M = C_gam_alp_beta{i,j,k};
            R = CRtemp_alphaa{k}.';

            if nMid ~= 1
                M = kron(speye(nMid),M);
            end

            if size(L,2) ~= size(M,1)
                error('mtimes: left/middle dimension mismatch in Cparams construction.');
            end

            if size(M,2) ~= size(R,1)
                error('mtimes: middle/right dimension mismatch in Cparams construction.');
            end

            Cparams{i} = Cparams{i} + L*M*R;
        end
    end
end


% -------------------------------------------------------------------------
% Construct C.
%
% Important: CZL is ordered according to Cvarsout, and CZR is ordered
% according to Cvarsin. Therefore Cvars/Cdom must use these same orders.
% -------------------------------------------------------------------------

[tf_in,idx_in] = ismember(Cvarsin,B.vars.in);
[tf_out,idx_out] = ismember(Cvarsout,A.vars.out);

if any(~tf_in)
    error('mtimes: final input variables do not match B input variables.');
end

if any(~tf_out)
    error('mtimes: final output variables do not match A output variables.');
end
% The semiseparable part above creates gamma indices only for variables that
% pass through the middle space (vs3a).  If B has an input-only variable with
% the same name as an A output-only variable, the final sopvar treats that
% name as a common variable.  The sequential operator is still a full integral
% in that variable, represented in 3-PI form by the sum of the lower and upper
% pieces and no multiplier piece.
Cvars_S3 = intersect(Cvarsin,Cvarsout);
extra_common = setdiff(Cvars_S3,vs3a);

if ~isempty(extra_common)
    [~,idx_extra_in] = ismember(extra_common,B.vars.in);
    [~,idx_extra_out] = ismember(extra_common,A.vars.out);

    if any(any(B.dom.in(idx_extra_in,:) ~= A.dom.out(idx_extra_out,:)))
        error('mtimes: matching final input/output variable names must have matching domains.');
    end
end

Cparams = expandCommonFullIntegrals(Cparams,vs3a,Cvars_S3);

Cvars = struct('in',{Cvarsin},'out',{Cvarsout});
Cdom  = struct('in',B.dom.in(idx_in,:), ...
               'out',A.dom.out(idx_out,:));

C = sopvar(Cparams,Cvars,CZR,CZL,Cdom,Cdims);
end

function Cparams_out = expandCommonFullIntegrals(Cparams_in,baseVars,finalVars)
% expandCommonFullIntegrals
%
% Cparams_in is indexed over baseVars.  finalVars may also contain variables
% that are common only because a B input-only variable has the same name as an
% A output-only variable.  For such variables the actual operation is a full
% integral, represented as gamma=2 plus gamma=3 and zero gamma=1.

baseVars = baseVars(:).';
finalVars = finalVars(:).';

if isequal(baseVars,finalVars)
    Cparams_out = Cparams_in;
    return
end

[tfBase,~] = ismember(baseVars,finalVars);
if any(~tfBase)
    error('expandCommonFullIntegrals: base variables must be final common variables.');
end

nBase = numel(baseVars);
nFinal = numel(finalVars);

finalIdx = paramMultiIdx_SS(nFinal);

if nFinal == 0
    Csize = [1,1];
elseif nFinal == 1
    Csize = [3,1];
else
    Csize = 3*ones(1,nFinal);
end

zeroParam = sparse(size(Cparams_in{1},1),size(Cparams_in{1},2));
Cparams_out = repmat({zeroParam},Csize);

[isBaseVar,basePos] = ismember(finalVars,baseVars);

for k = 1:size(finalIdx,1)
    gamma = finalIdx(k,:);

    % Extra common variables are full integrals, i.e. lower + upper only.
    if any((~isBaseVar) & (gamma == 1))
        Cparams_out{k} = zeroParam;
        continue
    end

    if nBase == 0
        Cparams_out{k} = Cparams_in{1};
    else
        gammaBase = zeros(1,nBase);
        gammaBase(basePos(isBaseVar)) = gamma(isBaseVar);
        linBase = paramMultiIdxToLinear_SS(gammaBase);
        Cparams_out{k} = Cparams_in{linBase};
    end
end
end

function idx = paramMultiIdx_SS(n)
% Parameter-cell multi-index order: first variable changes fastest.
if n == 0
    idx = zeros(1,0);
else
    idx = fliplr(dec2base(0:3^n-1,3,n)-'0') + 1;
end
end

function lin = paramMultiIdxToLinear_SS(idx)
idx = idx(:).';
lin = 1;
for i = 1:numel(idx)
    lin = lin + (idx(i)-1)*3^(i-1);
end
end

function [C_gam_alp_beta,ZL,ZR] = int_semisep(G,idxbeta,idxalpha,lims,Csize)
% int_semisep
%
% Converts
%
%   G(s3a) = (I_g1 \otimes ZG(s3a)')*CG
%
% into semiseparable coefficient matrices
%
%   I_{gamma,beta,alpha}G
%      = (I_g1 \otimes ZL(s3a)') ...
%            C_gam_alp_beta{gamma,beta,alpha} ...
%        (I_g2 \otimes ZR(s3a_dum))
%
% INPUTS
%   G.C      : coefficient matrix CG, size g1*NG by g2
%   G.Z      : 1-by-ns3a cell array of exponent vectors
%   idxbeta  : nbeta-by-ns3a array, beta_b indices
%   idxalpha : nalpha-by-ns3a array, alpha_a indices
%   lims     : ns3a-by-2 domain array
%   Csize    : size of final C.params over gamma, used for compatibility
%
% OUTPUTS
%   C_gam_alp_beta : cell array of size 3^ns3a by nbeta by nalpha
%   ZL             : left monomial basis in s3a
%   ZR             : right monomial basis in s3a_dum
%
% Each C_gam_alp_beta{k,i,j} has size
%
%   g1*prod(cellfun(@numel,ZL)) by g2*prod(cellfun(@numel,ZR))

ZG = G.Z(:).';
CG = G.C;

ns3a = numel(ZG);

% Normalize empty index arrays
if isempty(idxbeta)
    idxbeta = zeros(1,ns3a);
end
if isempty(idxalpha)
    idxalpha = zeros(1,ns3a);
end

nbeta  = size(idxbeta,1);
nalpha = size(idxalpha,1);

% Basic checks
if size(idxbeta,2) ~= ns3a
    error('int_semisep: idxbeta must have one column per s3a variable.');
end

if size(idxalpha,2) ~= ns3a
    error('int_semisep: idxalpha must have one column per s3a variable.');
end

if size(lims,1) ~= ns3a || size(lims,2) ~= 2
    error('int_semisep: lims must be numel(G.Z)-by-2.');
end

if any(idxbeta(:) < 1) || any(idxbeta(:) > 3) || ...
        any(idxbeta(:) ~= round(idxbeta(:)))
    error('int_semisep: idxbeta entries must be 1, 2, or 3.');
end

if any(idxalpha(:) < 1) || any(idxalpha(:) > 3) || ...
        any(idxalpha(:) ~= round(idxalpha(:)))
    error('int_semisep: idxalpha entries must be 1, 2, or 3.');
end

NG = prod(cellfun(@numel,ZG));

if NG == 0
    NG = 1;
end

if mod(size(CG,1),NG) ~= 0
    error('int_semisep: size(G.C,1) must be divisible by prod(numel(G.Z)).');
end

g1 = size(CG,1)/NG;
g2 = size(CG,2);

% No common pass-through variables.
% Then there is no semiseparable integration to perform.
if ns3a == 0
    C_gam_alp_beta = cell(1,nbeta,nalpha);

    for i = 1:nbeta
        for j = 1:nalpha
            C_gam_alp_beta{1,i,j} = CG;
        end
    end

    ZL = {};
    ZR = {};
    return
end

if nargin >= 5 && prod(Csize) ~= 3^ns3a
    error('int_semisep: Csize is inconsistent with the number of s3a variables.');
end

a = lims(:,1);
b = lims(:,2);

% Build common left/right bases.
% These must support:
%   G(s), G(t), int G from a to s/t, int G from s/t to b,
%   and int G between s and t.
ZL = cell(1,ns3a);

for i = 1:ns3a
    E = ZG{i}(:);

    if any(E == -1)
        error('int_semisep: exponent -1 cannot be integrated by polynomial monomial rules.');
    end

    ZL{i} = unique([0; E; E+1]);
end

ZR = ZL;

NL = prod(cellfun(@numel,ZL));
NR = NL;

% Gamma multi-indices in the same linear order as 3-by-3-by-... parameter cells.
gamIdx = fliplr(dec2base(0:3^ns3a-1,3,ns3a)-'0') + 1;

C_gam_alp_beta = cell(3^ns3a,nbeta,nalpha);

% Main loop over beta, alpha, gamma
for i = 1:nbeta

    beta = idxbeta(i,:);

    for j = 1:nalpha

        alpha = idxalpha(j,:);

        for k = 1:3^ns3a

            gamma = gamIdx(k,:);

            % Csep maps original ZG(eta) into the separated mixed basis
            % kron_i (ZL_i(s_i) \otimes ZR_i(t_i)).
            %
            % Csep has size
            %
            %   NG by NL*NR
            Csep = 1;

            for ell = 1:ns3a

                E = ZG{ell}(:);
                ZE = ZL{ell}(:);

                nE = numel(E);
                nZ = numel(ZE);

                % Local basis is kron(ZL_ell(s),ZR_ell(t)).
                % Column index is:
                %
                %   col = idx_t + (idx_s-1)*nZ
                %
                % so Ci has size nE by nZ^2.
                Ci = zeros(nE,nZ^2);

                key = 100*gamma(ell) + 10*beta(ell) + alpha(ell);

                for r = 1:nE

                    e = E(r);

                    switch key

                        % Evaluation at s
                        case {111,212,313}
                            coeffs = 1;
                            exp_s  = e;
                            exp_t  = 0;

                        % Evaluation at t
                        case {221,331}
                            coeffs = 1;
                            exp_s  = 0;
                            exp_t  = e;

                        % int_t^s eta^e deta
                        case 222
                            coeffs = [1/(e+1), -1/(e+1)];
                            exp_s  = [e+1, 0];
                            exp_t  = [0,   e+1];

                        % int_s^b eta^e deta
                        case 232
                            coeffs = [b(ell)^(e+1)/(e+1), -1/(e+1)];
                            exp_s  = [0, e+1];
                            exp_t  = [0, 0];

                        % int_t^b eta^e deta
                        case 332
                            coeffs = [b(ell)^(e+1)/(e+1), -1/(e+1)];
                            exp_s  = [0, 0];
                            exp_t  = [0, e+1];

                        % int_a^t eta^e deta
                        case 223
                            coeffs = [-a(ell)^(e+1)/(e+1), 1/(e+1)];
                            exp_s  = [0, 0];
                            exp_t  = [0, e+1];

                        % int_a^s eta^e deta
                        case 323
                            coeffs = [-a(ell)^(e+1)/(e+1), 1/(e+1)];
                            exp_s  = [0, e+1];
                            exp_t  = [0, 0];

                        % int_s^t eta^e deta
                        case 333
                            coeffs = [-1/(e+1), 1/(e+1)];
                            exp_s  = [e+1, 0];
                            exp_t  = [0,   e+1];

                        % Empty interval / incompatible ordering
                        otherwise
                            coeffs = [];
                            exp_s  = [];
                            exp_t  = [];
                    end

                    for h = 1:numel(coeffs)
                        is = find(ZE == exp_s(h),1);
                        it = find(ZE == exp_t(h),1);

                        if isempty(is) || isempty(it)
                            error('int_semisep: internal basis construction error.');
                        end

                        col = it + (is-1)*nZ;
                        Ci(r,col) = Ci(r,col) + coeffs(h);
                    end
                end

                Csep = kron(Csep,Ci);
            end

            if ~isequal(size(Csep),[NG,NL*NR])
                error('int_semisep: internal Csep size mismatch.');
            end

            % Convert
            %
            %   (I_g1 \otimes Zmix')*(I_g1 \otimes Csep')*CG
            %
            % into
            %
            %   (I_g1 \otimes ZL')*Cnew*(I_g2 \otimes ZR).
            Cmul = kron(speye(g1),Csep.');
            Cnew = rearrangeCoef(Cmul,CG,g1,g2,ZL);

            if ~isequal(size(Cnew),[g1*NL,g2*NR])
                error('int_semisep: output coefficient size mismatch.');
            end

            C_gam_alp_beta{k,i,j} = Cnew;
        end
    end
end

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function Cout = rearrangeCoef(Cmul,C,p,q,Zint)
% rearrangeCoef
%
% Converts
%
%   (I_p \otimes Zmix(s,t)') * Cmul*C
%
% where
%
%   Zmix(s,t) = kron_i (Z_i(s_i) \otimes Z_i(t_i))
%
% into
%
%   (I_p \otimes Z(s)') * Cout * (I_q \otimes Z(t)).
%
% INPUTS
%   Cmul : typically kron(I_p,Csep')
%   C    : size p*NG by q
%   p    : left matrix dimension
%   q    : right matrix dimension
%   Zint : 1-by-n cell basis for both s and t
%
% OUTPUT
%   Cout : size p*NZ by q*NZ
%
% where NZ = prod_i numel(Zint{i}).

Zint = Zint(:).';
ns = cellfun(@numel,Zint);

if isempty(ns)
    Cout = Cmul*C;
    return
end

NZ = prod(ns);

D = full(Cmul*C);

if ~isequal(size(D),[p*NZ^2,q])
    error('rearrangeCoef: input dimensions are inconsistent.');
end

n = numel(ns);

% D has rows ordered as:
%
%   [mixed monomial index, p-index]
%
% where the mixed monomial basis is
%
%   kron_i (Z_i(s_i) \otimes Z_i(t_i)).
%
% For each variable i, the local mixed basis is ordered as
%
%   kron(Z_i(s_i),Z_i(t_i)),
%
% so its local dimensions are [nt_i, ns_i] in MATLAB reshape order.
%
% First reshape into:
%
%   t_n, s_n, ..., t_1, s_1, p, q
%
% because MATLAB / kron ordering reverses dimensions under reshape.

D = reshape(D,[repelem(fliplr(ns),2), p, q]);

% Separate all s-indices from all t-indices.
%
% Current dimensions:
%   t_n, s_n, t_{n-1}, s_{n-1}, ..., t_1, s_1, p, q
%
% Desired dimensions:
%   s_n, ..., s_1, p, t_n, ..., t_1, q
%
% This gives row block (s,p) and column block (t,q).

perm = [2:2:2*n, 2*n+1, 1:2:2*n, 2*n+2];

D = permute(D,perm);

Cout = reshape(D,[NZ*p, NZ*q]);

end





%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [Cidx, Zint] = int_monomial(Z,idxAll,lims)
% int_monomial
%
% Performs one-sided PI monomial integration.
%
% For each row idx = idxAll(k,:), compute Cidx{k} such that
%
%       I_idx Z(s) = Cidx{k} * Zint(s)
%
% where Zint is a common expanded basis for all rows of idxAll.
%
% idx(i) = 1 : no integration
% idx(i) = 2 : int_{a_i}^{s_i} eta_i^Z deta_i
% idx(i) = 3 : int_{s_i}^{b_i} eta_i^Z deta_i
%
% INPUTS
%   Z      : 1xd cell array of exponent vectors
%   idxAll : nidx-by-d array with entries in {1,2,3}
%   lims   : d-by-2 domain array
%
% OUTPUTS
%   Cidx   : 1-by-nidx cell array
%            each Cidx{k} has size N-by-M
%   Zint   : 1xd cell array common expanded basis
%
% where
%   N = prod_i numel(Z{i})
%   M = prod_i numel(Zint{i})

Z = Z(:).';
d = numel(Z);

% Empty tensor-product basis
if d == 0
    if isempty(idxAll)
        idxAll = zeros(1,0);
    end

    Cidx = cell(1,size(idxAll,1));
    for k = 1:numel(Cidx)
        Cidx{k} = 1;
    end

    Zint = {};
    return
end

% Basic checks
if size(idxAll,2) ~= d
    error('int_monomial: idxAll must have one column per variable in Z.');
end

if size(lims,1) ~= d || size(lims,2) ~= 2
    error('int_monomial: lims must be numel(Z)-by-2.');
end

if any(idxAll(:) < 1) || any(idxAll(:) > 3) || any(idxAll(:) ~= round(idxAll(:)))
    error('int_monomial: idxAll entries must be 1, 2, or 3.');
end

a = lims(:,1);
b = lims(:,2);

% Build common expanded output basis.
% This is the important correction: Zint must work for every row of idxAll.
Zint = cell(1,d);

for i = 1:d
    E = Z{i}(:);

    if any((idxAll(:,i)==2) | (idxAll(:,i)==3)) && any(E == -1)
        error('int_monomial: exponent -1 cannot be integrated by polynomial monomial rules.');
    end

    Ei = [];

    % Needed for idx = 1
    if any(idxAll(:,i)==1)
        Ei = [Ei; E];
    end

    % Needed for idx = 2 or 3
    if any((idxAll(:,i)==2) | (idxAll(:,i)==3))
        Ei = [Ei; 0; E+1];
    end

    Zint{i} = unique(Ei);
end

N = prod(cellfun(@numel,Z));
M = prod(cellfun(@numel,Zint));

Cidx = cell(1,size(idxAll,1));

% Build coefficient matrix for each row of idxAll
for k = 1:size(idxAll,1)

    idx = idxAll(k,:);
    C = 1;

    for i = 1:d
        E  = Z{i}(:);
        Ei = Zint{i}(:);

        nE = numel(E);
        nI = numel(Ei);

        Ci = zeros(nE,nI);

        switch idx(i)

            case 1
                % No integration: s^E
                for j = 1:nE
                    col = find(Ei == E(j),1);
                    Ci(j,col) = 1;
                end

            case 2
                % int_a^s eta^E deta
                % = s^(E+1)/(E+1) - a^(E+1)/(E+1)
                col0 = find(Ei == 0,1);

                for j = 1:nE
                    colp = find(Ei == E(j)+1,1);

                    Ci(j,col0) = Ci(j,col0) - a(i).^(E(j)+1)./(E(j)+1);
                    Ci(j,colp) = Ci(j,colp) + 1./(E(j)+1);
                end

            case 3
                % int_s^b eta^E deta
                % = b^(E+1)/(E+1) - s^(E+1)/(E+1)
                col0 = find(Ei == 0,1);

                for j = 1:nE
                    colp = find(Ei == E(j)+1,1);

                    Ci(j,col0) = Ci(j,col0) + b(i).^(E(j)+1)./(E(j)+1);
                    Ci(j,colp) = Ci(j,colp) - 1./(E(j)+1);
                end
        end

        C = kron(C,Ci);
    end

    if ~isequal(size(C),[N,M])
        error('int_monomial: internal coefficient size mismatch.');
    end

    Cidx{k} = C;
end

end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [Z2aout, G3aout, Z3bout, PL, PR] = ...
            int_2b(ZL, ZR, Zvar, s2a, s2b, s3a, s3b, lims)
% int_2b
% This performs the factorization
%
%   int_{s2b} ZL*ZR' ds2b
%     = PL*(Im \otimes Z2aout')*G3aout(s3a)*(In \otimes Z3bout)*PR
%
% where
%
%   G3aout(s3a) = (I_{m*nz2a} \otimes Z3aout')*G3aout.C
%
% INPUTS
%   ZL   : cell array of exponent vectors for A.ZR, ordered by Zvar
%   ZR   : cell array of exponent vectors for B.ZL, ordered by Zvar
%   Zvar : variables shared by A input and B output
%   s2a  : variables in A.vars_S3 but not common pass-through of C
%   s2b  : variables integrated out completely
%   s3a  : variables that pass through C
%   s3b  : variables in B.vars_S3 but not common pass-through of C
%   lims : numel(s2b)-by-2 domain array for s2b
%
% OUTPUTS
%   Z2aout : condensed monomial basis for products in s2a
%   G3aout : struct with fields
%              G3aout.C : coefficient matrix
%              G3aout.Z : condensed monomial basis in s3a
%   Z3bout : condensed monomial basis for products in s3b
%   PL     : permutation matrix for left monomial ordering
%   PR     : permutation matrix for right monomial ordering
%
% Dimensions:
%   m       = prod_i length(ZL{i})
%   n       = prod_i length(ZR{i})
%   nz2a    = prod_i length(Z2aout{i})
%   nz3a    = prod_i length(G3aout.Z{i})
%   nz3b    = prod_i length(Z3bout{i})
%
%   PL                  : m x m
%   PR                  : n x n
%   G3aout.C            : (m*nz2a*nz3a) x (n*nz3b)
%   G3aout(s3a)         : (m*nz2a) x (n*nz3b)
%   integrated bridge   : m x n

% Normalize cell-array orientation
ZL   = ZL(:).';
ZR   = ZR(:).';
Zvar = Zvar(:).';

s2a = s2a(:).';
s2b = s2b(:).';
s3a = s3a(:).';
s3b = s3b(:).';

Zvar_ordered = [s2a, s2b, s3a, s3b];

% Basic consistency checks
if numel(ZL) ~= numel(Zvar) || numel(ZR) ~= numel(Zvar)
    error('int_2b: ZL, ZR, and Zvar must have the same number of variables.');
end

if numel(Zvar_ordered) ~= numel(Zvar) || ...
        any(~ismember(Zvar_ordered, Zvar)) || any(~ismember(Zvar, Zvar_ordered))
    error('int_2b: [s2a,s2b,s3a,s3b] must be a permutation of Zvar.');
end

if size(lims,1) ~= numel(s2b) || size(lims,2) ~= 2
    error('int_2b: lims must be numel(s2b)-by-2.');
end

% Locate variable groups in the original Zvar ordering
[~, loc2a] = ismember(s2a, Zvar);
[~, loc2b] = ismember(s2b, Zvar);
[~, loc3a] = ismember(s3a, Zvar);
[~, loc3b] = ismember(s3b, Zvar);

ZL2a = ZL(loc2a);   ZR2a = ZR(loc2a);
ZL2b = ZL(loc2b);   ZR2b = ZR(loc2b);
ZL3a = ZL(loc3a);   ZR3a = ZR(loc3a);
ZL3b = ZL(loc3b);   ZR3b = ZR(loc3b);

% Monomial dimensions in ordered variable grouping
dimL2a = cellfun(@numel, ZL2a);
dimL2b = cellfun(@numel, ZL2b);
dimL3a = cellfun(@numel, ZL3a);
dimL3b = cellfun(@numel, ZL3b);

dimR2a = cellfun(@numel, ZR2a);
dimR2b = cellfun(@numel, ZR2b);
dimR3a = cellfun(@numel, ZR3a);
dimR3b = cellfun(@numel, ZR3b);

dimsLord = [dimL2a, dimL2b, dimL3a, dimL3b];
dimsRord = [dimR2a, dimR2b, dimR3a, dimR3b];

m = prod(dimsLord);
n = prod(dimsRord);

% -------------------------------------------------------------------------
% Build PL and PR.
%
% The middle factor is built in ordered variable order:
%   [s2a,s2b,s3a,s3b].
%
% PL and PR convert it back to the original Zvar ordering:
%
%   ZL_original = PL*ZL_ordered
%   ZR_original = SR*ZR_ordered
%
% Hence:
%
%   ZL_original*ZR_original'
%     = PL*(ZL_ordered*ZR_ordered')*SR'
%
% and we return:
%
%   PR = SR'
% -------------------------------------------------------------------------

dimLold = cellfun(@numel, ZL);
dimRold = cellfun(@numel, ZR);

% [~, pL] = ismember(Zvar_ordered, Zvar);   % ordered vars = old vars(pL)
% [~, pR] = ismember(Zvar_ordered, Zvar);

% if numel(dimLold) <= 1
%     PLvec = (1:m).';
% else
%     invpL = zeros(size(pL));
%     invpL(pL) = 1:numel(pL);
%     PLvec = reshape(permute(reshape(1:m, dimLold(pL)), invpL), [], 1);
% end
% 
% if numel(dimRold) <= 1
%     PRvec = (1:n).';
% else
%     invpR = zeros(size(pR));
%     invpR(pR) = 1:numel(pR);
%     PRvec = reshape(permute(reshape(1:n, dimRold(pR)), invpR), [], 1);
% end

% Monomial vectors in PIETOOLS are ordered with the first variable outermost
% and the last variable fastest.  A reshape/permute on dimLold uses MATLAB's
% first-dimension-fastest convention and gives the wrong permutation for
% nontrivial variable reordering.  Build the permutation from explicit
% tensor-product multi-indices instead.
PL = kronPermutationMatrix_SS(Zvar,Zvar_ordered,dimLold);
SR = kronPermutationMatrix_SS(Zvar,Zvar_ordered,dimRold);
PR = SR.';

% -------------------------------------------------------------------------
% Condense product monomial bases in s2a, s3a, and s3b.
%
% For each variable si:
%   ZL_i(si)*ZR_i(si)' has exponent matrix ZL_i + ZR_i'
% We store only unique exponents, and keep index maps back into that basis.
% -------------------------------------------------------------------------

Z2anew = cell(1, numel(s2a));
idx2a  = cell(1, numel(s2a));
for i = 1:numel(s2a)
    E = ZL2a{i}(:) + ZR2a{i}(:).';
    [Z2anew{i}, ~, idx] = unique(E(:));
    idx2a{i} = reshape(idx, size(E));
end

Z3anew = cell(1, numel(s3a));
idx3a  = cell(1, numel(s3a));
for i = 1:numel(s3a)
    E = ZL3a{i}(:) + ZR3a{i}(:).';
    [Z3anew{i}, ~, idx] = unique(E(:));
    idx3a{i} = reshape(idx, size(E));
end

Z3bnew = cell(1, numel(s3b));
idx3b  = cell(1, numel(s3b));
for i = 1:numel(s3b)
    E = ZL3b{i}(:) + ZR3b{i}(:).';
    [Z3bnew{i}, ~, idx] = unique(E(:));
    idx3b{i} = reshape(idx, size(E));
end

% Integrate the s2b monomial products over their full domains
C2b = cell(1, numel(s2b));
a = lims(:,1);
b = lims(:,2);

for i = 1:numel(s2b)
    E = ZL2b{i}(:) + ZR2b{i}(:).';
    C2b{i} = (b(i).^(E+1) - a(i).^(E+1))./(E+1);
end

nz2a = prod(cellfun(@numel, Z2anew));
nz3a = prod(cellfun(@numel, Z3anew));
nz3b = prod(cellfun(@numel, Z3bnew));

% -------------------------------------------------------------------------
% Build multi-index tables for ordered Kronecker bases.
%
% IL(i,k) is the local monomial index of the i-th ordered ZL monomial
% in the k-th ordered variable.
%
% IR(j,k) is the analogous table for ZR.
% -------------------------------------------------------------------------

if isempty(dimsLord)
    IL = zeros(1,0);
elseif numel(dimsLord) == 1
    IL = (1:m).';
else
    tmp = cell(1, numel(dimsLord));
    [tmp{:}] = ind2sub(fliplr(dimsLord), 1:m);
    IL = flipud(vertcat(tmp{:})).';
end

if isempty(dimsRord)
    IR = zeros(1,0);
elseif numel(dimsRord) == 1
    IR = (1:n).';
else
    tmp = cell(1, numel(dimsRord));
    [tmp{:}] = ind2sub(fliplr(dimsRord), 1:n);
    IR = flipud(vertcat(tmp{:})).';
end

% Ordered group offsets
n2a = numel(s2a);
n2b = numel(s2b);
n3a = numel(s3a);
n3b = numel(s3b);

L2a = 1:n2a;
L2b = n2a + (1:n2b);
L3a = n2a+n2b + (1:n3a);
L3b = n2a+n2b+n3a + (1:n3b);

R2a = L2a;
R2b = L2b;
R3a = L3a;
R3b = L3b;

% -------------------------------------------------------------------------
% Assemble C3a.
%
% For each ordered pair of monomials i,j:
%
%   ZL_i * ZR_j integrated over s2b
%
% contributes to exactly one basis index in s2a, one in s3a, and one in
% s3b. The coefficient from s2b integration is scalar c.
%
% C3a is arranged so that:
%
%   G3aout(s3a)
%     = (I_{m*nz2a} \otimes Z3anew')*C3a
%
% and therefore:
%
%   (I_m \otimes Z2anew')*G3aout(s3a)*(I_n \otimes Z3bnew)
%
% has size m-by-n.
% -------------------------------------------------------------------------

II = zeros(m*n, 1);
JJ = zeros(m*n, 1);
VV = zeros(m*n, 1);

cnt = 0;

dims2a = cellfun(@numel, Z2anew);
dims3a = cellfun(@numel, Z3anew);
dims3b = cellfun(@numel, Z3bnew);

for i = 1:m
    for j = 1:n

        % Scalar contribution from full integration over s2b
        c = 1;
        for k = 1:n2b
            c = c * C2b{k}(IL(i,L2b(k)), IR(j,R2b(k)));
        end

        % Condensed basis index for s2a products
        if n2a == 0
            k2a = 1;
        else
            idx = zeros(1,n2a);
            for k = 1:n2a
                idx(k) = idx2a{k}(IL(i,L2a(k)), IR(j,R2a(k)));
            end

            if n2a == 1
                k2a = idx;
            else
                idxcell = num2cell(fliplr(idx));
                k2a = sub2ind(fliplr(dims2a), idxcell{:});
            end
        end

        % Condensed basis index for s3a products
        if n3a == 0
            k3a = 1;
        else
            idx = zeros(1,n3a);
            for k = 1:n3a
                idx(k) = idx3a{k}(IL(i,L3a(k)), IR(j,R3a(k)));
            end

            if n3a == 1
                k3a = idx;
            else
                idxcell = num2cell(fliplr(idx));
                k3a = sub2ind(fliplr(dims3a), idxcell{:});
            end
        end

        % Condensed basis index for s3b products
        if n3b == 0
            k3b = 1;
        else
            idx = zeros(1,n3b);
            for k = 1:n3b
                idx(k) = idx3b{k}(IL(i,L3b(k)), IR(j,R3b(k)));
            end

            if n3b == 1
                k3b = idx;
            else
                idxcell = num2cell(fliplr(idx));
                k3b = sub2ind(fliplr(dims3b), idxcell{:});
            end
        end

        % Row of C3a:
        %   k3a inside block (i,k2a)
        rowBlock = k2a + (i-1)*nz2a;
        row = k3a + (rowBlock-1)*nz3a;

        % Column of C3a:
        %   k3b inside block j
        col = k3b + (j-1)*nz3b;

        cnt = cnt + 1;
        II(cnt) = row;
        JJ(cnt) = col;
        VV(cnt) = c;
    end
end

C3a = sparse(II, JJ, VV, m*nz2a*nz3a, n*nz3b);

% Outputs
Z2aout = Z2anew;
Z3bout = Z3bnew;

G3aout = struct();
G3aout.C = C3a;
G3aout.Z = Z3anew;

end
function P = kronPermutationMatrix_SS(oldVars,newVars,dimsOld)
% kronPermutationMatrix_SS returns P such that Z_old = P*Z_new for tensor
% monomial vectors ordered with the first variable outermost and the last
% variable fastest.

oldVars = oldVars(:).';
newVars = newVars(:).';
dimsOld = dimsOld(:).';

if numel(oldVars) ~= numel(newVars) || numel(dimsOld) ~= numel(oldVars)
    error('kronPermutationMatrix_SS: inconsistent inputs.');
end

[tf,loc] = ismember(newVars,oldVars);
if any(~tf) || numel(unique(loc)) ~= numel(loc)
    error('kronPermutationMatrix_SS: newVars must be a permutation of oldVars.');
end

N = prod(dimsOld);
if isempty(dimsOld)
    N = 1;
end

if numel(dimsOld) <= 1
    P = speye(N);
    return
end

dimsNew = dimsOld(loc);
idxNew = kronMultiIdx_SS(dimsNew);
idxOld = zeros(N,numel(dimsOld));
idxOld(:,loc) = idxNew;
oldLin = kronSub2ind_SS(dimsOld,idxOld);

P = sparse(oldLin,(1:N).',1,N,N);
end

function idx = kronMultiIdx_SS(dims)
% Multi-indices for the PIETOOLS tensor order: first variable outermost,
% last variable fastest.

dims = dims(:).';
N = prod(dims);
if isempty(dims)
    idx = zeros(1,0);
elseif numel(dims) == 1
    idx = (1:N).';
else
    tmp = cell(1,numel(dims));
    [tmp{:}] = ind2sub(fliplr(dims),1:N);
    idx = flipud(vertcat(tmp{:})).';
end
end

function lin = kronSub2ind_SS(dims,idx)
% Inverse of kronMultiIdx_SS.

dims = dims(:).';
if isempty(dims)
    lin = ones(size(idx,1),1);
elseif numel(dims) == 1
    lin = idx(:,1);
else
    subs = num2cell(fliplr(idx),1);
    lin = sub2ind(fliplr(dims),subs{:});
    lin = lin(:);
end
end