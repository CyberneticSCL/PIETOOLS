function C = mtimes(A,B)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% C = MTIMES(A,B) returns the composition C = A*B of two 'sopvar' objects,
% the output space of B being the input space of A, or the product of a
% 'sopvar' object with a scalar or with a matrix of matching dimension.
%
% With B: L_2^p[S1,S3] -> L_2^q[S2,S3] and A: L_2^q[S4,S5] -> L_2^r[S6,S5],
% where [S4,S5] and [S2,S3] are the same variables, the composition is
%   C = A*B: L_2^p[S7,S8] -> L_2^r[S9,S8],
%   S7 = S1 and S3 cap S4 (input only), S8 = S3 cap S5 (shared),
%   S9 = S6 and S2 cap S5 (output only);
% the variables S2 cap S4, created by B and integrated out by A, do not
% appear in C. The kernels of C are the semiseparable integrals of the
% products of the kernels of A and B (sopvar document Sec. 6).
%
% INPUT
% - A, B:   'sopvar' objects with A.vars.in the same set of variables as
%           B.vars.out and A.dims(2) equal to B.dims(1); or one of them a
%           'double', a scalar or a matrix (m x A.dims(1), or B.dims(2) x n)
%           acting on the output or input components; or a 'dpvar' object,
%           in which case the call is routed to 'dpvar_op_copvar';
%
% OUTPUT
% - C:      'sopvar' object representing A*B, of dims [A.dims(1),B.dims(2)],
%           on the monomial bases the composition reaches, pruned of the
%           degrees that carry no coefficient; an 'sdopvar' object when a
%           factor is a 'dpvar'.
%
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
% MMP, 09/06/2026: A numeric matrix factor is lifted over the monomial index
%                  and the new dimension recorded; before, a matrix applied
%                  straight to the coefficients failed for degree >= 1 and
%                  left 'dims' wrong for degree 0.
% MMP, 09/07/2026: Variable lists compared with 'isequal', not
%                  'any(~strcmp(...))', which is false for lists of
%                  different length.
% MMP, 09/09/2026: B's output side is aligned onto A.vars.in by NAME: the
%                  canonical orders of the two lists differ (for A'*A the
%                  adjoint's vars.in is [S3,S2], A's vars.out [S2,S3]), so
%                  the variable SETS are compared and B.ZL, B.dom.out and
%                  the row index of each parameter are permuted onto A's
%                  order before the positional contraction below.
% MMP, 09/10/2026: int_2b assembles C3a vectorized over the (i,j) grid
%                  ('condenseIdx'); the scalar loop called num2cell, sub2ind
%                  and fliplr per pair. Exact: the triplet positions are
%                  bijections, so 'sparse' sums no duplicates.
% MMP, 09/10/2026: PART 4 forms only the nMid diagonal blocks of
%                  LeftMat*CMat*RightMat (sopvar document Sec. 5.4), not the
%                  full product and a dense mask.
% MMP, 09/11/2026: The composed bases CZL/CZR are pruned of the degrees that
%                  hold no coefficient before C is built
%                  ('prune_zero_monomials'): PART 3 lists the degrees the
%                  integrals can reach, most of which are empty and were
%                  redundant rows in every later equality constraint.
% MMP, 09/29/2026: A dpvar factor is routed to 'dpvar_op_copvar' (a scalar
%                  gives the sdopvar (a + b'd)*P, a matrix acts on the
%                  components); gam*P failed in the dims check before.
% DJ, 11/01/2026: Bugfix for scalar multiplication with zero blocks
%                   specified as [] or scalar 0.
% MMP, 11/06/2026: Comment pass (maintainer, 11/06/2026): help block in the
%                  house format, log entries condensed to defect and fix,
%                  commented-out code of earlier versions and debug lines
%                  removed, notes tightened, INPUTS/OUTPUTS for the local
%                  functions; code unchanged.
%

% dpvar factor: legacy semantics, before any work; the composition below    % MMP, 09/29/2026
% is untouched.                                                             % MMP, 09/29/2026
if isa(A,'dpvar') || isa(B,'dpvar')                                         % MMP, 09/29/2026
    C = dpvar_op_copvar('mtimes',A,B);                                      % MMP, 09/29/2026
    return                                                                  % MMP, 09/29/2026
end                                                                         % MMP, 09/29/2026

% Separately deal with scalar multiplication
if isnumeric(A)
    if ~isscalar(A) && size(A,2)~=size(B,1)
        error("Inner dimensions of matrix-valued objects must match for multiplication.")
    end
    C = B;
    if ~isscalar(A)                                                         % MMP, 09/06/2026
        % A matrix factor acts on the row index lifted over ZL:             % MMP, 09/06/2026
        %   M*(I_m kron ZL') = (I_p kron ZL')*kron(M,I_NL).                 % MMP, 09/06/2026
        C.dims(1) = size(A,1);                                              % MMP, 09/06/2026
        A = kron(sparse(A),speye(prod([cellfun(@numel,B.ZL),1])));          % MMP, 09/06/2026
    end                                                                     % MMP, 09/06/2026
    % Zero-blocks declared as [] or scalar 0 stay as is.
    Bparams = B.params;
    is_zeroB = cellfun(@(b) isempty(b) || (isscalar(b) && b==0), Bparams);  % DJ, 11/01/2026
    Bparams(~is_zeroB) = cellfun(@(b) A*b, Bparams(~is_zeroB),'UniformOutput',false);
    C.params = Bparams;
    return
elseif isnumeric(B)
    if ~isscalar(B) && size(A,2)~=size(B,1)
        error("Inner dimensions of matrix-valued objects must match for multiplication.")
    end
    C = A;
    if ~isscalar(B)                                                         % MMP, 09/06/2026
        C.dims(2) = size(B,2);                                              % MMP, 09/06/2026
        B = kron(sparse(B),speye(prod([cellfun(@numel,A.ZR),1])));          % MMP, 09/06/2026
    end                                                                     % MMP, 09/06/2026
    % Zero-blocks declared as [] or scalar 0 stay as is.
    Aparams = A.params;
    is_zeroA = cellfun(@(a) isempty(a) || (isscalar(a) && a==0), Aparams);  % DJ, 11/01/2026
    Aparams(~is_zeroA) = cellfun(@(a) a*B, Aparams(~is_zeroA),'UniformOutput',false); 
    C.params = Aparams;
    return
end

% Error handling: Checks to ensure A and B are compatible
if A.dims(2)~=B.dims(1)
    error('number of output components of B is different than number of input components of A');
end
% Legality is a question about the variable SETS. The composition body      % MMP, 09/09/2026
% below contracts A.ZR against B.ZL positionally, so the two lists must     % MMP, 09/09/2026
% then agree in ORDER -- and they need not, because each operand is         % MMP, 09/09/2026
% stored in its own canonical order: for A'*A the adjoint's vars.in is      % MMP, 09/09/2026
% [S3,S2] while A's vars.out is [S2,S3]. So B is aligned onto A.vars.in     % MMP, 09/09/2026
% by name below rather than the composition being refused.                  % MMP, 09/09/2026
if ~isempty(setxor(A.vars.in,B.vars.out))                                   % MMP, 09/09/2026
%if ~isequal(A.vars.in(:),B.vars.out(:))                                    % MMP, 09/09/2026 (was)
%if any(~strcmp(A.vars.in,B.vars.out))                                      % MMP, 09/07/2026 (was)
    error('names of B output variables differ from names of A input variables');
end

% Permute B's output side onto A's input order: reorder B.ZL, B.dom.out     % MMP, 09/09/2026
% and the ROW index of every parameter, the row being (matrix row outer,    % MMP, 09/09/2026
% ZL monomial inner). Done in place on this local copy, NOT through the     % MMP, 09/09/2026
% constructor, which would canonicalize B straight back. The permutation    % MMP, 09/09/2026
% touches only rows, so the canonical multiplier form -- a condition on     % MMP, 09/09/2026
% ZR columns -- is untouched, and vars_S1/S2/S3 are sorted and so are       % MMP, 09/09/2026
% unchanged by a reordering.                                                % MMP, 09/09/2026
if ~isequal(A.vars.in(:),B.vars.out(:))                                     % MMP, 09/09/2026
    [~,ordB] = ismember(A.vars.in,B.vars.out);                              % MMP, 09/09/2026
    pB   = monomial_gather(B.ZL,ordB);                                      % MMP, 09/09/2026
    NLb  = prod([cellfun(@numel,B.ZL),1]);                                  % MMP, 09/09/2026
    rowB = reshape(pB(:)+(0:B.dims(1)-1)*NLb,[],1);                         % MMP, 09/09/2026
    for ii = 1:numel(B.params)                                              % MMP, 09/09/2026
        B.params{ii} = B.params{ii}(rowB,:);                                % MMP, 09/09/2026
    end                                                                     % MMP, 09/09/2026
    B.ZL      = B.ZL(ordB);                                                 % MMP, 09/09/2026
    B.dom.out = B.dom.out(ordB,:);                                          % MMP, 09/09/2026
    B.vars.out = B.vars.out(ordB);                                          % MMP, 09/09/2026
end                                                                         % MMP, 09/09/2026

if any(any(A.dom.in~=B.dom.out))
    error('output domain of B differs from input domain of A');
end

% Partition of the middle variables (the output of B, input of A):
%   vs3a  pass through both operators: the shared variables of C;
%   vs3b  pass through B, integrated out by A: input-only variables of C;
%   vs2a  created by B, pass through A: output-only variables of C;
%   vs2b  created by B, integrated out by A: absent from C.
% The input-only variables of B and the output-only variables of A stay.
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
Cparams = repmat({sparse(Cdims(1),Cdims(2))}, Csize);
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

    if ~isempty(colsbetaa)
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

    if ~isempty(colsalphaa)
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

% Domains of the four variable groups
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

%%%%%%%%%%%%%%%%%%%%%%%%%% PART 1 %%%%%%%%%%%%%%%%%%%%%%%%%%
[CMat,barZ_3aL,barZ_3aR] = int_semisep(G_s3a,beta_b,alpha_a,dom3a,Csize,'packed'); 

%%%%%%%%%%%%%%%%%%%%%%%%%% PART 2 %%%%%%%%%%%%%%%%%%%%%%%%%%
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


CRtemp_alphaatransposed = cell(size(CRtemp_alphaa));
for k = 1:numel(CRtemp_alphaa)
    CRtemp_alphaatransposed{k} = CRtemp_alphaa{k}.';
end





%%%%%%%%%%%%%%%%%%%%%%%%%% PART 3 %%%%%%%%%%%%%%%%%%%%%%%%%%
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
% We need to compute 
% sum_{j, k} L_j [I_{nMid} otimes C_gam_alp_beta{i, j, k}] R_k
% 
% If we have L_j = [L_{j1}, ... L_{j,nMid}] and 
%            R_k = [R_{k1}; ...; R_{k, nMid}], then
%
% sum_{j, k} L_j [I_{nMid} otimes C_gam_alp_beta{i, j, k}] R_k=
% sum_{j, k, d} L_{jd} C_gam_alp_beta{i, j, k} R_{kd}
%
% To compute the expression, we construct block matrices
%
% [LeftMat]_{d, j} = L_{j, d} and
% [RightMat]_{k, d} = R_{k, d}
% [CMat]_{j, k} = C_gam_alp_beta{i, j, k} (computed in int_semisep_AT)
% Then, 
%
% sum_{j, k} L_j [I_{nMid} otimes C_gam_alp_beta{i, j, k}] R_k
% is a block diagonal of LeftMat*CMat*RightMat,
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
% construct LeftMat, RightMat, CMat


nMid = A.dims(2);
Left_CLtemp_betab = cell2mat(CLtemp_betab); %[CL{1}, ... ,CL{end}]
Right_CRtemp_alphaaT = cell2mat(CRtemp_alphaatransposed'); %[CR{1}; ... ;CR{end}]

% % % % STEP 1 % % % % 
% Left_CLtemp=> LeftMat
% Left_CLtemp has the form
% [L_{11}, .., L_{1d}, ... , L_{n1}, .. L_{nd}] = [L_{jk}];
% we need [L_{11}, ..., L_{n1};
%          ... 
%         [L_{1d}, ..., L_{nd}]


[ind_r, ind_c, vals] = find(Left_CLtemp_betab);

sz_CL = size(Left_CLtemp_betab);
sz_LeftMat = [sz_CL(1)*nMid, floor(sz_CL(2)/nMid)];
BLCK_SIZE = sz_CL(2)/(nMid*length(CLtemp_betab));

ind_inBLC = mod(ind_c - 1, BLCK_SIZE) + 1; % indices in L_{id} 
ij_indx = floor((ind_c - 1)/BLCK_SIZE) + 1;
i_indx = floor((ij_indx - 1)/nMid) + 1; 
d_indx = mod(ij_indx - 1, nMid) + 1;


new_ind_c = ind_inBLC + (i_indx-1)*BLCK_SIZE; %mod(ind_c - 1, sz_LeftMat(2)) + 1;
new_ind_r = ind_r + (d_indx - 1)*sz_CL(1); 

LeftMat = sparse(new_ind_r, new_ind_c, vals, sz_LeftMat(1), sz_LeftMat(2));


% % % % STEP 2 % % % % 
% Right_CRtemp=> RightMat
% Left_CLtemp has the form
% [R_{11}; .. ; R_{1d} ; ... ; R_{n1}, .. R_{nd}] = [R_{jk}];
% we need [R_{11}, ..., R_{1d};
%          ... 
%         [R_{n1}, ..., R_{nd}]

[ind_r, ind_c, vals] = find(Right_CRtemp_alphaaT);

sz_CR = size(Right_CRtemp_alphaaT);
sz_RMat = [floor(sz_CR(1)/nMid), sz_CR(2)*nMid];

BLCK_SIZE = sz_CR(1)/(nMid*length(CLtemp_betab));

ind_inBLC = mod(ind_r - 1, BLCK_SIZE) + 1; % indices in R_{id} 
ij_indx = floor((ind_r - 1)/BLCK_SIZE) + 1;
i_indx = floor((ij_indx - 1)/nMid) + 1; 
d_indx = mod(ij_indx - 1, nMid) + 1;



new_ind_r = ind_inBLC + (i_indx-1)*BLCK_SIZE;
new_ind_c = ind_c + (d_indx - 1)*sz_CR(2); 


RightMat = sparse(new_ind_r, new_ind_c, vals, sz_RMat(1), sz_RMat(2));

% % % % STEP 3 % % % %
% Only the nMid DIAGONAL blocks of LeftMat*CMat*RightMat are wanted (sopvar % MMP, 09/10/2026
% document Sec. 5.4): block row d of LeftMat is [L_{1,d},...,L_{n,d}],      % MMP, 09/10/2026
% block column d of RightMat is [R_{1,d};...;R_{n,d}], and the sum over d   % MMP, 09/10/2026
% of Ld{d}*CMat*Rd{d} is formed directly, without the full product.         % MMP, 09/10/2026
right_BLCK_size = sz_CR(2);
left_BLCK_size  = sz_CL(1);
%Mask = kron(speye(nMid), ones(left_BLCK_size, right_BLCK_size));           % MMP, 09/10/2026 (was)
%LeftMult  = kron(ones(nMid, 1), speye(left_BLCK_size));                    % MMP, 09/10/2026 (was)
%RightMult = kron(ones(nMid, 1), speye(right_BLCK_size));                   % MMP, 09/10/2026 (was)
Ld = cell(1,nMid);   Rd = cell(1,nMid);                                     % MMP, 09/10/2026
for d = 1:nMid                                                              % MMP, 09/10/2026
    Ld{d} = LeftMat((d-1)*left_BLCK_size+(1:left_BLCK_size),:);             % MMP, 09/10/2026
    Rd{d} = RightMat(:,(d-1)*right_BLCK_size+(1:right_BLCK_size));          % MMP, 09/10/2026
end                                                                         % MMP, 09/10/2026


%%%%%%%%%%%%%%%%%%%%%%%%%% PART 4 %%%%%%%%%%%%%%%%%%%%%%%%%%
% Construct Cparams
for i = 1:numel(Cparams)
    
    Cparams{i} = sparse(size(CLtemp_betab{1},1),size(CRtemp_alphaa{1},1)); 
    for d = 1:nMid                                                          % MMP, 09/10/2026
        Cparams{i} = Cparams{i} + Ld{d}*CMat{i}*Rd{d};                      % MMP, 09/10/2026
    end                                                                     % MMP, 09/10/2026
%   CLCMCR = LeftMat*CMat{i}*RightMat;                                      % MMP, 09/10/2026 (was)
%   PROD2 = Mask.*CLCMCR;                                                   % MMP, 09/10/2026 (was)
%   PROD3 =LeftMult'*PROD2*RightMult;                                       % MMP, 09/10/2026 (was)
%   Cparams{i} = Cparams{i} + PROD3;                                        % MMP, 09/10/2026 (was)
end


%%%%%%%%%%%%%%%%%%%% PART 5 %%%%%%%%%%%%%%%%%%%%%%
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

% Drop the composed-basis degrees that hold no coefficient: CZL/CZR list    % MMP, 09/11/2026
% the degrees the integrals can REACH, not the ones attained. An empty      % MMP, 09/11/2026
% degree is the zero function, so the operator is unchanged, and a          % MMP, 09/11/2026
% redundant row in every later equality constraint is avoided.              % MMP, 09/11/2026
NLc  = prod([cellfun(@numel,CZL),1]);                                       % MMP, 09/11/2026
NRc  = prod([cellfun(@numel,CZR),1]);                                       % MMP, 09/11/2026
occL = false(NLc,1);                                                        % MMP, 09/11/2026
occR = false(NRc,1);                                                        % MMP, 09/11/2026
for ii = 1:numel(Cparams)                                                   % MMP, 09/11/2026
    % [] and scalar 0 are the zero-block shorthand and hold no degree.      % MMP, 09/11/2026
    if numel(Cparams{ii})>1                                                 % MMP, 09/11/2026
        [ri,ci] = find(Cparams{ii});                                        % MMP, 09/11/2026
        % Rows are (matrix row outer, monomial inner), columns likewise.    % MMP, 09/11/2026
        occL(mod(ri-1,NLc)+1) = true;                                       % MMP, 09/11/2026
        occR(mod(ci-1,NRc)+1) = true;                                       % MMP, 09/11/2026
    end                                                                     % MMP, 09/11/2026
end                                                                         % MMP, 09/11/2026
[CZL,CZR,keepL,keepR] = prune_zero_monomials(CZL,CZR,occL,occR);            % MMP, 09/11/2026
if numel(keepL)<NLc || numel(keepR)<NRc                                     % MMP, 09/11/2026
    rowK = reshape(keepL(:)+(0:Cdims(1)-1)*NLc,[],1);                       % MMP, 09/11/2026
    colK = reshape(keepR(:)+(0:Cdims(2)-1)*NRc,[],1);                       % MMP, 09/11/2026
    for ii = 1:numel(Cparams)                                               % MMP, 09/11/2026
        if numel(Cparams{ii})>1                                             % MMP, 09/11/2026
            Cparams{ii} = Cparams{ii}(rowK,colK);                           % MMP, 09/11/2026
        end                                                                 % MMP, 09/11/2026
    end                                                                     % MMP, 09/11/2026
end                                                                         % MMP, 09/11/2026

C = sopvar(Cparams,Cvars,CZL,CZR,Cdom,Cdims);
end

function Cparams_out = expandCommonFullIntegrals(Cparams_in,baseVars,finalVars)
% CPARAMS_OUT = EXPANDCOMMONFULLINTEGRALS(CPARAMS_IN,BASEVARS,FINALVARS)
% re-indexes the parameter cells from the shared variables BASEVARS of the
% composition to the shared variables FINALVARS of the result. A variable
% of FINALVARS not in BASEVARS is shared only because an input-only
% variable of B has the same name as an output-only variable of A; the
% operation along it is a full integral, stored as the lower (2) plus the
% upper (3) integral with a zero multiplier (1).
%
% INPUTS
% - Cparams_in: cell array of parameters indexed over BASEVARS;
% - baseVars:   cellstr, the shared variables the cells are indexed by;
% - finalVars:  cellstr, the shared variables of the result (a superset).
%
% OUTPUT
% - Cparams_out: 3 x ... x 3 cell array indexed over FINALVARS.

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
% IDX = PARAMMULTIIDX_SS(N) lists the multi-indices of the 3^N parameter
% cells in linear-index order: first variable changes fastest.
% INPUT   n: number of shared variables.
% OUTPUT  idx: 3^n x n array with entries 1, 2, 3.
if n == 0
    idx = zeros(1,0);
else
    idx = fliplr(dec2base(0:3^n-1,3,n)-'0') + 1;
end
end

function lin = paramMultiIdxToLinear_SS(idx)
% LIN = PARAMMULTIIDXTOLINEAR_SS(IDX) returns the linear cell index of one
% multi-index, lin = 1 + sum_i (idx(i)-1)*3^(i-1) (first variable fastest).
% INPUT   idx: 1 x n multi-index with entries 1, 2, 3.
% OUTPUT  lin: linear index into the 3 x ... x 3 cell array.
idx = idx(:).';
lin = 1;
for i = 1:numel(idx)
    lin = lin + (idx(i)-1)*3^(i-1);
end
end

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

        Ci = sparse(nE,nI);

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

dims2a = cellfun(@numel, Z2anew);
dims3a = cellfun(@numel, Z3anew);
dims3b = cellfun(@numel, Z3bnew);

% Vectorized over the whole m-by-n grid: every per-direction lookup depends % MMP, 09/10/2026
% on (i,j) only through IL(i,:) and IR(j,:), so each is a linear-index      % MMP, 09/10/2026
% gather. (i,k2a,k3a) -> row and (j,k3b) -> col are bijections, so no two   % MMP, 09/10/2026
% (i,j) collide and 'sparse' sums no duplicates: exact, not approximate.    % MMP, 09/10/2026
iv = reshape(repmat((1:m).',1,n),[],1);                                     % MMP, 09/10/2026
jv = reshape(repmat(1:n,m,1),[],1);                                         % MMP, 09/10/2026

% Scalar contribution from full integration over s2b.                       % MMP, 09/10/2026
VV = ones(m*n,1);                                                           % MMP, 09/10/2026
for k = 1:n2b                                                               % MMP, 09/10/2026
    Ck = C2b{k};                                                            % MMP, 09/10/2026
    % reshape, not just indexing: when Ck is a ROW vector MATLAB returns a  % MMP, 09/10/2026
    % row for a column index, and 'VV .* row' would implicitly expand to an % MMP, 09/10/2026
    % m*n by m*n matrix instead of failing loudly.                          % MMP, 09/10/2026
    VV = VV .* reshape(Ck(IL(iv,L2b(k)) + (IR(jv,R2b(k))-1)*size(Ck,1)),[],1);% MMP, 09/10/2026
end                                                                         % MMP, 09/10/2026

% Condensed basis indices for the s2a, s3a and s3b products.                % MMP, 09/10/2026
k2a = condenseIdx(idx2a,IL,IR,iv,jv,L2a,R2a,dims2a,n2a);                    % MMP, 09/10/2026
k3a = condenseIdx(idx3a,IL,IR,iv,jv,L3a,R3a,dims3a,n3a);                    % MMP, 09/10/2026
k3b = condenseIdx(idx3b,IL,IR,iv,jv,L3b,R3b,dims3b,n3b);                    % MMP, 09/10/2026

% Row: k3a inside block (i,k2a).  Column: k3b inside block j.               % MMP, 09/10/2026
II = k3a + (k2a + (iv-1)*nz2a - 1)*nz3a;                                    % MMP, 09/10/2026
JJ = k3b + (jv-1)*nz3b;                                                     % MMP, 09/10/2026

C3a = sparse(II, JJ, VV, m*nz2a*nz3a, n*nz3b);

% The condensed bases and the coefficient structure
Z2aout = Z2anew;
Z3bout = Z3bnew;

G3aout = struct();
G3aout.C = C3a;
G3aout.Z = Z3anew;

end
function kk = condenseIdx(T,IL,IR,iv,jv,Lk,Rk,dims,nk)                      % MMP, 09/10/2026
% KK = CONDENSEIDX(T,IL,IR,IV,JV,LK,RK,DIMS,NK) condenses the per-direction % MMP, 09/10/2026
% indices of one variable group into a single basis index, for every (i,j)  % MMP, 09/10/2026
% of the grid at once: the LAST direction runs fastest, so direction k      % MMP, 09/10/2026
% carries stride prod(dims(k+1:end)); kk = 1 for nk == 0.                   % MMP, 09/10/2026
% INPUTS  T: cell of per-direction lookup tables (ZL_i + ZR_j -> index);
%         IL, IR: multi-index tables of the ordered bases; iv, jv: the
%         (i,j) pairs of the grid; Lk, Rk: the columns of IL, IR of the
%         group; dims: the condensed basis sizes; nk: number of directions.
% OUTPUT  kk: numel(iv) x 1 condensed basis index of each pair.
kk = ones(numel(iv),1);                                                     % MMP, 09/10/2026
if nk == 0                                                                  % MMP, 09/10/2026
    return                                                                  % MMP, 09/10/2026
end                                                                         % MMP, 09/10/2026
str = ones(1,nk);                                                           % MMP, 09/10/2026
for k = nk-1:-1:1                                                           % MMP, 09/10/2026
    str(k) = str(k+1)*dims(k+1);                                            % MMP, 09/10/2026
end                                                                         % MMP, 09/10/2026
for k = 1:nk                                                                % MMP, 09/10/2026
    Tk = T{k};                                                              % MMP, 09/10/2026
    a  = Tk(IL(iv,Lk(k)) + (IR(jv,Rk(k))-1)*size(Tk,1));                    % MMP, 09/10/2026
    kk = kk + (a(:)-1)*str(k);                                              % MMP, 09/10/2026
end                                                                         % MMP, 09/10/2026
end                                                                         % MMP, 09/10/2026


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function P = kronPermutationMatrix_SS(oldVars,newVars,dimsOld)
% P = KRONPERMUTATIONMATRIX_SS(OLDVARS,NEWVARS,DIMSOLD) returns the
% permutation matrix P with Z_old = P*Z_new for tensor monomial vectors
% ordered with the first variable outermost and the last variable fastest.
% INPUTS  oldVars, newVars: cellstr, the same variables in two orders;
%         dimsOld: numbers of monomials per variable, in the old order.
% OUTPUT  P: sparse prod(dimsOld) x prod(dimsOld) permutation matrix.

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
% IDX = KRONMULTIIDX_SS(DIMS) lists the multi-indices of a tensor basis in
% the PIETOOLS order: first variable outermost, last variable fastest.
% INPUT   dims: numbers of monomials per variable.
% OUTPUT  idx: prod(dims) x numel(dims) array of 1-based positions.

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
% LIN = KRONSUB2IND_SS(DIMS,IDX) returns the linear positions of the
% multi-indices IDX in the tensor basis, the inverse of kronMultiIdx_SS.
% INPUTS  dims: numbers of monomials per variable; idx: rows of 1-based
%         per-variable positions.
% OUTPUT  lin: size(idx,1) x 1 linear positions.

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
