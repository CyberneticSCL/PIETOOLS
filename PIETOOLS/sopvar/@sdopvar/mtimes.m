function C = mtimes(A,B)
% C = MTIMES(A,B) returns the product A*B of an 'sdopvar' object and a
% scalar, a matrix or a fixed 'sopvar' object (composition), affine in the
% decision variables of the 'sdopvar' factor.
%
% INPUTS
% - A, B:   the factors, one of which is an 'sdopvar'. The other is
%           - a numeric scalar (scales the operator),
%           - a numeric matrix of compatible size (multiplies the output or
%             input components),
%           - a fixed 'sopvar' (composition; the product stays affine in
%             the decision variables), or
%           - a 'dpvar' (routed to 'dpvar_op_copvar').
%           A product of two 'sdopvar' objects is refused, being quadratic
%           in the decision variables;
%
% OUTPUTS
% - C:      'sdopvar' object representing the product A*B;
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PIETOOLS - mtimes
%
% Copyright (C) 2026 PIETOOLS Team
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
% MMP, 09/22/2026: Assembly loop of 'mtimesSdopvarAndSopvar' (29% of a
%                  three-variable composition): (1) the decision coefficient
%                  is accumulated as triplets and built with one 'sparse'
%                  per gamma, instead of a sparse '+' per (j,k) that rebuilt
%                  the structure 19683 times (54% of the routine); (2)
%                  vec(L*MR) = (MR' kron I)*vec(L) is applied as a reshape,
%                  and kron(MR,I) is built directly rather than as the
%                  transpose of kron(MR.',I); (3) a gamma cell whose middle
%                  factor is all zero is skipped before the kron (93.2% of
%                  the cells at three variables). The product is unchanged
%                  to one ulp (relative difference 6.2e-17) but not bit for
%                  bit: 'sparse' sums duplicates in another order, so 0.6%
%                  more entries, at eps, are stored (not pruned: a tolerance
%                  could drop a small coefficient). Composition 3.1x faster
%                  at three variables, 1.5x at two, unchanged at one.
% MMP, 09/11/2026: The composed monomial bases are pruned before C is built,
%                  as in '@sopvar/mtimes': CZL/CZR hold every degree the
%                  semiseparable integrals can reach (sopvar.pdf S4.1) and
%                  most hold no coefficient (31 carried against 13 occupied
%                  per direction on the 2-D heat equation); the empty
%                  degrees are columns of every 'B' and redundant equality
%                  rows later. A degree is dropped only when neither 'A' nor
%                  any row of 'B' has a nonzero there, so the operator is
%                  unchanged.
% MMP, 09/09/2026: B's output side is aligned onto A.vars.in by NAME, as in
%                  '@sopvar/mtimes'. The lists are stored in canonical order
%                  (vars.out = [S2,S3], vars.in = [S3,S1]), so a legal
%                  composition such as A'*A presents them in different
%                  orders, and the old guard refused every composition whose
%                  intermediate space had both a pass-through and a
%                  non-pass-through variable.
% MMP, 09/29/2026: A dpvar factor is routed at the top to 'dpvar_op_copvar'
%                  (legacy @dopvar/mtimes semantics): a constant dpvar
%                  scales, one with decision variables is refused (the
%                  product is quadratic), a matrix multiplies the output or
%                  input components. Before, dpvar*P left C unassigned.
% MMP, 11/06/2026: Comment pass (maintainer, 11/06/2026): what-line, log
%                  entries condensed to defect and fix, step comments through
%                  the composition, what/INPUT/OUTPUT headers for the local
%                  functions, commented-out calls of the retired
%                  'MatrixMultiply' and a commented-out reshape/permute
%                  construction of PL/PR removed; code unchanged.
% MMP, 11/09/2026: Inlining pass (maintainer, 11/09/2026: locals that are
%                  very short and called once go into the body):
%                  parameterMultiIndexToLinear into
%                  expandCommonFullIntegrals_SD, kronSub2ind_SS into
%                  kronPermutationMatrix_SS, each with its header as the
%                  comment of the block.
% DJ, 11/10/2026: The decision data of the left factor is no longer
%                  mapped onto all nRowsC*nColsL columns before the
%                  assembly. Only the columns reached from a
%                  nonzero column of B are formed, and the rows of
%                  kron(MR,I) on those are gathered from MR.'. Same
%                  triplets as before. Replaces the kron(MR,speye(nRowsC))
%                  form of the decision coefficient (MMP, 09/22/2026, (2)).


% dpvar factor: legacy semantics, before any work; the composition below    % MMP, 09/29/2026
% is untouched.                                                             % MMP, 09/29/2026
if isa(A,'dpvar') || isa(B,'dpvar')                                         % MMP, 09/29/2026
    C = dpvar_op_copvar('mtimes',A,B);                                      % MMP, 09/29/2026
    return                                                                  % MMP, 09/29/2026
end                                                                         % MMP, 09/29/2026
if isa(A,'sdopvar') 
    if isnumeric(B)
        % Scalar: scale every coefficient
        if isscalar(B)
            C = A;
            for ii=1:numel(C.params.A)
                C.params.A{ii} = C.params.A{ii}*B;
                C.params.B{ii} = C.params.B{ii}*B;
            end
            return
        end
        if A.dims(2)~=size(B,1)
            error('Inner dimensions do not match.');
        end
        % Matrix on the input side: every kernel is right-multiplied by
        % (B kron I_nR), the columns being (component outer, ZR inner)
        nL = prod(cellfun(@numel,A.ZL));
        nR = prod(cellfun(@numel,A.ZR));
        
        L = speye(A.dims(1)*nL);
        R = kron(B,speye(nR));
        
        params = transformParams(A.params,L,R);
        C = sdopvar(params,A.vars,A.Zd,A.ZL,A.ZR,A.dom,[A.dims(1),size(B,2)]);
    elseif isa(B,'sopvar')
        % Composition with a fixed operator on the right
        C = mtimesSdopvarAndSopvar(A,B);
    elseif isa(B,'sdopvar')
        error(['Composition of two sdopvar objects is not supported because ' ...
               'it is generally quadratic in the decision variables.']);
    end
elseif isa(B,'sdopvar')
    if isnumeric(A)
        if isscalar(A)
            C = B;
            for ii=1:numel(C.params.A)
                C.params.A{ii} = A*C.params.A{ii};
                C.params.B{ii} = A*C.params.B{ii};
            end
            return
        end
        if size(A,2)~=B.dims(1)
            error('Inner dimensions do not match.');
        end
        % Matrix on the output side: every kernel is left-multiplied by
        % (A kron I_nL), the rows being (component outer, ZL inner)
        nL = prod(cellfun(@numel,B.ZL));
        nR = prod(cellfun(@numel,B.ZR));
        
        L = kron(A,speye(nL));
        R = speye(B.dims(2)*nR);
        
        params = transformParams(B.params,L,R);
        C = sdopvar(params,B.vars,B.Zd,B.ZL,B.ZR,B.dom,[size(A,1),B.dims(2)]);
    elseif isa(A,'sopvar')
        % Fixed operator on the left: A*B = (B'*A')', computed by the
        % branch above
        Ct = mtimes(B',A');
        C = Ct';
    end
end
end
function C = mtimesSdopvarAndSopvar(A,B)
% C = MTIMESSDOPVARANDSOPVAR(A,B) composes the 'sdopvar' A with the fixed
% 'sopvar' B (sopvar.pdf sec. 6). The intermediate variables
% A.vars.in = B.vars.out are partitioned into those passing through both
% operators (vs3a), through A only (vs2a), through B only (vs3b) and
% through neither (vs2b, integrated over the full domain); the product of
% the two monomial bases is integrated over vs2b ('int_2b'), the one-sided
% integrals of A over vs2a and of B over vs3b are applied ('int_monomial'),
% the semiseparable composition in vs3a is tabulated ('int_semisep'), and
% the result is assembled cell by cell, affine in the decision variables.
% INPUTS  A: 'sdopvar' object; B: 'sopvar' object with B.dims(1)==A.dims(2)
%         and the same intermediate variables and domains.
% OUTPUT  C: 'sdopvar' object A*B on the decision variables A.Zd.
if A.dims(2) ~= B.dims(1)
    error('mtimes:innerDimensions', ...
        'The output dimension of B must equal the input dimension of A.');
end

% Legality is a question about the variable SETS. The body below contracts  % MMP, 09/09/2026
% A.ZR against B.ZL positionally, so the two lists must then agree in       % MMP, 09/09/2026
% ORDER -- and they need not, because each operand is stored in its own     % MMP, 09/09/2026
% canonical order: for A'*A the adjoint's vars.in is [S3,S2] while A's      % MMP, 09/09/2026
% vars.out is [S2,S3]. So B is aligned onto A.vars.in by name below         % MMP, 09/09/2026
% rather than the composition being refused.                                % MMP, 09/09/2026
if ~isempty(setxor(A.vars.in,B.vars.out))                                   % MMP, 09/09/2026
%if ~isequal(A.vars.in(:),B.vars.out(:))                                    % MMP, 09/09/2026 (was)
    error('mtimes:intermediateVariables', ...
        'B output variable names differ from A input variable names.');
end

% Permute B's output side onto A's input order: reorder B.ZL, B.dom.out     % MMP, 09/09/2026
% and the ROW index of every parameter, the row being (matrix row outer,    % MMP, 09/09/2026
% ZL monomial inner). Done in place on this local copy, NOT through the     % MMP, 09/09/2026
% constructor, which would canonicalize B straight back. The permutation    % MMP, 09/09/2026
% touches only rows, so the canonical multiplier form -- a condition on     % MMP, 09/09/2026
% ZR columns -- is untouched. B is an sopvar here, so its parameters are    % MMP, 09/09/2026
% plain matrices; A's decision-variable coefficients are not touched.       % MMP, 09/09/2026
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

if ~isequal(A.dom.in,B.dom.out)
    error('mtimes:intermediateDomains', ...
        'B output domains must match A input domains.');
end

% Partition the intermediate variables A.vars.in = B.vars.out: S1 of A is
% integrated out by A and S3 of A passes through A; S2 of B is produced by
% B and S3 of B passes through B. Hence
%   vs3a = pass through both (a pass-through variable of C),
%   vs3b = pass through B, integrated out by A,
%   vs2a = pass through A, produced by B,
%   vs2b = integrated out by A, produced by B: a full integral.
AvarsS1 = setdiff(A.vars.in,A.vars.out);
AvarsS3 = intersect(A.vars.in,A.vars.out);

BvarsS2 = setdiff(B.vars.out,B.vars.in);
BvarsS3 = intersect(B.vars.in,B.vars.out);

vs3a = intersect(AvarsS3,BvarsS3);
vs3b = setdiff(BvarsS3,vs3a);
vs2a = setdiff(AvarsS3,vs3a);
vs2b = setdiff(AvarsS1,vs3b);

% Final component dimensions.
Cdims = [A.dims(1),B.dims(2)];

% The cells of C are indexed by vs3a, one direction each; further
% pass-through variables of C are added by 'expandCommonFullIntegrals_SD'.
n3a = numel(vs3a);

if n3a == 0
    Csize = [1,1];
elseif n3a == 1
    Csize = [3,1];
else
    Csize = 3*ones(1,n3a);
end

nalpha = numel(BvarsS3);
nbeta  = numel(AvarsS3);

if numel(A.params.A) ~= 3^nbeta || ...
        numel(A.params.B) ~= 3^nbeta
    error('mtimes:malformedSdopvar', ...
        'A must have one affine parameter for every S3 index.');
end

if numel(B.params) ~= 3^nalpha
    error('mtimes:malformedSopvar', ...
        'B must have one parameter for every S3 index.');
end

% Parameter multi-indices. The first variable changes fastest.
alphaIdx = parameterMultiIndex(nalpha);
betaIdx  = parameterMultiIndex(nbeta);

% Split the cell index beta of A (over AvarsS3 = vs2a u vs3a) into beta_a
% over vs2a and beta_b over vs3a, and alpha of B (over BvarsS3 = vs3a u
% vs3b) into alpha_a over vs3a and alpha_b over vs3b.
[found,colsbetaa] = ismember(vs2a,AvarsS3);
if any(~found)
    error('mtimes:betaAIndexing', ...
        'vs2a was not found in the common variables of A.');
end

[found,colsbetab] = ismember(vs3a,AvarsS3);
if any(~found)
    error('mtimes:betaBIndexing', ...
        'vs3a was not found in the common variables of A.');
end

% Split alpha into [alpha_a,alpha_b].
[found,colsalphaa] = ismember(vs3a,BvarsS3);
if any(~found)
    error('mtimes:alphaAIndexing', ...
        'vs3a was not found in the common variables of B.');
end

[found,colsalphab] = ismember(vs3b,BvarsS3);
if any(~found)
    error('mtimes:alphaBIndexing', ...
        'vs3b was not found in the common variables of B.');
end

% The distinct sub-indices, in order of first appearance
beta_a = unique( ...
    betaIdx(:,colsbetaa),'rows','stable');

beta_b = unique( ...
    betaIdx(:,colsbetab),'rows','stable');

alpha_a = unique( ...
    alphaIdx(:,colsalphaa),'rows','stable');

alpha_b = unique( ...
    alphaIdx(:,colsalphab),'rows','stable');

% Package the sdopvar coefficients as affine structures.
Acoeff = cell(size(A.params.A));

for k = 1:numel(Acoeff)
    Acoeff{k} = struct( ...
        'A',A.params.A{k}, ...
        'B',A.params.B{k});
end

% Regroup the cells as 2-D arrays CA{beta_a,beta_b} and CB{alpha_a,alpha_b}.
CA = regroupParameterCells( ...
    Acoeff,betaIdx,colsbetaa,colsbetab);

CB = regroupParameterCells( ...
    B.params,alphaIdx,colsalphaa,colsalphab);

% vs2b consists of output-only variables of B.
dom2b = selectDomain( ...
    vs2b,B.vars.out,B.dom.out,'vs2b');

% vs2a consists of common variables of A.
dom2a = selectDomain( ...
    vs2a,A.vars.in,A.dom.in,'vs2a');

% vs3a and vs3b are common variables of B.
dom3a = selectDomain( ...
    vs3a,B.vars.in,B.dom.in,'vs3a');

dom3b = selectDomain( ...
    vs3b,B.vars.in,B.dom.in,'vs3b');

p = prod(cellfun(@numel,A.ZR));
q = prod(cellfun(@numel,B.ZL));

% Integrate the product A.ZR(s)*B.ZL(s)' over vs2b and condense the product
% monomials in vs2a, vs3a and vs3b: the bridge between the two kernels,
% factored as PL*(I kron Zhat_s2a')*G_s3a(s3a)*(I kron Zhat_s3b)*PR
[Zhat_s2a,G_s3a,Zhat_s3b,PL,PR] = int_2b( ...
    A.ZR,B.ZL,A.vars.in, ...
    vs2a,vs2b,vs3a,vs3b,dom2b);

% In a vs3b direction B integrates its input variable th over [a,s] (2) or
% [s,b] (3) and A integrates s over [a,b]; exchanging the order, s runs
% over [th,b] or [a,th], so the monomials in s take the opposite one-sided
% integral: 2 <-> 3
alpha_b_int = alpha_b;
alpha_b_int(alpha_b == 2) = 3;
alpha_b_int(alpha_b == 3) = 2;

% One-sided integrals of A over vs2a (index beta_a) and of B over vs3b
% (index alpha_b), each as a coefficient matrix onto an expanded basis
[Cs2a_betaa,KZhat_s2a] = ...
    int_monomial(Zhat_s2a,beta_a,dom2a);

[Cs3b_alphab,KZhat_s3b] = ...
    int_monomial(Zhat_s3b,alpha_b_int,dom3b);

% Composition in the vs3a directions: cell gamma of C receives from the
% cells (beta_b of A, alpha_a of B) the coefficient C_gam_alp_beta
% {gamma,beta_b,alpha_a} on the bases barZ_3aL (left) and barZ_3aR (right)
[C_gam_alp_beta,barZ_3aL,barZ_3aR] = ...
    int_semisep(G_s3a,beta_b,alpha_a,dom3a,Csize);

% Left factor: the columns of every cell of A (component outer, ZR inner)
% are permuted by PL into the variable order [vs2a,vs2b,vs3a,vs3b] of
% 'int_2b'
nAL = prod(cellfun(@numel,A.ZL));
nAR = prod(cellfun(@numel,A.ZR));

nRowsA = A.dims(1)*nAL;
nColsA = A.dims(2)*nAR;

if ~isequal(size(PL),[nAR,nAR])
    error('mtimes:PLDimensions', ...
        'PL has dimensions inconsistent with A.ZR.');
end

PLbig = kron(speye(A.dims(2)),PL);
Ileft = speye(nRowsA);

nBetaA = size(beta_a,1);
nBetaB = size(beta_b,1);

CA_PL = cell(nBetaA,nBetaB);

for ia = 1:nBetaA
    for ib = 1:nBetaB
        param = CA{ia,ib};

        if numel(param.A) ~= nRowsA*nColsA || ...
                size(param.B,2) ~= nRowsA*nColsA
            error('mtimes:affineCoefficientDimensions', ...
                'An affine coefficient of A has inconsistent dimensions.');
        end

        CA_PL{ia,ib} = lrmultiply( ...
            Ileft,param,PLbig);
    end
end


qL = size(PLbig,2);

% Apply the vs2a integrals of A on the column side of each cell, then move
% their monomials (in the vs2a variables, output variables of C) from the
% column side onto the row basis of A: ZLtemp over varLtemp
CB_s2a = cell(nBetaA,1);

for ia = 1:nBetaA
    CB_s2a{ia} = kron( ...
        speye(qL),Cs2a_betaa{ia}.');
end

[varLtemp,ZLtemp,CLtemp_betaab] = ...
    leftShiftMonomials_SD( ...
        A.vars.out,A.ZL,CA_PL, ...
        vs2a,KZhat_s2a,CB_s2a);

% C has no cell index in vs2a: sum over beta_a
CLtemp = cell(1,nBetaB);

for ib = 1:nBetaB
    CLtemp{ib} = CLtemp_betaab{1,ib};

    for ia = 2:nBetaA
        CLtemp{ib}.A = CLtemp{ib}.A ...
                     + CLtemp_betaab{ia,ib}.A;

        CLtemp{ib}.B = CLtemp{ib}.B ...
                     + CLtemp_betaab{ia,ib}.B;
    end
end


% Right factor, built on the transposed cells of B: the rows of every cell
% (component outer, ZL inner) are permuted by PR into int_2b's order, the
% vs3b integrals of B applied, and their monomials (in the vs3b variables,
% input variables of C) moved onto the column basis of B: ZRtemp over
% varRtemp; then sum over alpha_b
nAlphaA = size(alpha_a,1);
nAlphaB = size(alpha_b,1);

CBt = cell(nAlphaB,nAlphaA);

for ia = 1:nAlphaA
    for ib = 1:nAlphaB
        CBt{ib,ia} = CB{ia,ib}.';
    end
end

nBL = prod(cellfun(@numel,B.ZL));

if ~isequal(size(PR),[nBL,nBL])
    error('mtimes:PRDimensions', ...
        'PR has dimensions inconsistent with B.ZL.');
end

PRbig = kron(speye(B.dims(1)),PR.');

CBt_PR = cell(size(CBt));

for k = 1:numel(CBt)
    CBt_PR{k} = CBt{k}*PRbig;
end
qR = size(PRbig,2);

CB_s3b = cell(nAlphaB,1);

for ib = 1:nAlphaB
    CB_s3b{ib} = kron( ...
        speye(qR),Cs3b_alphab{ib}.');
end
[varRtemp,ZRtemp,CRtemp_alphaba] = ...
    leftShiftMonomials_SS( ...
        B.vars.in,B.ZR,CBt_PR, ...
        vs3b,KZhat_s3b,CB_s3b);
CRtemp = cell(nAlphaA,1);

for ia = 1:nAlphaA
    CRtemp{ia} = CRtemp_alphaba{1,ia};

    for ib = 2:nAlphaB
        CRtemp{ia} = CRtemp{ia} ...
                   + CRtemp_alphaba{ib,ia};
    end
end

% The vs3a monomials of the semiseparable composition, barZ_3aL on the left
% and barZ_3aR on the right, are moved onto the row basis of the left factor
% and the column basis of the right factor (identity coefficient, the bases
% only merge): CZL over Cvarsout, CZR over Cvarsin
nZLtemp = prod(cellfun(@numel,ZLtemp));
if isempty(ZLtemp)
    nZLtemp = 1;
end

nRowsL = A.dims(1)*nZLtemp;
nCoefL = numel(CLtemp{1}.A);

if mod(nCoefL,nRowsL) ~= 0
    error('mtimes:leftCoefficientDimensions', ...
        'CLtemp cannot be reshaped consistently with ZLtemp.');
end

qLbar = nCoefL/nRowsL;

nBarL = prod(cellfun(@numel,barZ_3aL));
if isempty(barZ_3aL)
    nBarL = 1;
end

CB_barL = {speye(qLbar*nBarL)};

% Field B of CLtemp_betab is returned unmapped, and TLtemp{1} is            % DJ (CC), 11/10/2026
% applied to its nonzero columns only, below, not to all of them.
%[Cvarsout,CZL,CLtemp_betab] = leftShiftMonomials_SD( ...                   % DJ (CC), 11/10/2026 (was)
%    varLtemp,ZLtemp,CLtemp,vs3a,barZ_3aL,CB_barL);                         
[Cvarsout,CZL,CLtemp_betab,TLtemp] = leftShiftMonomials_SD( ...             % DJ (CC), 11/10/2026
     varLtemp,ZLtemp,CLtemp,vs3a,barZ_3aL,CB_barL);                                             
CRtemp_row = CRtemp.';
qRbar = size(CRtemp_row{1},2);

nBarR = prod(cellfun(@numel,barZ_3aR));
if isempty(barZ_3aR)
    nBarR = 1;
end

CB_barR = {speye(qRbar*nBarR)};
[Cvarsin,CZR,CRtemp_alphaa] = ...
    leftShiftMonomials_SS( ...
        varRtemp,ZRtemp,CRtemp_row, ...
        vs3a,barZ_3aR,CB_barR);

nCZL = prod(cellfun(@numel,CZL));
nCZR = prod(cellfun(@numel,CZR));

if isempty(CZL)
    nCZL = 1;
end
if isempty(CZR)
    nCZR = 1;
end

nRowsC = A.dims(1)*nCZL;
nColsC = B.dims(2)*nCZR;
nDecision = numel(A.Zd);
nCoefC = nRowsC*nColsC;

Cparams.A = repmat( ...
    {sparse(nCoefC,1)},Csize);

Cparams.B = repmat( ...
    {sparse(nDecision,nCoefC)},Csize);

nMid = A.dims(2);

% For B_L = B*TLtemp{1}.', only the columns that TLtemp{1} reaches from a   % DJ (CC), 11/10/2026
% nonzero column of B can hold data. Keep those as BLc{j}, the position 
% (c,l) of each in crow, lrow: column c+nRowsC*(l-1) of B_L.
BLc  = cell(size(beta_b,1),1);                                              % DJ (CC), 11/10/2026 (start)
lrow = cell(size(BLc));
crow = cell(size(BLc));
for j = 1:size(beta_b,1)
    nzB = find(any(CLtemp_betab{j}.B,1));
    [ti,tj,tv] = find(TLtemp{1}(:,nzB));
    [rr,~,ri] = unique(ti(:));
    BLc{j} = CLtemp_betab{j}.B(:,nzB) ...
           *sparse(ri,tj(:),tv(:),numel(rr),numel(nzB)).';
    lrow{j} = ceil(rr/nRowsC);
    crow{j} = rr-nRowsC*(lrow{j}-1);
end                                                                         % DJ (CC), 11/10/2026 (end)

% Assemble each cell gamma of C as the sum over (beta_b, alpha_a) of
% L_j * (I_nMid kron M_{gamma,j,k}) * R_k, with L_j the affine left factor,
% M the vs3a coefficient table and R_k the right factor
for igamma = 1:numel(Cparams.A)

    Aaccum = sparse(nCoefC,1);
    % The decision coefficient is accumulated as triplets and built with    % MMP, 09/22/2026
    % one 'sparse' per gamma ('sparse' sums duplicate subscripts); a sparse % MMP, 09/22/2026
    % '+' per (j,k) rebuilt the structure on every addition.                % MMP, 09/22/2026
    nAcc = 0;                                                               % MMP, 09/22/2026
    accI = cell(size(beta_b,1)*size(alpha_a,1),1);                          % MMP, 09/22/2026
    accJ = cell(size(accI));                                                % MMP, 09/22/2026
    accV = cell(size(accI));                                                % MMP, 09/22/2026
%   Baccum = sparse(nDecision,nCoefC);                                      % MMP, 09/22/2026 (was)

    for j = 1:size(beta_b,1)

        Lparam = CLtemp_betab{j};

        if mod(numel(Lparam.A),nRowsC) ~= 0
            error('mtimes:leftAssemblyDimensions', ...
                'The left affine coefficient has inconsistent dimensions.');
        end

        nColsL = numel(Lparam.A)/nRowsC;
        % The left constant term as a matrix, hoisted out of the k loop: it % MMP, 09/22/2026
        % does not depend on k, and with it the A accumulation can use      % MMP, 09/22/2026
        % vec(L*MR) = (MR' kron I)*vec(L) directly and never form the kron. % MMP, 09/22/2026
        Lmat = reshape(Lparam.A,nRowsC,nColsL);                             % MMP, 09/22/2026

%       if size(Lparam.B,1) ~= nDecision || ...
%               size(Lparam.B,2) ~= numel(Lparam.A)                         % DJ (CC), 11/10/2026 (was)
        if size(Lparam.B,1) ~= nDecision || ...
                size(Lparam.B,2) ~= size(TLtemp{1},2) || ...
                size(TLtemp{1},1) ~= numel(Lparam.A)                        % DJ (CC), 11/10/2026
            error('mtimes:leftAffineDimensions', ...
                'The left affine decision data has inconsistent dimensions.');
        end

        for k = 1:size(alpha_a,1)

            M = C_gam_alp_beta{igamma,j,k};

            if nMid ~= 1
                M = kron(speye(nMid),M);
            end

            R = CRtemp_alphaa{k}.';

            if size(M,1) ~= nColsL
                error('mtimes:leftMiddleDimensions', ...
                    'The left and middle factors have incompatible dimensions.');
            end

            if size(M,2) ~= size(R,1)
                error('mtimes:middleRightDimensions', ...
                    'The middle and right factors have incompatible dimensions.');
            end

            MR = M*R;

            if size(MR,2) ~= nColsC
                error('mtimes:outputCoefficientDimensions', ...
                    'The assembled right factor has an unexpected size.');
            end

            % An all-zero middle factor contributes nothing: skip the kron  % MMP, 09/22/2026
            % and the multiply (93.2% of the cells at three variables).      % MMP, 09/22/2026
            % Every dimension check above still runs.                       % MMP, 09/22/2026
            if nnz(MR)==0                                                   % MMP, 09/22/2026
                continue                                                    % MMP, 09/22/2026
            end                                                             % MMP, 09/22/2026

            % vec(L*MR) = (MR' kron I)*vec(L), applied as a reshape rather  % MMP, 09/22/2026
            % than by forming the Kronecker product.                        % MMP, 09/22/2026
            Aaccum = Aaccum + reshape(Lmat*MR,[],1);                        % MMP, 09/22/2026

            % kron(MR,I) is kron(MR.',I).', so it is built in that form     % MMP, 09/22/2026
            % directly and the transpose of a q-column sparse is avoided.   % MMP, 09/22/2026
%           [bi,bj,bv] = find(Lparam.B*kron(MR,speye(nRowsC)));             % DJ (CC), 11/10/2026 (was)
            % Compute the product B_L*kron(MR,I) only for the retained      % DJ (CC), 11/10/2026
            % columns. The kronecker product kron(MR,I) is not formed,
            % row (c,l) of kron(MR,I) is row l of MR at column 
            % c+nRowsC*(c'-1), gathered from MR.';
            if isempty(lrow{j})                                             % DJ (CC), 11/10/2026 (start)
                continue
            end
            MRt = MR.';
            [cc,pp,vv] = find(MRt(:,lrow{j}));
            % find returns rows if size(MR,2) = 1: force columns
            cc = cc(:);  pp = pp(:);  vv = vv(:);
            [bi,bj,bv] = find(BLc{j}*sparse(pp,crow{j}(pp)+nRowsC*(cc-1), ...
                                            vv,numel(lrow{j}),nCoefC));     % DJ (CC), 11/10/2026 (end)
            nAcc = nAcc+1;                                                  % MMP, 09/22/2026
            accI{nAcc} = bi(:);     accJ{nAcc} = bj(:);                     % MMP, 09/22/2026
            accV{nAcc} = bv(:);                                             % MMP, 09/22/2026
%           T = kron(MR.',speye(nRowsC));                                   % MMP, 09/22/2026 (was)
%           Aaccum = Aaccum + T*Lparam.A;                                   % MMP, 09/22/2026 (was)
%           Baccum = Baccum + Lparam.B*T.';                                 % MMP, 09/22/2026 (was)
        end
    end

    Cparams.A{igamma} = Aaccum;
    Cparams.B{igamma} = sparse(vertcat(accI{1:nAcc}), ...                   % MMP, 09/22/2026
                               vertcat(accJ{1:nAcc}), ...                   % MMP, 09/22/2026
                               vertcat(accV{1:nAcc}),nDecision,nCoefC);     % MMP, 09/22/2026
%   Cparams.B{igamma} = Baccum;                                             % MMP, 09/22/2026 (was)
end
% The variables of C: Cvarsout from the output side of A, Cvarsin from the
% input side of B, with their domains
[foundIn,idxIn] = ismember(Cvarsin,B.vars.in);
[foundOut,idxOut] = ismember(Cvarsout,A.vars.out);

if any(~foundIn)
    error('mtimes:finalInputVariables', ...
        'Final input variables do not match the input variables of B.');
end

if any(~foundOut)
    error('mtimes:finalOutputVariables', ...
        'Final output variables do not match the output variables of A.');
end

% A variable produced by A (its S2) and consumed by B (its S1) is a
% pass-through of C that neither factor passed through: its two domains
% must agree, and C carries a full integral in that direction
CvarsS3 = intersect(Cvarsin,Cvarsout);
extraCommon = setdiff(CvarsS3,vs3a);

if ~isempty(extraCommon)
    [foundExtraIn,idxExtraIn] = ...
        ismember(extraCommon,B.vars.in);

    [foundExtraOut,idxExtraOut] = ...
        ismember(extraCommon,A.vars.out);

    if any(~foundExtraIn) || any(~foundExtraOut)
        error('mtimes:extraCommonVariables', ...
            'An additional common variable could not be located.');
    end

    if ~isequal( ...
            B.dom.in(idxExtraIn,:), ...
            A.dom.out(idxExtraOut,:))
        error('mtimes:extraCommonDomains', ...
            'Matching final common variables must have matching domains.');
    end
end

% Expand the cell array of C from vs3a to all its pass-through variables
% (full integral = lower + upper in each extra direction)
Cparams = expandCommonFullIntegrals_SD( ...
    Cparams,vs3a,CvarsS3);

Cvars = struct( ...
    'in',{Cvarsin}, ...
    'out',{Cvarsout});

Cdom = struct( ...
    'in',B.dom.in(idxIn,:), ...
    'out',A.dom.out(idxOut,:));

% Drop the composed-basis degrees that hold no coefficient, as in           % MMP, 09/11/2026
% '@sopvar/mtimes' ('prune_zero_monomials'); here that also narrows the     % MMP, 09/11/2026
% column axis of every 'B'. A degree is empty only when neither 'A' nor     % MMP, 09/11/2026
% any row of 'B' has a nonzero there, so the operator is unchanged. 'find'  % MMP, 09/11/2026
% with two outputs reads the nonzero columns of 'B' without densifying.     % MMP, 09/11/2026
occL = false(nCZL,1);                                                       % MMP, 09/11/2026
occR = false(nCZR,1);                                                       % MMP, 09/11/2026
for ii = 1:numel(Cparams.A)                                                 % MMP, 09/11/2026
    if numel(Cparams.A{ii})~=nCoefC                                         % MMP, 09/11/2026
        continue                                                            % MMP, 09/11/2026
    end                                                                     % MMP, 09/11/2026
    [~,cB] = find(Cparams.B{ii});                                           % MMP, 09/11/2026
    % Positions into vec(C) of an nRowsC x nColsC block, column major.      % MMP, 09/11/2026
    lin = [find(Cparams.A{ii}); cB(:)];                                     % MMP, 09/11/2026
    rr  = mod(lin-1,nRowsC);                                                % MMP, 09/11/2026
    cc  = floor((lin-1)/nRowsC);                                            % MMP, 09/11/2026
    % Rows are (matrix row outer, monomial inner), columns likewise.        % MMP, 09/11/2026
    occL(mod(rr,nCZL)+1) = true;                                            % MMP, 09/11/2026
    occR(mod(cc,nCZR)+1) = true;                                            % MMP, 09/11/2026
end                                                                         % MMP, 09/11/2026
[CZL,CZR,keepL,keepR] = prune_zero_monomials(CZL,CZR,occL,occR);            % MMP, 09/11/2026
if numel(keepL)<nCZL || numel(keepR)<nCZR                                   % MMP, 09/11/2026
    rowK = reshape(keepL(:)+(0:A.dims(1)-1)*nCZL,[],1);                     % MMP, 09/11/2026
    colK = reshape(keepR(:)+(0:B.dims(2)-1)*nCZR,[],1);                     % MMP, 09/11/2026
    vecK = reshape(rowK(:)+(colK(:).'-1)*nRowsC,[],1);                      % MMP, 09/11/2026
    for ii = 1:numel(Cparams.A)                                             % MMP, 09/11/2026
        if numel(Cparams.A{ii})~=nCoefC                                     % MMP, 09/11/2026
            continue                                                        % MMP, 09/11/2026
        end                                                                 % MMP, 09/11/2026
        Cparams.A{ii} = Cparams.A{ii}(vecK);                                % MMP, 09/11/2026
        Cparams.B{ii} = Cparams.B{ii}(:,vecK);                              % MMP, 09/11/2026
    end                                                                     % MMP, 09/11/2026
end                                                                         % MMP, 09/11/2026

C = sdopvar( ...
    Cparams,Cvars,A.Zd,CZL,CZR,Cdom,Cdims);
end


% helper functions
function paramsOut = expandCommonFullIntegrals_SD( ...
    paramsIn,baseVars,finalVars)
% PARAMSOUT = EXPANDCOMMONFULLINTEGRALS_SD(PARAMSIN,BASEVARS,FINALVARS)
% expands the cell array of affine parameters from the pass-through
% variables BASEVARS to the superset FINALVARS: in each extra direction the
% operator is a full integral, so the multiplier cell is zero and the lower
% and upper cells both hold the source coefficient.
% INPUTS  paramsIn: struct with cell arrays A, B of 3^numel(baseVars)
%         cells; baseVars, finalVars: cellstr, baseVars a subset.
% OUTPUT  paramsOut: the same struct over 3^numel(finalVars) cells.

baseVars = baseVars(:).';
finalVars = finalVars(:).';

if isequal(baseVars,finalVars)
    paramsOut = paramsIn;
    return
end

if any(~ismember(baseVars,finalVars))
    error('expandCommonFullIntegrals_SD:variableMismatch', ...
        'Base variables must be contained in the final variables.');
end

nBase = numel(baseVars);
nFinal = numel(finalVars);

finalIdx = parameterMultiIndex(nFinal);

if nFinal == 0
    finalSize = [1,1];
elseif nFinal == 1
    finalSize = [3,1];
else
    finalSize = 3*ones(1,nFinal);
end

zeroA = sparse(size(paramsIn.A{1},1),1);
zeroB = sparse( ...
    size(paramsIn.B{1},1), ...
    size(paramsIn.B{1},2));

paramsOut.A = repmat({zeroA},finalSize);
paramsOut.B = repmat({zeroB},finalSize);

[isBaseVariable,basePosition] = ...
    ismember(finalVars,baseVars);

for k = 1:size(finalIdx,1)
    gamma = finalIdx(k,:);

    % An extra common variable represents a full integral:
    %
    %   full integral = lower integral + upper integral.
    %
    % Therefore its multiplier entry, gamma == 1, is zero.
    if any(~isBaseVariable & gamma == 1)
        continue
    end

    if nBase == 0
        sourceIndex = 1;
    else
        gammaBase = zeros(1,nBase);

        gammaBase(basePosition(isBaseVariable)) = ...
            gamma(isBaseVariable);

%       sourceIndex = parameterMultiIndexToLinear(gammaBase);               % MMP, 11/09/2026 (was)
        % the linear cell index 1 + sum_k (gammaBase(k)-1)*3^(k-1) of the   % MMP, 11/09/2026
        % multi-index gammaBase, the map of multiIndexToLinear below (the   % MMP, 11/09/2026
        % former local PARAMETERMULTIINDEXTOLINEAR)                         % MMP, 11/09/2026
        sourceIndex = 1;                                                    % MMP, 11/09/2026
        for kb = 1:numel(gammaBase)                                         % MMP, 11/09/2026
            sourceIndex = sourceIndex + (gammaBase(kb)-1)*3^(kb-1);         % MMP, 11/09/2026
        end                                                                 % MMP, 11/09/2026
    end

    paramsOut.A{k} = paramsIn.A{sourceIndex};
    paramsOut.B{k} = paramsIn.B{sourceIndex};
end
end
function idx = parameterMultiIndex(n)
% IDX = PARAMETERMULTIINDEX(N) lists the 3^N cell multi-indices in {1,2,3}^N
% as rows, in linear cell order (first variable fastest).
% INPUT   n: number of pass-through variables.
% OUTPUT  idx: 3^n x n array.

if n == 0
    idx = zeros(1,0);
else
    idx = fliplr( ...
        dec2base(0:3^n-1,3,n)-'0') + 1;
end
end

function selectedDom = selectDomain(selectedVars,allVars,allDom,label)
% SELECTEDDOM = SELECTDOMAIN(SELECTEDVARS,ALLVARS,ALLDOM,LABEL) returns the
% rows of the domain table ALLDOM (one row [a,b] per entry of ALLVARS) for
% the variables SELECTEDVARS, in that order; LABEL names the group in the
% error message when a variable is missing.
% INPUTS  selectedVars, allVars: cellstr; allDom: numel(allVars) x 2.
% OUTPUT  selectedDom: numel(selectedVars) x 2.

[found,idx] = ismember(selectedVars,allVars);

if any(~found)
    error('mtimes:missingVariables', ...
        '%s variables were not found in the expected variable list.',label);
end

selectedDom = allDom(idx,:);
end
function grouped = regroupParameterCells( ...
    paramCells,multiIdx,firstCols,secondCols)
% GROUPED = REGROUPPARAMETERCELLS(PARAMCELLS,MULTIIDX,FIRSTCOLS,SECONDCOLS)
% rearranges the linear cell array PARAMCELLS, whose k-th cell has the
% multi-index MULTIIDX(k,:), into a 2-D cell array indexed by the linear
% index of the sub-multi-index over the columns FIRSTCOLS and by that over
% SECONDCOLS.
% INPUTS  paramCells: 3^n cells; multiIdx: 3^n x n; firstCols, secondCols:
%         a partition of 1:n.
% OUTPUT  grouped: 3^numel(firstCols) x 3^numel(secondCols) cell array.

nFirst  = 3^numel(firstCols);
nSecond = 3^numel(secondCols);

grouped = cell(nFirst,nSecond);

for k = 1:numel(paramCells)
    i = multiIndexToLinear(multiIdx(k,firstCols));
    j = multiIndexToLinear(multiIdx(k,secondCols));

    grouped{i,j} = paramCells{k};
end
end

function idx = multiIndexToLinear(alpha)
% IDX = MULTIINDEXTOLINEAR(ALPHA) is the linear cell index of the
% multi-index ALPHA, 1 + sum_k (alpha(k)-1)*3^(k-1).
% INPUT   alpha: 1 x n multi-index in {1,2,3}^n.
% OUTPUT  idx: integer in 1..3^n.

idx = 1;

for k = 1:numel(alpha)
    idx = idx + (alpha(k)-1)*3^(k-1);
end
end
function paramsOut = transformParams(paramsIn,L,R)
% PARAMSOUT = TRANSFORMPARAMS(PARAMSIN,L,R) applies C -> L*C*R to the
% coefficient matrix of every cell, constant and decision parts, through
% 'lrmultiply'.
% INPUTS  paramsIn: struct with cell arrays A (vec of C) and B (nd x
%         numel(C)); L, R: matrices conformable with C.
% OUTPUT  paramsOut: struct of the same form.

paramsOut.A = cell(size(paramsIn.A));
paramsOut.B = cell(size(paramsIn.B));

for k = 1:numel(paramsIn.A)
    param = struct( ...
        'A',paramsIn.A{k}, ...
        'B',paramsIn.B{k});

    param = lrmultiply(L,param,R);

    paramsOut.A{k} = param.A;
    paramsOut.B{k} = param.B;
end
end
function [Cidx, Zint] = int_monomial(Z,idxAll,lims)
% [CIDX,ZINT] = INT_MONOMIAL(Z,IDXALL,LIMS) applies the one-sided integrals
% of a cell multi-index to a tensor monomial basis.
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

% Build the common expanded output basis, which must serve every row of
% idxAll
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
% [Z2AOUT,G3AOUT,Z3BOUT,PL,PR] = INT_2B(ZL,ZR,ZVAR,S2A,S2B,S3A,S3B,LIMS)
% integrates the product of two tensor monomial bases over the variables
% S2B and condenses the product monomials in the others: the factorization
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

% Monomial vectors are ordered with the first variable outermost and the
% last fastest, the opposite of MATLAB's first-dimension-fastest reshape,
% so the permutation is built from explicit multi-indices
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

% OUTPUTS
Z2aout = Z2anew;
Z3bout = Z3bnew;

G3aout = struct();
G3aout.C = C3a;
G3aout.Z = Z3anew;

end
function P = kronPermutationMatrix_SS(oldVars,newVars,dimsOld)
% P = KRONPERMUTATIONMATRIX_SS(OLDVARS,NEWVARS,DIMSOLD) returns the
% permutation matrix with Z_old = P*Z_new for tensor monomial vectors
% ordered with the first variable outermost and the last fastest.
% INPUTS  oldVars: cellstr; newVars: a permutation of oldVars; dimsOld:
%         number of monomials per variable, in the order of oldVars.
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
% oldLin = kronSub2ind_SS(dimsOld,idxOld);                                  % MMP, 11/09/2026 (was)
% oldLin: the linear position of each multi-index row of idxOld in the      % MMP, 11/09/2026
% tensor order of dimsOld, the inverse of kronMultiIdx_SS (the former local % MMP, 11/09/2026
% KRONSUB2IND_SS)                                                           % MMP, 11/09/2026
dimsRow = dimsOld(:).';                                                     % MMP, 11/09/2026
if isempty(dimsRow)                                                         % MMP, 11/09/2026
    oldLin = ones(size(idxOld,1),1);                                        % MMP, 11/09/2026
elseif numel(dimsRow) == 1                                                  % MMP, 11/09/2026
    oldLin = idxOld(:,1);                                                   % MMP, 11/09/2026
else                                                                        % MMP, 11/09/2026
    subs = num2cell(fliplr(idxOld),1);                                      % MMP, 11/09/2026
    oldLin = sub2ind(fliplr(dimsRow),subs{:});                              % MMP, 11/09/2026
    oldLin = oldLin(:);                                                     % MMP, 11/09/2026
end                                                                         % MMP, 11/09/2026

P = sparse(oldLin,(1:N).',1,N,N);
end

function idx = kronMultiIdx_SS(dims)
% IDX = KRONMULTIIDX_SS(DIMS) lists the multi-indices of a tensor monomial
% vector, first variable outermost and last fastest, one row per monomial.
% INPUT   dims: 1 x d number of monomials per variable.
% OUTPUT  idx: prod(dims) x d array.

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

