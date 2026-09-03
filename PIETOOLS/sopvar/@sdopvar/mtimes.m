function C = mtimes(A,B)
% Multiplication for sdopvar with numeric matrices/scalars.
if isa(A,'sdopvar') 
    if isnumeric(B)
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
        % If B is matrix
        % left_size = A.dims(1)*prod(cellfun(@numel,A.ZL));
        % R = kron(B,eye(prod(cellfun(@numel,A.ZR))));
        % params = MatrixMultiply(A.params,speye(left_size),R);
        % A is sdopvar and B is numeric
        nL = prod(cellfun(@numel,A.ZL));
        nR = prod(cellfun(@numel,A.ZR));
        
        L = speye(A.dims(1)*nL);
        R = kron(B,speye(nR));
        
        params = transformParams(A.params,L,R);
        C = sdopvar(params,A.vars,A.Zd,A.ZL,A.ZR,A.dom,[A.dims(1),size(B,2)]);
    elseif isa(B,'sopvar')
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
        % If A is matrix
        % right_size = B.dims(2)*prod(cellfun(@numel,B.ZR));
        % L = kron(A,eye(prod(cellfun(@numel,B.ZL))));
        % params = MatrixMultiply(B.params,L,speye(right_size));
        % A is numeric and B is sdopvar
        nL = prod(cellfun(@numel,B.ZL));
        nR = prod(cellfun(@numel,B.ZR));
        
        L = kron(A,speye(nL));
        R = speye(B.dims(2)*nR);
        
        params = transformParams(B.params,L,R);
        C = sdopvar(params,B.vars,B.Zd,B.ZL,B.ZR,B.dom,[size(A,1),B.dims(2)]);
    elseif isa(A,'sopvar')
        Ct = mtimes(B',A');
        C = Ct';
    end
end
end
function C = mtimesSdopvarAndSopvar(A,B)
if A.dims(2) ~= B.dims(1)
    error('mtimes:innerDimensions', ...
        'The output dimension of B must equal the input dimension of A.');
end

if ~isequal(A.vars.in(:),B.vars.out(:))
    error('mtimes:intermediateVariables', ...
        'B output variables must match A input variables in the same order.');
end

if ~isequal(A.dom.in,B.dom.out)
    error('mtimes:intermediateDomains', ...
        'B output domains must match A input domains.');
end

% Derive S1, S2, and S3 partitions from public properties.
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

% Cparams is indexed by the common variables that pass through both
% operators.
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

% Split beta into [beta_a,beta_b].
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

% Reorganize the parameter cells according to the split indices.
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

[Zhat_s2a,G_s3a,Zhat_s3b,PL,PR] = int_2b( ...
    A.ZR,B.ZL,A.vars.in, ...
    vs2a,vs2b,vs3a,vs3b,dom2b);

alpha_b_int = alpha_b;
alpha_b_int(alpha_b == 2) = 3;
alpha_b_int(alpha_b == 3) = 2;

[Cs2a_betaa,KZhat_s2a] = ...
    int_monomial(Zhat_s2a,beta_a,dom2a);

[Cs3b_alphab,KZhat_s3b] = ...
    int_monomial(Zhat_s3b,alpha_b_int,dom3b);

[C_gam_alp_beta,barZ_3aL,barZ_3aR] = ...
    int_semisep(G_s3a,beta_b,alpha_a,dom3a,Csize);

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

CB_s2a = cell(nBetaA,1);

for ia = 1:nBetaA
    CB_s2a{ia} = kron( ...
        speye(qL),Cs2a_betaa{ia}.');
end

[varLtemp,ZLtemp,CLtemp_betaab] = ...
    leftShiftMonomials_SD( ...
        A.vars.out,A.ZL,CA_PL, ...
        vs2a,KZhat_s2a,CB_s2a);

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


% Each row corresponds to one pair (CA row, decision variable).

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

[Cvarsout,CZL,CLtemp_betab] = ...
    leftShiftMonomials_SD( ...
        varLtemp,ZLtemp,CLtemp, ...
        vs3a,barZ_3aL,CB_barL);
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

for igamma = 1:numel(Cparams.A)

    Aaccum = sparse(nCoefC,1);
    Baccum = sparse(nDecision,nCoefC);

    for j = 1:size(beta_b,1)

        Lparam = CLtemp_betab{j};

        if mod(numel(Lparam.A),nRowsC) ~= 0
            error('mtimes:leftAssemblyDimensions', ...
                'The left affine coefficient has inconsistent dimensions.');
        end

        nColsL = numel(Lparam.A)/nRowsC;

        if size(Lparam.B,1) ~= nDecision || ...
                size(Lparam.B,2) ~= numel(Lparam.A)
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

            % vec(L*MR) = (MR' kron I)*vec(L)
            T = kron(MR.',speye(nRowsC));

            Aaccum = Aaccum + T*Lparam.A;
            Baccum = Baccum + Lparam.B*T.';
        end
    end

    Cparams.A{igamma} = Aaccum;
    Cparams.B{igamma} = Baccum;
end
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

Cparams = expandCommonFullIntegrals_SD( ...
    Cparams,vs3a,CvarsS3);

Cvars = struct( ...
    'in',{Cvarsin}, ...
    'out',{Cvarsout});

Cdom = struct( ...
    'in',B.dom.in(idxIn,:), ...
    'out',A.dom.out(idxOut,:));

C = sdopvar( ...
    Cparams,Cvars,A.Zd,CZL,CZR,Cdom,Cdims);
end


% helper functions
function paramsOut = expandCommonFullIntegrals_SD( ...
    paramsIn,baseVars,finalVars)

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

        sourceIndex = parameterMultiIndexToLinear(gammaBase);
    end

    paramsOut.A{k} = paramsIn.A{sourceIndex};
    paramsOut.B{k} = paramsIn.B{sourceIndex};
end
end
function idx = parameterMultiIndex(n)

if n == 0
    idx = zeros(1,0);
else
    idx = fliplr( ...
        dec2base(0:3^n-1,3,n)-'0') + 1;
end
end

function linearIndex = parameterMultiIndexToLinear(idx)

linearIndex = 1;

for k = 1:numel(idx)
    linearIndex = linearIndex ...
                + (idx(k)-1)*3^(k-1);
end
end
function selectedDom = selectDomain(selectedVars,allVars,allDom,label)

[found,idx] = ismember(selectedVars,allVars);

if any(~found)
    error('mtimes:missingVariables', ...
        '%s variables were not found in the expected variable list.',label);
end

selectedDom = allDom(idx,:);
end
function grouped = regroupParameterCells( ...
    paramCells,multiIdx,firstCols,secondCols)

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

idx = 1;

for k = 1:numel(alpha)
    idx = idx + (alpha(k)-1)*3^(k-1);
end
end
function paramsOut = transformParams(paramsIn,L,R)

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