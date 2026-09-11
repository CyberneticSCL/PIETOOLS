function [varsC,ZC,CC] = leftShiftMonomials_SS(varsA,ZA,CA,varsB,ZB,CB)
% leftShiftMonomials_SS
%
% Computes
%
%   ((I_p \otimes ZA')*CA{I,J})*((I_q \otimes ZB')*CB{I})
%       = (I_p \otimes ZC')*CC{I,J}
%
% for each cell index I,J.
%
% INPUTS
%   varsA : variables for ZA
%   ZA    : 1-by-nA cell array of exponent vectors
%   CA    : nI-by-nJ cell array, CA{I,J} has size p*NA by q
%   varsB : variables for ZB
%   ZB    : 1-by-nB cell array of exponent vectors
%   CB    : 1-by-nI or nI-by-1 cell array, CB{I} has size q*NB by r
%
% OUTPUTS
%   varsC : merged variable list
%   ZC    : merged/product monomial basis
%   CC    : nI-by-nJ cell array, CC{I,J} has size p*NC by r
%
% MMP, 09/10/2026: Hoist the 'kron(speye(p),Qmat)' lift out of the I,J loops,
%                  where it was rebuilt nI*nJ times identically; check the
%                  cell shapes once instead of inside them; preallocate CC.
%                  Taking Sec. 5.1's partial sum over I inside this routine
%                  was tried and rejected as slower -- see the note at the
%                  coefficient loop.

% % % % % % % % % % COMPUTATION v3 % % % % % % % %
% We need to compute (I_p otimes Z1^T) A (I_q otimes Z^2T) B
% we have that A (I_q otimes Z2^T) = (A otimes 1)(I_q otimes Z2) =
%              A otimes Z2^T =(I_{p*Na} otimes Z2^T)(A otimes I_Nb)
%
%
% Also, we have
%      (I_p otimes Z1^T)(I_{p*Na} otimes Z2^T) =
%     (I_p otimes Z1^T)(I_{p} otimes I_{Na} otimes Z2^T) =
%      I_p otimes (Z1^T (I_{Na} otimes Z2^T)) = 
%      I_p otimes (Z1^T otimes 1)(I_{Na} otimes Z2^T) =
%      I_p otimes (Z1^T otimes Z2^T) 
% Given (Z1^T otimes Z2^T) = Z3^T Q, we have
%      I_p otimes (Z1^T otimes Z2^T) = I_p otimes (Z3^T Q) =
%      (I_p otimes Z3^T)(I_p otimes Q)
%
% Then, 
% Then, (I_p otimes Z1^T) A (I_q otimes Z^2T) B =
%       (I_p otimes Z3^T)(I_p otimes Q)(A otimes I_Nb) B
%          p x           Q is Nc x Na*Nb 

varsA = varsA(:).';
varsB = varsB(:).';
ZA = ZA(:).';
ZB = ZB(:).';

NA = prod(cellfun(@numel,ZA));
NB = prod(cellfun(@numel,ZB));

if isempty(ZA)
    NA = 1;
end
if isempty(ZB)
    NB = 1;
end

[nI,nJ] = size(CA);

if numel(CB) ~= nI
    error('leftShiftMonomials_SS: CB must have one cell for each row of CA.');
end

[mA,q] = size(CA{1,1});
[mB,r] = size(CB{1});

if mod(mA,NA) ~= 0
    error('leftShiftMonomials_SS: size(CA{1},1) must be divisible by NA.');
end

if mB ~= q*NB
    error('leftShiftMonomials_SS: size(CB{1},1) must equal size(CA{1},2)*NB.');
end

p = mA/NA;

% Merge variables. This uses MATLAB sorted union, matching current convention.
varsC = union(varsA,varsB);
nC = numel(varsC);

% Build output basis.
ZC = cell(1,nC);

for t = 1:nC
    ia = find(strcmp(varsA,varsC{t}),1);
    ib = find(strcmp(varsB,varsC{t}),1);

    if isempty(ia)
        ZC{t} = ZB{ib}(:);
    elseif isempty(ib)
        ZC{t} = ZA{ia}(:);
    else
        ea = ZA{ia}(:);
        eb = ZB{ib}(:);
        ZC{t} = unique(reshape(ea + eb.',[],1));
    end
end

NC = prod(cellfun(@numel,ZC));
if isempty(ZC)
    NC = 1;
end

% Expand ZA exponents into varsC coordinates.
EAc = zeros(NA,nC);

for t = 1:numel(varsA)
    k = find(strcmp(varsC,varsA{t}),1);

    lenA = cellfun(@numel,ZA);
    left  = prod(lenA(1:t-1));
    right = prod(lenA(t+1:end));

    EAc(:,k) = kron(ones(left,1),kron(ZA{t}(:),ones(right,1)));
end

% Expand ZB exponents into varsC coordinates.
EBc = zeros(NB,nC);

for t = 1:numel(varsB)
    k = find(strcmp(varsC,varsB{t}),1);

    lenB = cellfun(@numel,ZB);
    left  = prod(lenB(1:t-1));
    right = prod(lenB(t+1:end));

    EBc(:,k) = kron(ones(left,1),kron(ZB{t}(:),ones(right,1)));
end

% Expand ZC exponents into varsC coordinates.
EC = zeros(NC,nC);

for t = 1:nC
    lenC = cellfun(@numel,ZC);
    left  = prod(lenC(1:t-1));
    right = prod(lenC(t+1:end));

    EC(:,t) = kron(ones(left,1),kron(ZC{t}(:),ones(right,1)));
end

% % % % % Computing Q such that ZA otimes ZB = Q ZC % % % % % 


% An empty merged variable list leaves nothing to permute, so Qmat is the   % MMP, 08/30/2026
% 1-by-1 identity. Without this branch ismember returns an empty index and  % MMP, 08/30/2026
% 'sparse' builds the 1-by-1 ZERO instead, which zeroed every composition   % MMP, 08/30/2026
% involving R^n. Matches the @sdopvar copy, which already had the guard.    % MMP, 08/30/2026
if nC == 0                                                                  % MMP, 08/30/2026
    Qmat = sparse(1,1,1,1,1);                                               % MMP, 08/30/2026
else                                                                        % MMP, 08/30/2026
    EC3dt = reshape(EAc, 1,  [],size(EAc, 2)) + reshape(EBc, [], 1, size(EBc, 2));
    % Now EC2d is a 2darray with EC2d(i, :) = EAc(*, :) + EBc(**, :)
    EC2dt = reshape(EC3dt, [], size(EC3dt, 3));
    [~, Idx_array] = ismember(EC2dt, EC, 'rows');
    Qmat = sparse(Idx_array, 1:(NA*NB), ones(size(Idx_array)), NC, NA*NB);
end                                                                         % MMP, 08/30/2026
% Qmat Z3 = Z12
% Idx_array = reshape(Idx_array, NA, NB); 
% Idx_array = reshape(I)
% Idx_array(i, j) = k if G(:, :, k) = ... + A_i B_j + ... 

% x = rand(1);
% Z12 = x.^EC2d;
% Z3 =  x.^EC;

% Cell shapes, checked once rather than inside the coefficient loops.       % MMP, 09/10/2026
for I = 1:nI                                                                % MMP, 09/10/2026
    if ~isequal(size(CB{I}),[q*NB,r])                                       % MMP, 09/10/2026
        error('leftShiftMonomials_SS: inconsistent CB cell dimensions.');   % MMP, 09/10/2026
    end                                                                     % MMP, 09/10/2026
    for J = 1:nJ                                                            % MMP, 09/10/2026
        if ~isequal(size(CA{I,J}),[p*NA,q])                                 % MMP, 09/10/2026
            error('leftShiftMonomials_SS: inconsistent CA cell dimensions.');% MMP, 09/10/2026
        end                                                                 % MMP, 09/10/2026
    end                                                                     % MMP, 09/10/2026
end                                                                         % MMP, 09/10/2026

% % % % % % % % % % NEW COMPUTATION OF G % % % % % % % %
% Then, (I_p otimes Z1^T) A (I_q otimes Z^2T) B =
%       (I_p otimes Z3^T)(I_p otimes Q)(A otimes I_Nb) B
%          p x           Q is Nc x Na*Nb

% Q is fixed by the bases alone, so lift it over the matrix dimension once. % MMP, 09/10/2026
% It was previously rebuilt inside both loops, nI*nJ times identically.     % MMP, 09/10/2026
IpQ = kron(speye(p),Qmat);                                                  % MMP, 09/10/2026

% Sec. 5.1 of the sopvar document takes the partial sum                     % MMP, 09/10/2026
%     L_betab = sum_betaa C^A_{betaa,betab} K_betaa                         % MMP, 09/10/2026
% on the common basis BEFORE the monomials are shifted, which would let     % MMP, 09/10/2026
% IpQ be applied once per J instead of once per (I,J). That was tried and   % MMP, 09/10/2026
% measured SLOWER (0.90x overall, 0.73x worst, at nBetaA = nAlphaB = 9):    % MMP, 09/10/2026
% Qmat compresses the monomial axis, NC <= NA*NB, so applying IpQ first     % MMP, 09/10/2026
% shrinks the rows from p*NA*NB to p*NC and the caller's sum afterwards     % MMP, 09/10/2026
% adds the SMALLER matrices. Pre-summing instead accumulates the larger     % MMP, 09/10/2026
% pre-shift matrix nI times. Left as the document has it only in spirit:    % MMP, 09/10/2026
% the sum stays with the caller.                                           % MMP, 09/10/2026
CC = cell(nI,nJ);                                                           % MMP, 09/10/2026
for I = 1:nI

    Bcoef = CB{I};

    for J = 1:nJ

        Acoef = CA{I,J};

        AINb= kron(Acoef, speye(NB));
        H = IpQ*AINb*Bcoef;
        CC{I, J} = H;
        % max(abs(H - HH))
    end
end

end