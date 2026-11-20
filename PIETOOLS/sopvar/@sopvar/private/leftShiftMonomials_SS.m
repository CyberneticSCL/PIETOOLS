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
%   varsA : 1-by-nA cellstr, variables for ZA
%   ZA    : 1-by-nA cell array of exponent vectors
%   CA    : nI-by-nJ cell array, CA{I,J} has size p*NA by q
%   varsB : 1-by-nB cellstr, variables for ZB
%   ZB    : 1-by-nB cell array of exponent vectors
%   CB    : 1-by-nI or nI-by-1 cell array, CB{I} has size q*NB by r
%
% OUTPUTS
%   varsC : merged variable list, the sorted union of varsA and varsB
%   ZC    : merged/product monomial basis, per variable the sums of the
%           exponents of ZA and ZB
%   CC    : nI-by-nJ cell array, CC{I,J} has size p*NC by r
%
% For support, contact M. Peet, Arizona State University at mpeet@asu.edu

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PIETOOLS - leftShiftMonomials_SS
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
% Initial coding SS - 04/15/2026

% % % % % % % % % % DERIVATION % % % % % % % %
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

% Variables and bases as rows; basis sizes NA, NB (1 without variables)
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

% Cell and block sizes: one CB cell per row of CA, p = rows of A per
% monomial block
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

    % A variable of one side keeps its exponents; a shared variable gets
    % every sum of an exponent of ZA and one of ZB
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


% An empty merged variable list leaves nothing to permute, so Qmat is the
% 1-by-1 identity (ismember would return an empty index, and 'sparse' the
% 1-by-1 zero).
if nC == 0
    Qmat = sparse(1,1,1,1,1);
else
    EC3dt = reshape(EAc, 1,  [],size(EAc, 2)) + reshape(EBc, [], 1, size(EBc, 2));
    % Now EC2d is a 2darray with EC2d(i, :) = EAc(*, :) + EBc(**, :)
    EC2dt = reshape(EC3dt, [], size(EC3dt, 3));
    [~, Idx_array] = ismember(EC2dt, EC, 'rows');
    Qmat = sparse(Idx_array, 1:(NA*NB), ones(size(Idx_array)), NC, NA*NB);
end

% Cell shapes, checked once rather than inside the coefficient loops.
for I = 1:nI
    if ~isequal(size(CB{I}),[q*NB,r])
        error('leftShiftMonomials_SS: inconsistent CB cell dimensions.');
    end
    for J = 1:nJ
        if ~isequal(size(CA{I,J}),[p*NA,q])
            error('leftShiftMonomials_SS: inconsistent CA cell dimensions.');
        end
    end
end

% The products: (I_p otimes Z1^T) A (I_q otimes Z2^T) B =
%               (I_p otimes Z3^T)(I_p otimes Q)(A otimes I_Nb) B.
% Q is fixed by the bases alone, so it is lifted over the matrix dimension
% once, not inside the loops.
IpQ = kron(speye(p),Qmat);

% One product per cell (I,J), the sum over I not taken here: Qmat shrinks
% the rows from p*NA*NB to p*NC, so summing after the shift adds the
% smaller matrices (measured faster).
CC = cell(nI,nJ);
for I = 1:nI

    Bcoef = CB{I};

    for J = 1:nJ

        Acoef = CA{I,J};

        AINb= kron(Acoef, speye(NB));
        H = IpQ*AINb*Bcoef;
        CC{I, J} = H;
    end
end

end