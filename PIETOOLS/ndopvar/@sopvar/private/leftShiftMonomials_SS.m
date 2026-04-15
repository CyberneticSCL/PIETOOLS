function [varsC,ZC,CC] = leftShiftMonomials_SS(varsA,ZA,CA,varsB,ZB,CB)
% leftShiftMonomials
%   ((I\otimes ZA')*CA)*((I \otimes ZB')*CB) = (I\otimes ZC')*CC
%
% INPUTS
%   varsA : cell array of variable names for ZA
%   ZA    : cell array of exponent vectors for ZA
%   CA    : (p*NA) x q coefficient matrix
%   varsB : cell array of variable names for ZB
%   ZB    : cell array of exponent vectors for ZB
%   CB    : (q*NB) x r coefficient matrix
%
% OUTPUTS
%   varsC : output variable list
%   ZC    : output tensor-product monomial basis
%   CC    : (p*NC) x r coefficient matrix
NA = prod(cellfun(@numel, ZA));  % length of ZA monomial vec
NB = prod(cellfun(@numel, ZB));  % length of ZB

[mA,q] = size(CA); %must be length(NA)*p x q   if (I\otimes ZA')*CA is pxq
[mB,r] = size(CB); %as above

p  = mA/NA;  % must be integer by default
qB = mB/NB;

varsC = union(varsA, varsB, 'stable');  % just merge the vars
nC = numel(varsC);

% Build output basis per variable
ZC = cell(1,nC);
for t = 1:nC
    ia = find(strcmp(varsA, varsC{t}), 1);
    ib = find(strcmp(varsB, varsC{t}), 1);

    if isempty(ia) % if not present in ZA, just pick ZB
        ZC{t} = ZB{ib}(:);
    elseif isempty(ib) % if not present in ZB, just pick ZA
        ZC{t} = ZA{ia}(:);
    else  % if present in both, multiply as kronecker product
        va = ZA{ia}(:);
        vb = ZB{ib}(:);
        ZC{t} = unique(reshape(va + vb.', [], 1), 'stable');
    end
end

NC = prod(cellfun('length', ZC));

% Expand ZA into exponent table aligned to varsC
EAc = zeros(NA, nC);
if ~isempty(varsA)
    lenA = cellfun(@numel, ZA);
    for t = 1:numel(varsA)
        k = find(strcmp(varsC, varsA{t}), 1);
        left  = prod(lenA(1:t-1));
        right = prod(lenA(t+1:end));
        EAc(:,k) = kron(ones(left,1), kron(ZA{t}(:), ones(right,1)));
    end
end

% Expand ZB into exponent table aligned to varsC
EBc = zeros(NB, nC);
if ~isempty(varsB)
    lenB = cellfun(@numel, ZB);
    for t = 1:numel(varsB)
        k = find(strcmp(varsC, varsB{t}), 1);
        left  = prod(lenB(1:t-1));
        right = prod(lenB(t+1:end));
        EBc(:,k) = kron(ones(left,1), kron(ZB{t}(:), ones(right,1)));
    end
end

% Expand ZC into exponent table
EC = zeros(NC, nC);
if ~isempty(varsC)
    lenC = cellfun(@numel, ZC);
    for t = 1:nC
        left  = prod(lenC(1:t-1));
        right = prod(lenC(t+1:end));
        EC(:,t) = kron(ones(left,1), kron(ZC{t}(:), ones(right,1)));
    end
end

% Reshape coefficients
% CA: (p*NA) x q -> A: p x q x NA
% CB: (q*NB) x r -> B: q x r x NB

idxA = size(CA,1); 
idxB = size(CA,2);
for i=1:idxA
    for j=1:idxB
A = permute(reshape(CA{i}{j}, [NA, p, q]), [2 3 1]);
B = permute(reshape(CB{i}, [NB, q, r]), [2 3 1]);

% Accumulate output coefficients
G = zeros(p, r, NC);
for i = 1:NA
    Ai = A(:,:,i);
    ei = EAc(i,:);
    for j = 1:NB
        e = ei+EBc(j,:);
        k = find(all(EC==e,2),1);
        G(:,:,k) = G(:,:,k)+Ai*B(:,:,j);
    end
end

% Pack back to (p*NC) x r
CC{i}{j} = reshape(permute(G, [3 1 2]), [NC*p, r]);
    end
end
end