function [varsC,ZC,CC] = leftShiftMonomials_SS(varsA,ZA,CA,varsB,ZB,CB)
% leftShiftMonomials
%   ((I\otimes ZA')*CA)*((I \otimes ZB')*CB) = (I\otimes ZC')*CC
%
% INPUTS
%   varsA : cell array of variable names for ZA
%   ZA    : cell array of exponent vectors for ZA
%   CA    : mxn cell of (p*NA) x q coefficient matrices
%   varsB : cell array of variable names for ZB
%   ZB    : cell array of exponent vectors for ZB
%   CB    : 1xm cell of (q*NB) x r coefficient matrices
%
% OUTPUTS
%   varsC : output variable list
%   ZC    : output tensor-product monomial basis
%   CC    : (p*NC) x r coefficient matrix
NA = prod(cellfun(@numel, ZA));  % length of ZA monomial vec
NB = prod(cellfun(@numel, ZB));  % length of ZB



% if varsA or varsB is empty, we include a dummy variable
if NA == 0
    varsA = {'<dummy_varA>'};
    ZA = {[0]};
    NA = 1;
end

if NB == 0
    varsB = {'<dummy_varB>'};
    ZB = {[0]};
    NB = 1;
end

[mA,q] = size(CA{1}); %must be length(NA)*p x q   if (I\otimes ZA')*CA is pxq
[~,r] = size(CB{1}); %as above

p  = mA/NA;  % must be integer by default

varsC = union(varsA, varsB);  % just merge the vars
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
        ZC{t} = unique(reshape(va + vb.', [], 1));
    end
end

NC = prod(cellfun('length', ZC));

% Expand ZA into exponent table aligned to varsC
% ith row of EAc is the ith monomial in ZA but padded with varsC exponents 
% i.e., varsA^{ZA} = varsC^EAc
EAc = zeros(NA, nC);
if ~isempty(varsA)
    lenA = cellfun(@numel, ZA); % number of exponents in each varsA
    for t = 1:numel(varsA)
        k = find(strcmp(varsC, varsA{t}), 1); % find location of varsA in varsC
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
for I=1:idxA
    for J=1:idxB
% A(:,:,i) = coefficient matrix of i-th monomial in A
% B(:,:,j) = coefficient matrix of j-th monomial in B
% in other words, 
% (I\otimes ZA')*CA = sum_i A(:,:,i)*varsC{i}^EAc(i,:)
% (I\otimes ZB')*CA = sum_i B(:,:,i)*varsC{i}^EBc(i,:)
A = permute(reshape(CA{I,J}, [NA, p, q]), [2 3 1]);
B = permute(reshape(CB{I}, [NB, q, r]), [2 3 1]);

% Accumulate output coefficients
G = zeros(p, r, NC); % G(:,:,k) = coefficient matrix of k-th monomial in C
% (I\otimes ZC')*CC = sum_i G(:,:,i)*varsC{i}^EC(i,:)
for i = 1:NA  % now multiply each monomial of left polynomial with monomial of right polynomial
    Ai = A(:,:,i);
    ei = EAc(i,:);
    for j = 1:NB
        e = ei+EBc(j,:); % we can add since all exponents are aligned to varsC
        k = find(all(EC==e,2),1); % find where in ZC these should go
        G(:,:,k) = G(:,:,k)+Ai*B(:,:,j);
    end
end

% Pack back to (p*NC) x r
CC{I,J} = reshape(permute(G, [3 1 2]), [NC*p, r]);
    end
end

% now we remove dummy variable if they exist
var_names_removed = cellfun( ...
    @(name) ~isempty(char(eraseBetween(string(name), "<", ">",'Boundaries','inclusive'))),...
    varsC, 'UniformOutput',false); % names to rename
var_names_removed = cell2mat(var_names_removed)
varsC = varsC(var_names_removed);
ZC = ZC(var_names_removed);

end