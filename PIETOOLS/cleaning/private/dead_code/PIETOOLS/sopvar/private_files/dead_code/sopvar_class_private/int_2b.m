function [Z2aout, G3aout, Z3bout] = int_2b(ZL, ZR, ZLvar, ZRvar, s2a, s2b, s3a, s3b, lims)
% This performs the factorization
%   int(ZL*ZR',s2b) = (Im \otimes Z2a') * G3a * (In \otimes Z3b)
%
% where G3a is returned as a struct with fields
%   G3a.C : coefficient matrix C3a
%   G3a.Z : monomial basis Z3anew
%
% so that
%   G3a = (I3a \otimes Z3anew') * C3a

% separate ZL and ZR into groups 2a, 2b, 3a, 3b
ZL2a = pick_group(ZL, ZLvar, s2a);
ZR2a = pick_group(ZR, ZRvar, s2a);

ZL2b = pick_group(ZL, ZLvar, s2b);
ZR2b = pick_group(ZR, ZRvar, s2b);

ZL3a = pick_group(ZL, ZLvar, s3a);
ZR3a = pick_group(ZR, ZRvar, s3a);

ZL3b = pick_group(ZL, ZLvar, s3b);
ZR3b = pick_group(ZR, ZRvar, s3b);

% form pairwise exponent matrices for each variable group
Z2a = pairwise_sum_cells(ZL2a, ZR2a);
Z2b = pairwise_sum_cells(ZL2b, ZR2b);
Z3a = pairwise_sum_cells(ZL3a, ZR3a);
Z3b = pairwise_sum_cells(ZL3b, ZR3b);

% integrate the 2b block
% int_a^b s^E ds = (b^(E+1) - a^(E+1)) / (E+1)
a = lims(:,1);
b = lims(:,2);

C2b = 1;
for i = 1:length(s2b)
    E = Z2b{i};
    Ci = (b(i).^(E+1) - a(i).^(E+1)) ./ (E+1);
    C2b = kron(C2b, Ci);
end

% condense 2a, 3a, 3b blocks
% Z2a matrix  -> (I2a  \otimes Z2anew') * K2a
% Z3a matrix  -> (I3a0 \otimes Z3anew') * K3a
% Z3b matrix  -> K3b * (I3b \otimes Z3bnew)

[Z2anew, K2a] = condense_left_factor(Z2a);
[Z3anew, K3a] = condense_left_factor(Z3a);
[Z3bnew, K3b] = condense_right_factor(Z3b);

% build center in the form (I3a \otimes Z3anew') * C3a
% We start from
%   K2a \otimes C2b \otimes ((I3a0 \otimes Z3anew')*K3a) \otimes K3b
%
% and rewrite it as
%   (I3a \otimes Z3anew') * C3a

A = kron(K2a, C2b);
B = K3b;

r3a  = prod(cellfun(@length, ZL3a));
nz3a = prod(cellfun(@length, Z3anew));

mA = size(A,1);
nA = size(A,2);
mB = size(B,1);
nB = size(B,2);

% perfect shuffle to move Z3anew' to the far left of the B block
% maps [q, p] ordering to [p, q]
q = nz3a;
p = nB;
perm3a = reshape(1:q*p, [q, p]);
perm3a = reshape(perm3a.', 1, []);
P3a = sparse(1:q*p, perm3a, 1, q*p, q*p);
C3a = kron(kron(kron(A, speye(r3a)), speye(mB)), speye(nz3a)) ...
    * kron(speye(nA * r3a), P3a) ...
    * kron(kron(speye(nA), K3a), speye(nB));

% absorb the 2a outer permutation into C3a
% We have
%   ((I2a \otimes Z2anew') \otimes Irest) = (Im \otimes Z2anew') * P2a
%
% and want to absorb the induced permutation into C3a:
%   P2a * (I3a \otimes Z3anew') = (I3anew \otimes Z3anew') * Q2a
% so that
%   C3a <- Q2a * C3a

n2a  = prod(cellfun(@length, Z2anew));
n3a  = prod(cellfun(@length, Z3anew));
nI2a = prod(cellfun(@length, ZL2a));
nI3a = size(C3a,1) / n3a;
nrest = nI3a;

% permutation for the outer-left factor
permBlock = reshape(1:n2a*nrest, [n2a, nrest]);
permBlock = reshape(permBlock.', [], 1);

% induced permutation on the row identity block of C3a
permOuter = reshape(1:nI2a*n2a*nrest, [n2a*nrest, nI2a]);
permOuter = permOuter(permBlock, :);
permOuter = permOuter(:);

% lift across untouched Z3anew rows
permCenter = reshape(1:nI3a*n3a, [nI3a, n3a]);
permCenter = permCenter(permOuter, :);
permCenter = permCenter(:);

C3a = C3a(permCenter, :);

% outputs
Z2aout = Z2anew;
Z3bout = Z3bnew;
G3aout = struct('C', C3a, 'Z', {Z3anew});

end


function Zg = pick_group(Z, Zvar, vars)
% Extract exponent vectors for the requested variable names.
% Missing variables are represented by scalar 0.

Zg = cell(1, length(vars));
for i = 1:length(vars)
    [tf, loc] = ismember(vars(i), Zvar);
    if tf
        Zg{i} = Z{loc};
    else
        Zg{i} = 0;
    end
end
end


function Zsum = pairwise_sum_cells(ZLgrp, ZRgrp)
% For each variable, form pairwise exponent sums ZL + ZR'.

Zsum = cell(1, length(ZLgrp));
for i = 1:length(ZLgrp)
    Zsum{i} = ZLgrp{i} + ZRgrp{i}';
end
end


function [Znew, K] = condense_left_factor(Z)
% Condense exponent matrix E into (I \otimes Znew') * K

Znew = cell(1, length(Z));
K = 1;

for i = 1:length(Z)
    E = Z{i};
    [u, ~, ic] = unique(E(:), 'stable');
    Znew{i} = u;

    Ki = sparse(ic, 1:numel(E), 1, length(u), numel(E));
    K = kron(K, Ki);
end
end


function [Znew, K] = condense_right_factor(Z)
% Condense exponent matrix E into K * (I \otimes Znew)

Znew = cell(1, length(Z));
K = 1;

for i = 1:length(Z)
    E = Z{i};
    [u, ~, ic] = unique(E(:), 'stable');
    Znew{i} = u;

    Ki = sparse(1:numel(E), ic, 1, numel(E), length(u));
    K = kron(K, Ki);
end
end