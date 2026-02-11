function [S,Sd,Z3] = monomial_outer_multiple(Z1,Z2,n1,n2)
%monomial_outer_multiple This function perform Z1(n1)*Z2(n2)^T
% Inputs:
%   Z1: cell array representing tensor monomial product,
%       Z1{n1(1)}\otimes..\otimes Z1{n1(end)}
%   Z2: cell array representing tensor monomial product,
%       Z2{n2(1)}\otimes..\otimes Z2{n2(end)}
%   n1: variable names for Z1 monomials
%   n2: variable names for Z2 monomials
%
% Outputs: (I \otimes Z1)*(I\otimes Z2') 
%          = (I\otimes Z3')(I\otimes Sd) 
%          = (I\otimes S)(I\otimes Z3)
%


nu = union(n1,n2);
n_nu = numel(nu);
% init lifted cells with constant-only basis
Z1u = repmat({0}, 1, n_nu);
Z2u = repmat({0}, 1, n_nu);

% map names -> positions in nu, then assign
[tf1, loc1] = ismember(n1, nu);
Z1u(loc1(tf1)) = Z1(tf1);

[tf2, loc2] = ismember(n2, nu);
Z2u(loc2(tf2)) = Z2(tf2);

S_cells  = cell(1,n);
Sd_cells = cell(1,n);
Z3       = cell(1,n);

% n2(i) = length(Z2i);  m3(i) = length(Z3i)
n2 = zeros(1,n);
m3 = zeros(1,n);

% --- 1) Per-variable maps ---
for i = 1:n
    [Si, Sdi, Z3i] = monomial_outer(Z1u{i}, Z2u{i});
    S_cells{i}  = sparse(Si);
    Sd_cells{i} = sparse(Sdi);
    Z3{i}       = Z3i;

    n2(i) = numel(Z2{i});
    m3(i) = size(Z3i,1);
end

% --- 2) Build big Kronecker products efficiently (indices-only) ---
Sbig  = kron_sparse_chain(S_cells);
Sdbig = kron_sparse_chain(Sd_cells);

% --- 3) Apply permutations as row reindexing (no permutation matrices) ---
% P_L is the same regrouping permutation from Thm(leftPermute)
pL = leftPermuteVec(n2, m3);        % PL*x = x(pL)
pR = invertPerm(pL);                % PR*x = x(pR), since PR = PL'

% Left-multiply by PL / PR = row permutation:
S  = Sbig(pL, :);
Sd = Sdbig(pR, :);
end

% ======================================================================
% Efficient sparse Kronecker product of a chain of sparse matrices
% without calling kron repeatedly on large sparse intermediates.
% ======================================================================
function K = kron_sparse_chain(Acells)
% K = A1 ⊗ A2 ⊗ ... ⊗ Ak
% Uses find() triplets and combines indices iteratively.

if isempty(Acells)
    K = sparse(1,1,1);
    return;
end

% Start with first factor's triplets
A = sparse(Acells{1});
[i, j, v] = find(A);
[r, c] = size(A);

for t = 2:numel(Acells)
    B = sparse(Acells{t});
    [ib, jb, vb] = find(B);
    [rb, cb] = size(B);

    % Kronecker index rule:
    % (iA,jA) with value vA and (iB,jB) with value vB produce:
    % iK = (iA-1)*rb + iB
    % jK = (jA-1)*cb + jB
    % vK = vA*vB

    na = numel(v);
    nb = numel(vb);

    % Combine indices without forming dense grids
    i = repelem((i-1)*rb, nb) + repmat(ib, na, 1);
    j = repelem((j-1)*cb, nb) + repmat(jb, na, 1);

    % Values
    v = repelem(v, nb) .* repmat(vb, na, 1);

    % Update sizes
    r = r * rb;
    c = c * cb;
end

% Build final sparse (duplicates summed automatically)
K = sparse(i, j, v, r, c);
end

function invp = invertPerm(p)
invp = zeros(size(p));
invp(p) = 1:numel(p);
end