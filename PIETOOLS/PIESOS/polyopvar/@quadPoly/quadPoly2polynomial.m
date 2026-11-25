function P = quadPoly2polynomial(Q)
%QUADPOLY2POLYNOMIAL Convert quadPoly to polynomial object.
% Static method, called as quadPoly.quadPoly2polynomial(Q).
%
% INPUTS
% - Q:      m x n 'quadPoly' object (or 'double' matrix);
%
% OUTPUTS
% - P:      m x n 'polynomial' object in the variables Q.ns and Q.nt.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PIETOOLS - quadPoly2polynomial
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
% For support, contact M. Peet, Arizona State University at mpeet@asu.edu
% SS, 02/11/2026: Initial coding;

    [m,n] = size(Q);
    if isdouble(Q)
        P = polynomial(Q);
        return
    elseif isempty(Q.C) || nnz(Q.C) == 0
        coeff  = sparse(1, m*n);   % 1 term, coefficients all zero
        degmat = sparse(1, 0);     % 1 term, 0 variables
        varname = {};              % no variables
        P = polynomial(coeff, degmat, varname, [m n]);
        return;
    end

    % Build unified variable list: [ns, nt]
    varname = [Q.ns(:); Q.nt(:)];
    V = numel(varname);

    % Determine tensor basis sizes
    ds = basisSize(Q.Zs);
    dt = basisSize(Q.Zt);

    % Expand tensor exponents into full monomial degree table
    % Ds: ds×|ns|, Dt: dt×|nt|
    Ds = fullDegreeTable(Q.Zs);
    Dt = fullDegreeTable(Q.Zt);

    % Term list corresponds to all pairs (ks, kt): total T = ds*dt
    T = ds * dt;

    % Build degmat (T×V): [Ds(ks,:), Dt(kt,:)]
    degmat = sparse(T, V);
    % Fill left degrees
    if ~isempty(Ds)
        % repeat Ds for each kt
        leftBlock = repelem(Ds, dt, 1); % (ds*dt)×|ns|
        degmat(:, 1:size(leftBlock,2)) = sparse(leftBlock);
    end
    % Fill right degrees
    if ~isempty(Dt)
        % tile Dt for each ks
        rightBlock = repmat(Dt, ds, 1); % (ds*dt)×|nt|
        degmat(:, (size(Ds,2)+1):V) = sparse(rightBlock);
    end

    % Build coefficient matrix: Tx(m*n) where column is entry (i + m*(j-1))
    coeff = sparse(T, m*n);

    % Q.C is (m*ds)×(n*dt). Each block (i,j) is ds×dt.
    % For each matrix entry (i,j), pull its ds×dt block, vectorize in the same
    % term ordering (ks outer, kt inner) used above.
    for j = 1:n
        colBlock = (j-1)*dt + (1:dt);  % columns in Q.C per kt for this j
        for i = 1:m
            rowBlock = (i-1)*ds + (1:ds); % rows in Q.C per ks for this i
            Bij = Q.C(rowBlock, colBlock); % ds×dt
            v = Bij.';                    % (ds*dt)×1, kt varies fastest (MATLAB)
            v = v(:);
            idxCol = i + m*(j-1);
            coeff(:, idxCol) = sparse(v);
        end
    end

    % Optionally drop identically-zero terms to keep polynomial small
    nzTerms = any(coeff ~= 0, 2);
    coeff = coeff(nzTerms, :);
    degmat = degmat(nzTerms, :);

    % Construct polynomial
    P = polynomial(coeff, degmat, varname, [m n]);
end

% ---------------- helpers ----------------
function d = basisSize(Zcell)
% D = BASISSIZE(ZCELL) returns the number of monomials of the tensor
% product basis defined by the exponent vectors in ZCELL (an empty
% element counts as the constant monomial).
%
% INPUTS
% - Zcell:  1 x k cell of exponent vectors;
%
% OUTPUTS
% - d:      positive integer, the product of their lengths (1 if empty).
    if isempty(Zcell), d = 1; return; end
    d = 1;
    for k = 1:numel(Zcell)
        z = Zcell{k};
        if isempty(z), z = 0; end
        d = d * numel(z);
    end
end

function D = fullDegreeTable(Zcell)
% Return full tensor degree table: (prod_i |Z{i}|) x (#vars)
% kron order: last variable varies fastest
%
% INPUTS
% - Zcell:  1 x k cell of exponent vectors (an empty one counts as degree
%           0);
%
% OUTPUTS
% - D:      prod_i |Z{i}| x k array, row j the degrees of monomial j
%           (1 x 0 for an empty Zcell).
    if isempty(Zcell)
        D = zeros(1,0);
        return;
    end
    k = numel(Zcell);
    sizes = cellfun(@(z) numel(z(:)), Zcell);
    N = prod(sizes);

    % subs: the N x k tensor subscripts (1-based) of all positions, the last
    % variable varying fastest
    subs = zeros(N,k);
    idx = (0:N-1).';
    for ii = k:-1:1
        si = sizes(ii);
        subs(:,ii) = mod(idx, si) + 1;
        idx = floor(idx / si);
    end
    D = zeros(N, k);
    for i = 1:k
        zi = Zcell{i}(:);
        if isempty(zi), zi = 0; end
        D(:,i) = zi(subs(:,i));
    end
end

