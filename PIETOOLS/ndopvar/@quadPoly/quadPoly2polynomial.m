function P = quadPoly2polynomial(Q)
%quadPoly2polynomial Convert quadPoly to polynomial object.

    m = Q.dim(1);
    n = Q.dim(2);
    
    if isempty(Q.C) || nnz(Q.C) == 0
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
            v = Bij(:);                    % (ds*dt)×1, kt varies fastest (MATLAB)
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
    if isempty(Zcell), d = 1; return; end
    d = 1;
    for k = 1:numel(Zcell)
        z = Zcell{k};
        if isempty(z), z = 0; end
        d = d * numel(z);
    end
end

function D = fullDegreeTable(Zcell)
% Return full tensor degree table: (prod_i |Z{i}|) × (#vars)
% kron order: last variable varies fastest (consistent with quadPoly_eval)
    if isempty(Zcell)
        D = zeros(1,0);
        return;
    end
    k = numel(Zcell);
    sizes = cellfun(@(z) numel(z(:)), Zcell);
    N = prod(sizes);

    subs = tensorDecompAll(N, sizes); % N×k positions
    D = zeros(N, k);
    for i = 1:k
        zi = Zcell{i}(:);
        if isempty(zi), zi = 0; end
        D(:,i) = zi(subs(:,i));
    end
end

function subs = tensorDecompAll(N, sizes)
% N = prod(sizes). Return N×k subscripts (1-based), last var varies fastest.
    k = numel(sizes);
    subs = zeros(N,k);
    idx = (0:N-1).';
    for i = k:-1:1
        si = sizes(i);
        subs(:,i) = mod(idx, si) + 1;
        idx = floor(idx / si);
    end
end
