function Q = Sym2quadPoly(Fsym, sVars, tVars)
%Sym2quadPoly Convert symbolic polynomial/matrix to quadPoly object.
%
%   Q = quadPoly_fromSym(Fsym, sVars, tVars)
%
% Inputs:
%   Fsym  : sym scalar or sym matrix (m x n)
%   sVars : sym vector [s1 ... sk] (may be empty)
%   tVars : sym vector [t1 ... tl] (may be empty)
%
% Output:
%   Q : quadPoly with fields C, Zs, Zt, dim, ns, nt
%
% Representation:
%   F(s,t) = (I_m ⊗ Z(s)^T) * C * (I_n ⊗ Z(t))
% where Z(s) = Zs1(s1) ⊗ ... ⊗ Zsk(sk), each Zsi is monomials of degrees in Zs{i}.
%
% Notes:
%   - Coefficients must be convertible to double. If not, error.
%   - Degree sets Zs{i}, Zt{j} are inferred from monomials present in Fsym.
%   - Ordering matches MATLAB kron ordering: last variable varies fastest.

    if ~isa(Fsym, 'sym')
        error('Fsym must be a symbolic (sym) scalar or matrix.');
    end

    if nargin < 2 || isempty(sVars), sVars = sym([]); end
    if nargin < 3 || isempty(tVars), tVars = sym([]); end

    sVars = sVars(:).'; % row
    tVars = tVars(:).'; % row

    [m, n] = size(Fsym);
    k = numel(sVars);
    l = numel(tVars);

    % Variable names
    ns = arrayfun(@char, sVars, 'UniformOutput', false);
    nt = arrayfun(@char, tVars, 'UniformOutput', false);

    % Collect per-variable degree sets appearing anywhere in Fsym
    degS = cell(1, k);
    degT = cell(1, l);
    for i = 1:k, degS{i} = []; end
    for j = 1:l, degT{j} = []; end

    % Helper to get exponent of monomial wrt a variable
    expOf = @(mon, var) double(feval(symengine, 'degree', mon, var));

    Fsym = expand(Fsym);

    for ii = 1:m
        for jj = 1:n
            p = Fsym(ii,jj);
            if isequal(p, sym(0))
                continue;
            end

            varsAll = [sVars, tVars];
            if isempty(varsAll)
                % constant only; no degrees to collect
                continue;
            end

            [~, mons] = coeffs(p, varsAll); % monomials present
            if isempty(mons)
                continue;
            end

            for r = 1:numel(mons)
                mon = mons(r);
                for a = 1:k
                    degS{a}(end+1) = expOf(mon, sVars(a)); %#ok<AGROW>
                end
                for b = 1:l
                    degT{b}(end+1) = expOf(mon, tVars(b)); %#ok<AGROW>
                end
            end
        end
    end

    % If variable exists but never appears, include degree 0 so basis is valid.
    Zs = cell(1,k);
    Zt = cell(1,l);
    for a = 1:k
        if isempty(degS{a})
            Zs{a} = 0;
        else
            Zs{a} = unique(degS{a}(:));
        end
        Zs{a} = sort(Zs{a});
        Zs{a} = Zs{a}(:); % column
    end
    for b = 1:l
        if isempty(degT{b})
            Zt{b} = 0;
        else
            Zt{b} = unique(degT{b}(:));
        end
        Zt{b} = sort(Zt{b});
        Zt{b} = Zt{b}(:); % column
    end

    % Basis sizes
    nS = ones(1,k);
    nT = ones(1,l);
    for a = 1:k, nS(a) = numel(Zs{a}); end
    for b = 1:l, nT(b) = numel(Zt{b}); end
    dS = prod(nS);
    dT = prod(nT);

    % Degree->position maps per variable
    mapS = cell(1,k);
    mapT = cell(1,l);
    for a = 1:k
        mapS{a} = containers.Map('KeyType','double','ValueType','double');
        for p = 1:numel(Zs{a})
            mapS{a}(double(Zs{a}(p))) = p;
        end
    end
    for b = 1:l
        mapT{b} = containers.Map('KeyType','double','ValueType','double');
        for p = 1:numel(Zt{b})
            mapT{b}(double(Zt{b}(p))) = p;
        end
    end

    % Multi-index -> linear index consistent with kron(Z1, Z2, ..., Zk)
    % last variable varies fastest
    linIndex = @(idxVec, nVec) 1 + sum((idxVec(:)'-1).*[cumprod(nVec(2:end)), 1]);

    % Preallocate sparse triplets for C
    % C size: (m*dS) x (n*dT)
    I = [];
    J = [];
    V = [];

    varsAll = [sVars, tVars];

    for ii = 1:m
        for jj = 1:n
            p = Fsym(ii,jj);
            if isequal(p, sym(0))
                continue;
            end

            if isempty(varsAll)
                % pure constant
                cst = double(p);
                if ~isfinite(cst), error('Non-finite constant.'); end
                row0 = (ii-1)*dS + 1;
                col0 = (jj-1)*dT + 1;
                I(end+1,1) = row0; %#ok<AGROW>
                J(end+1,1) = col0; %#ok<AGROW>
                V(end+1,1) = cst;  %#ok<AGROW>
                continue;
            end

            [cs, mons] = coeffs(p, varsAll);
            if isempty(mons)
                continue;
            end

            for r = 1:numel(mons)
                mon = mons(r);
                cr  = cs(r);

                % Extract exponent positions for s and t
                idxS = ones(1,k);
                idxT = ones(1,l);

                for a = 1:k
                    ea = expOf(mon, sVars(a));
                    if ~isKey(mapS{a}, ea)
                        error('Internal: missing degree in Zs for %s.', ns{a});
                    end
                    idxS(a) = mapS{a}(ea);
                end
                for b = 1:l
                    eb = expOf(mon, tVars(b));
                    if ~isKey(mapT{b}, eb)
                        error('Internal: missing degree in Zt for %s.', nt{b});
                    end
                    idxT(b) = mapT{b}(eb);
                end

                aLin = 1; bLin = 1;
                if k > 0, aLin = linIndex(idxS, nS); end
                if l > 0, bLin = linIndex(idxT, nT); end

                % Numeric coefficient
                cnum = double(cr);
                if ~isfinite(cnum)
                    error('Coefficient not convertible to finite double (entry %d,%d).', ii, jj);
                end

                row = (ii-1)*dS + aLin;
                col = (jj-1)*dT + bLin;

                I(end+1,1) = row; %#ok<AGROW>
                J(end+1,1) = col; %#ok<AGROW>
                V(end+1,1) = cnum; %#ok<AGROW>
            end
        end
    end

    C = sparse(I, J, V, m*dS, n*dT);

    % Build quadPoly object
    Q = quadPoly(C,Zs,Zt,[m,n],ns,nt);
end
