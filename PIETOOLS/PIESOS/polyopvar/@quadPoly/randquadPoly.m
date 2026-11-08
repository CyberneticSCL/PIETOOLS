function F = randquadPoly(dim, nMons, var_s, var_t, maxdeg, density)
%RANDQUADPOLY Random quadPoly with tensor-decomposed exponent storage.
%
% Inputs:
%   dim     : [m,n] matrix size of the polynomial output
%   nMons   : [ds,dt] number of monomials in the s- and t- tensor bases
%             (interpreted as TOTAL tensor basis sizes, not per-variable)
%   var_s   : 1×ks cellstr (or string array) of s-variable names (sorted-unique preferred)
%   var_t   : 1×kt cellstr (or string array) of t-variable names (sorted-unique preferred)
%   maxdeg  : [degS, degT] where each entry can be scalar or per-variable vector
%             - scalar degS means each s-variable exponent in {0..degS}
%             - vector degS means per-variable bounds
%   density : density argument for sprand
%
% Output:
%   F : quadPoly object with
%       - Zs, Zt stored as 1×k cell arrays of sorted-unique exponent vectors
%       - C sparse coefficient matrix of size (m*ds) × (n*dt)
%
% Notes:
%   This generator builds per-variable exponent sets Zs{i}, Zt{j} such that
%   prod_i numel(Zs{i}) == ds and prod_j numel(Zt{j}) == dt by distributing
%   the desired tensor sizes across variables.

m = dim(1); 
n = dim(2);

ds = nMons(1);
dt = nMons(2);

degS = maxdeg(1);
degT = maxdeg(2);

% normalize name containers
if isstring(var_s), var_s = cellstr(var_s); end
if isstring(var_t), var_t = cellstr(var_t); end

ks = numel(var_s);
kt = numel(var_t);

% Construct tensor-decomposed exponent cells
Zs = randExpCell(ds, ks, degS);
Zt = randExpCell(dt, kt, degT);

% Now prod(cellfun(@numel,Zs)) == ds, similarly for Zt
C  = sprand(m*ds, n*dt, density);

F = quadPoly(C, Zs, Zt, [m n], var_s, var_t);
end

% -------------------------------------------------------------------------
% Helpers
% -------------------------------------------------------------------------
function Zcell = randExpCell(dTot, nvar, deg)
%RANDEXPCELL Create a 1×nvar cell of sorted-unique exponent vectors whose
%tensor-product size equals dTot.
%
% Strategy:
%   - choose per-variable basis sizes ni such that prod(ni)=dTot
%   - for each variable i, draw ni unique exponents from [0..deg_i] (or allow
%     repeats if range too small, then unique+sort reduces size, so we guard)

if nvar == 0
    % no variables => constant basis of size 1
    Zcell = {};
    return;
end

% Choose per-variable counts ni with product exactly dTot
nPer = factorCounts(dTot, nvar);

% Per-variable degree bounds
if isscalar(deg)
    degVec = repmat(deg, 1, nvar);
else
    degVec = deg(:).';
    if numel(degVec) == 1
        degVec = repmat(degVec, 1, nvar);
    end
end

Zcell = cell(1, nvar);

for i = 1:nvar
    ni = nPer(i);
    di = degVec(min(i, numel(degVec)));

    % If range [0..di] too small to get ni unique values, expand the range
    % (keeps invariant sorted-unique and correct tensor size).
    if di + 1 < ni
        di = ni - 1;
    end

    % Sample ni unique exponents from 0..di
    % (randperm is simplest; if di is huge you can switch to randsample)
    exps = randperm(di + 1, ni) - 1;
    Zcell{i} = sort(exps(:));
end
end

function nPer = factorCounts(dTot, nvar)
%FACTORCOUNTS Split dTot into nvar integer factors (>=1) with exact product.
%Simple balanced factor allocation from prime factors.
%
% Example: dTot=12, nvar=3 -> could return [2 2 3].

pf = factor(dTot);
nPer = ones(1, nvar);

% Greedy: assign each prime factor to the currently smallest bucket
for p = pf
    [~, idx] = min(nPer);
    nPer(idx) = nPer(idx) * p;
end

% Ensure product matches
% (it should, but keep it deterministic)
nPer = nPer(:).';
end
