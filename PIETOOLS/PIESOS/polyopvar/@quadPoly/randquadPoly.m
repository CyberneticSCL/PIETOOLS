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

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PIETOOLS - randquadPoly
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
% Initial coding SS - 01/16/2026

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
%
% INPUT
% - dTot:   positive integer, the size of the tensor-product basis;
% - nvar:   nonnegative integer, the number of variables;
% - deg:    scalar or 1 x nvar array of nonnegative integers, the maximal
%           degree deg_i of each variable, raised to ni-1 where needed;
%
% OUTPUT
% - Zcell:  1 x nvar cell of sorted unique exponent column vectors, with
%           prod_i numel(Zcell{i}) = dTot ({} for nvar = 0).

if nvar == 0
    % no variables => constant basis of size 1
    Zcell = {};
    return;
end

% Choose per-variable counts ni with product exactly dTot
% Per-variable counts nPer with product exactly dTot: each prime factor of
% dTot goes to the currently smallest bucket (dTot = 12, nvar = 3 gives
% [2 2 3])
pf = factor(dTot);
nPer = ones(1, nvar);
for p = pf
    [~, idx] = min(nPer);
    nPer(idx) = nPer(idx) * p;
end
nPer = nPer(:).';

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

