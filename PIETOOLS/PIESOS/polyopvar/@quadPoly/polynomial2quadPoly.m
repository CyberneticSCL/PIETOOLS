function Q = polynomial2quadPoly(P, sNames, tNames)
%POLYNOMIAL2QUADPOLY Convert polynomial object to quadPoly.
%
% Inputs:
%   P      : polynomial object
%   sNames : cellstr of left variables (ns)
%   tNames : cellstr of right variables (nt)
%
% Output:
%   Q      : quadPoly with tensor-product bases per variable.
%
% Notes:
%   - Builds Zs/Zt as unique degrees per variable seen in P.degmat.
%   - Requires that P.varname contains only sNames∪tNames (constants allowed).

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PIETOOLS - polynomial2quadPoly
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
% Initial coding SS - 02/11/2026

if isa(P,'double')
    P = polynomial(P);
elseif ~isa(P,'polynomial')
    error("Input must be of type 'polynomial'.")
end

    m = P.matdim(1);
    n = P.matdim(2);

    var = P.varname(:);
    deg = full(P.degmat);          % T×V (can be big but manageable)
    coeff = sparse(P.coefficient); % T×(m*n)

    % Map polynomial variables into left/right positions
    [tfS, posS] = ismember(sNames(:), var);
    [tfT, posT] = ismember(tNames(:), var);

    % Build Zs/Zt degree lists (unique degrees seen per variable)
    Zs = cell(1, numel(sNames));
    for i = 1:numel(sNames)
        if posS(i) == 0
            Zs{i} = 0; % missing variable -> only degree 0
        else
            Zs{i} = unique(deg(:, posS(i)), 'sorted');
            if isempty(Zs{i}), Zs{i} = 0; end
        end
    end

    Zt = cell(1, numel(tNames));
    for i = 1:numel(tNames)
        if posT(i) == 0
            Zt{i} = 0; % missing variable -> only degree 0
        else
            Zt{i} = unique(deg(:, posT(i)), 'sorted');
            if isempty(Zt{i}), Zt{i} = 0; end
        end
    end

    ds = basisSize(Zs);
    dt = basisSize(Zt);

    % Precompute maps from degree tuples -> tensor index for s and t
    % We use ismember per variable and then linearize subscripts.
    sSizes = cellfun(@numel, Zs);
    tSizes = cellfun(@numel, Zt);
    sStride = kronStrides(sSizes);
    tStride = kronStrides(tSizes);

    % For each term row in P, compute ks (1..ds) and kt (1..dt)
    T = size(deg,1);
    ks = ones(T,1);
    kt = ones(T,1);

    for i = 1:numel(sNames)
        if posS(i) == 0
            % missing variable -> degree is implicitly 0 -> loc is always 1 -> add nothing
            continue;
        end
        [ok, loc] = ismember(deg(:, posS(i)), Zs{i});
        if any(~ok), error('polynomial_to_quadPoly:degreeMissingS','Degree missing in Zs at %s', sNames{i}); end
        ks = ks + (loc-1) * sStride(i);
    end

    for i = 1:numel(tNames)
        if posT(i) == 0
            continue;
        end
        [ok, loc] = ismember(deg(:, posT(i)), Zt{i});
        if any(~ok), error('polynomial_to_quadPoly:degreeMissingT','Degree missing in Zt at %s', tNames{i}); end
        kt = kt + (loc-1) * tStride(i);
    end

    % Now place coefficients into Q.C:
    % Q.C row = (i-1)*ds + ks(term)
    % Q.C col = (j-1)*dt + kt(term)
    %
    % We build by accumulating sparse triplets.
    [termIdx, entryIdx, val] = find(coeff); % coeff(term, entry)

    % entryIdx corresponds to matrix entry col in polynomial format:
    % entry = i + m*(j-1)  => i = mod(entry-1,m)+1, j = floor((entry-1)/m)+1
    iMat = mod(entryIdx-1, m) + 1;
    jMat = floor((entryIdx-1)/m) + 1;

    row = (iMat-1) * ds + ks(termIdx);
    col = (jMat-1) * dt + kt(termIdx);

    C = sparse(row, col, val, m*ds, n*dt);

    % Construct quadPoly (property assignment is fine if no custom ctor)
    Q = quadPoly(C,Zs,Zt,[m,n],sNames(:).', tNames(:).');
end

% ---------------- helpers ----------------
function d = basisSize(Zcell)
% D = BASISSIZE(ZCELL) returns the number of monomials of the tensor
% product basis defined by the exponent vectors in ZCELL (an empty
% element counts as the constant monomial).
%
% INPUT
% - Zcell:  1 x k cell of exponent vectors;
%
% OUTPUT
% - d:      positive integer, the product of their lengths (1 if empty).
    if isempty(Zcell), d = 1; return; end
    d = 1;
    for k = 1:numel(Zcell)
        z = Zcell{k};
        if isempty(z), z = 0; end
        d = d * numel(z);
    end
end

function stride = kronStrides(sizes)
% STRIDE = KRONSTRIDES(SIZES) returns the strides of the linear index of a
% tensor-product basis with factor sizes SIZES,
%   stride(i) = prod(sizes(i+1:end)).
%
% INPUT
% - sizes:  1 x k array of positive integers;
%
% OUTPUT
% - stride: 1 x k array, with stride(k) = 1.
    % last variable varies fastest
    stride = cumprod(sizes,2,"reverse");
    stride = [stride(2:end),1];
end
