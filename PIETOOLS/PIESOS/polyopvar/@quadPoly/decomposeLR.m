function [L, R, info] = decomposeLR(F, tol)
% DECOMPOSELR  Decompose F(s,t) into L(s)*R(t) with minimal inner dim.
%
%   F(s,t) = (I_m ⊗ Z(s)^T) * C * (I_n ⊗ Z(t))
%   C = U*V  (rank r, as small as possible numerically)
%   L(s) = (I_m ⊗ Z(s)^T) * U     depends only on s
%   R(t) = V * (I_n ⊗ Z(t))       depends only on t
%
% INPUTS
%   F       m x n quadPoly object
%   tol     (optional) relative tolerance for the numerical rank:
%           singular values below tol*max(svals) are dropped
%           (default 1e-14)
%
% OUTPUTS
%   L, R    m x r and r x n quadPoly objects such that F = L*R (in the
%           quadratic-form sense)
%   info    'struct' with fields r (the rank), svals (the singular values
%           computed), relerr_est (the relative Frobenius norm of C-U*V)
%           and method ('svd(full)' or 'svds(grow)')
%
% NOTES
%   - For L(s) we set the t-side to constant monomial 1 (Zt = {0},
%     nt = {}).
%   - For R(t) we set the s-side to constant monomial 1 (Zs = {0},
%     ns = {}).

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PIETOOLS - decomposeLR
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

    if nargin < 2 || isempty(tol),     tol = 1e-14; end

    C = F.C;
    [M, N] = size(C);

    maxRank = min(M, N);

    % --- Compute low-rank factorization of C with minimal numerical rank ---
    % Use full svd for moderate sizes; otherwise svds with a growing k.
    useFullSVD = (maxRank <= 2000) && (numel(C) <= 5e6);  % simple heuristic

    if useFullSVD
        info.method = 'svd(full)';
        [U0, S0, V0] = svd(full(C), 'econ');
        s = diag(S0);
    else
        info.method = 'svds(grow)';
        % Grow k until tail singular values fall below tol*max(s)
        k = min([50, maxRank]);   % start
        s = [];
        U0 = []; S0 = []; V0 = [];
        while true
            [Uk, Sk, Vk] = svds(C, k);   % largest k singular values
            sk = diag(Sk);
            sk = sort(sk, 'descend');    % svds may not return sorted
            smax = sk(1);

            % Determine numerical rank within the computed spectrum
            rk = find(sk >= tol*smax, 1, 'last');
            if isempty(rk), rk = 0; end

            % If rk < k or we've hit maxRank, accept.
            if rk < k || k >= maxRank
                % Recompute once with k = max(rk,1) to get consistent factors
                rfinal = max(rk, 1);
                [U0, S0, V0] = svds(C, rfinal);
                s = diag(S0);
                s = sort(s, 'descend');
                break;
            end

            % otherwise increase k
            kNew = min([2*k, maxRank]);
            if kNew == k
                break;
            end
            k = kNew;
        end
    end

    if isempty(s)
        % In case something odd happens; treat as zero
        r = 0;
        s = 0;
        U = sparse(M,0);
        V = sparse(0,N);
    else
        smax = max(s);
        r = find(s >= tol*smax, 1, 'last');
        if isempty(r), r = 0; end

        if r == 0
            U = sparse(M,0);
            V = sparse(0,N);
        else
            % Build U*V using symmetric sqrt(S) split for better conditioning:
            % C ≈ (U0(:,1:r)*sqrt(S)) * (sqrt(S)*V0(:,1:r)')
            if useFullSVD
                U0r = U0(:,1:r);
                V0r = V0(:,1:r);
                sr  = s(1:r);
            else
                % svds returned exactly rfinal; keep top-r (already small)
                % sort not aligned with columns, but for typical svds it's OK.
                % If you want strict alignment, skip sorting and use diag(S0).
                U0r = U0(:,1:r);
                V0r = V0(:,1:r);
                sr  = diag(S0);
                sr  = sr(1:r);
            end

            sqrtSr = spdiags(sqrt(sr(:)), 0, r, r);
            U = sparse(U0r) * sqrtSr;        % (M x r)
            V = sqrtSr * sparse(V0r)';       % (r x N)
        end
    end

    % --- Package as quadPoly objects ---
    m = F.dim(1);
    n = F.dim(2);

    % L(s): dim = [m, r], depends on s only
    L = quadPoly(U,F.Zs,{0},[m,r],F.ns,{});

    % R(t): dim = [r, n], depends on t only
    R = quadPoly(V,{0},F.Zt,[r,n],{},F.nt);

    % --- Diagnostics ---
    info.r = r;
    info.svals = s(:);
    if r > 0
        % crude relative reconstruction error estimate using Fro norm
        denom = norm(C, 'fro');
        if denom > 0
            info.relerr_est = norm(C - U*V, 'fro') / denom;
        else
            info.relerr_est = 0;
        end
    else
        info.relerr_est = norm(C, 'fro') / max(norm(C, 'fro'), 1);
    end
end
