function [K,Kop,info] = getController(P,Z,opts)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% [K,KOP,INFO] = GETCONTROLLER(P,Z[,OPTS]) returns the state-feedback
% gain K, u = K*x, x in R^k x L2^m[s], of the solved H-infinity or H2
% controller synthesis LPI, an approximation of
%   K = Z*P^{-1}:
% by default (method 'qp') the K of bounded degree that minimizes
% ||Z - K*P||^2, with no inverse of P; with method 'inv', K = Z*P^{-1}, the
% inverse the Gohberg-Krein inverse of the L2 block, its kernels fitted by
% polynomials, and the finite-dimensional Schur complement for an R^k
% space, and the product the composition of the containers. The
% coefficients of K are not truncated.
%
% INPUT
% - P:      'opvar' object, the
%           solved Lyapunov operator, with the same spaces and dimensions in
%           and out. P may be an N-D opvar, in any number of spatial
%           variables and with any spaces; the inverse (method 'inv') exists
%           in one spatial variable for R^k x L2^m[s] with at most one R^k
%           and one L2 space, and in several for a separable P; otherwise
%           the method 'qp' is required;
% - Z:      'opvar' object, the solved free operator, from the
%           spaces of P to R^nu, R^k x L2^m[s] -> R^nu in 1-D; for the method
%           'qp' also N-D (the spatial variables of P and Z);
% - opts:   (optional) 'struct', each field optional. The options of the
%           inverse: N (nodes of the RK4 grid, default 101), tol (target
%           relative RMS residual of the polynomial fit, 1e-8), deg0 (first
%           fit degree, 4), degmax (largest fit degree, 16), deg (one fixed
%           fit degree, in place of deg0, degmax, tol), tol_rank (relative
%           SVD cut of the kernel ranks, 1e-13). And the method:
%             method     'qp' (default), the gain K of bounded degree
%                        minimizing ||Z - K*P||^2 by a quadratic program in
%                        the coefficients of K, with no inverse of P (P must
%                        be self-adjoint), for any number of spatial
%                        variables and required for a non-separable P in
%                        several; or 'inv', the gain K = Z*P^{-1} above (in
%                        several variables an error for a non-separable P).
%                        The options of the 'qp' method: qp_deg (maximal
%                        degree of K in the spatial variable), qp_basis
%                        ('monomial', 'legendre' or 'chebyshev'), qp_assembly
%                        ('container' or 'apply'), qp_solver ('quadprog' or
%                        'backslash'), qp_lambda (>= 0, Tikhonov term),
%                        qp_verbose; each unset one at the default of the
%                        quadratic program.
%           A number, the TOL of the PIETOOLS 2025 call getController(P,Z,
%           tol), stops with the identifier getController:call2025.
%
% OUTPUT
% - K:      'opvar' object, R^k x L2^m[s] -> R^nu;
% - Kop:    [], kept so that [K,~,INFO] calls keep the position of INFO; the
%           gain is K;
% - info:   for the method 'inv', 'struct' of the inverse, with fields N, h
%           (grid step), r1, r2 (kernel ranks), n (= r1 + r2), rcondU22,
%           condR0 (largest condition number of R0 on the grid), d (fit
%           degree), relrms (1x3 relative RMS residuals of the fits), note,
%           and condT (condition number of the Schur complement) when P has
%           both an R^k and an L2 space; a 'struct' with no fields when P has
%           no L2 space. For the method 'qp', a 'struct' with the field qp,
%           a 1 x nu cell, the diagnostics of the quadratic program for each
%           of the nu columns of Z' (without its matrices H and f).
%
% See also GETOBSERVER, OPVAR/INV, QP_OPVAR, CLOSEDLOOPPIE.
%
% For support, contact M. Peet, Arizona State University at mpeet@asu.edu

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PIETOOLS - getController
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
% Initial coding MMP, 10/09/2026
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% A number is the TOL of the PIETOOLS 2025 call getController(P,Z,tol), which
% isfield below would ignore: stop with the current call
if nargin>=3 && isnumeric(opts) && ~isempty(opts)
    error('getController:call2025',['getController(P,Z,TOL) is the call of ' ...
        'PIETOOLS 2025; the third input is now a struct of options. Call ' ...
        '[K,Kop,info] = getController(P,Z) for the method ''qp'', or ' ...
        'getController(P,Z,struct(''method'',''inv'')) for K = Z*P^{-1}; ' ...
        'no coefficient of K is truncated.']);
end
% Default options
if nargin<3,    opts = struct();    end
% The method: 'qp' (default) or 'inv'
method = 'qp';
if isfield(opts,'method'),  method = opts.method;   end
if ~any(strcmpi(method,{'inv','qp'}))
    error('getController:method','opts.method must be ''inv'' or ''qp''.');
end
if strcmpi(method,'qp')
    % K = argmin ||Z - K*P||^2 over K of bounded degree. P is self-adjoint, so
    % K' = argmin ||Z' - P K'||^2: for each column j of Z' (input R^nu) a least
    % squares problem with A = P and B = Z' e_j, a B with one input. No
    % inverse of P is formed. The degree of K is bounded by opts.qp_deg; the
    % true K is not a polynomial, which limits the accuracy.
    qo = struct();
    for f = {'deg','basis','assembly','solver','lambda','verbose'}
        if isfield(opts,['qp_' f{1}]),  qo.(f{1}) = opts.(['qp_' f{1}]);    end
    end
    nu = Z.dim_out(1);
    Zt = Z';
    Xt = [];    info = struct('qp',{cell(1,nu)});
    for j = 1:nu
        % column selector e_j: R^1 -> R^nu
        % as a container on R, with the variable registry of Z
        Ej = mat2opvar(full(sparse(j,1,1,nu,1)),struct('out',nu,'in',1),...
                 struct('out',{{{}}},'in',{{{}}}),struct('vars',{Z.vars},'dom',Z.dom));
        [Xj,info.qp{j}] = qp_opvar(P,Zt*Ej,qo);
        info.qp{j}.H = [];   info.qp{j}.f = [];     % keep the diagnostics, not the n x n matrix
        Ejt = Ej';
        if j==1,    Xt = Xj*Ejt;    else,   Xt = Xt + Xj*Ejt;   end
    end
    K = Xt';                                % K = (sum_j X_j e_j')'
else
    % In several variables the inverse exists for a separable P
    % K = Z*P^{-1}
    [Pinv,info] = inv(P,opts);
    K = Z*Pinv;
end

% No operator copy of K: the gain is the container K
Kop = [];
end
