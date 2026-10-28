function [L,Lop,info] = getObserver_sop(P,Z,opts)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% [L,LOP,INFO] = GETOBSERVER_SOP(P,Z[,OPTS]) returns the observer gain
%   L = P^{-1}*Z
% of the solved H-infinity or H2 estimator LPI, with P^{-1} the inverse of
% P by COPVAR/INV (the Gohberg-Krein inverse of the L2 block, its kernels
% fitted by polynomials, and the finite-dimensional Schur complement for an
% R^k space) and the product the composition of the containers.
%
% INPUT
% - P:      'copvar' object (or 'opvar', converted by OPVAR2COPVAR), the
%           solved Lyapunov operator, with the same spaces and dimensions in
%           and out. P may be an N-D copvar, in any number of spatial
%           variables and with any spaces; but the inverse (method 'inv')
%           exists for one spatial variable only, R^k x L2^m[s] with at most
%           one R^k and one L2 space. If P or Z depends on more than one
%           spatial variable, the method 'qp' is required;
% - Z:      'copvar' object (or 'opvar'), the solved free operator, from R^ny
%           to the spaces of P, R^ny -> R^k x L2^m[s] in 1-D; for the method
%           'qp' also N-D (the spatial variables of P and Z);
% - opts:   (optional) 'struct', each field optional. The options of the
%           inverse: N (nodes of the RK4 grid, default 101), tol (target
%           relative RMS residual of the polynomial fit, 1e-8), deg0 (first
%           fit degree, 4), degmax (largest fit degree, 16), deg (one fixed
%           fit degree, in place of deg0, degmax, tol), tol_rank (relative
%           SVD cut of the kernel ranks, 1e-13). And the method:
%             method     'inv' (default), the gain L = P^{-1}*Z above, for
%                        one spatial variable only (an error for an N-D P or
%                        Z), or 'qp', the gain L of bounded degree minimizing
%                        ||P*L - Z||^2 by QP_COPVAR, with no inverse of P
%                        (P must be self-adjoint), for any number of spatial
%                        variables and required when there are several. The options of the 'qp'
%                        method, passed to QP_COPVAR: qp_deg (maximal degree
%                        of L in the spatial variable), qp_basis
%                        ('monomial', 'legendre' or 'chebyshev'), qp_assembly
%                        ('container' or 'apply'), qp_solver ('quadprog' or
%                        'backslash'), qp_lambda (>= 0, Tikhonov term),
%                        qp_verbose; each defaults to that of QP_COPVAR.
%
% OUTPUT
% - L:      'copvar' object, R^ny -> R^k x L2^m[s];
% - Lop:    'opvar' object, L converted by COPVAR2OPVAR ([] when fewer than
%           two outputs are requested, or, with a warning, when L depends on
%           more than one spatial variable: there is no 'opvar' form then);
% - info:   for the method 'inv', 'struct' of the inverse, with fields N, h
%           (grid step), r1, r2 (kernel ranks), n (= r1 + r2), rcondU22,
%           condR0 (largest condition number of R0 on the grid), d (fit
%           degree), relrms (1x3 relative RMS residuals of the fits), note,
%           and condT (condition number of the Schur complement) when P has
%           both an R^k and an L2 space; a 'struct' with no fields when P has
%           no L2 space. For the method 'qp', a 'struct' with the field qp,
%           a 1 x ny cell, the INFO of QP_COPVAR for each of the ny columns
%           of Z (without its matrices H and f).
%
% See also GETCONTROLLER_SOP, COPVAR/INV, QP_COPVAR, GETOBSERVER, GETOBSERVER_2D.
%
% For support, contact M. Peet, Arizona State University at mpeet@asu.edu

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PIETOOLS - getObserver_sop
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
% Initial coding MMP, 10/09/2026.
% MMP, 10/28/2026: Header to the header protocol (what the function does,
%                inputs, outputs), step comments in the body; code
%                unchanged. The deleted header text was a comparison with
%                the stock getObserver, the uses of Lop elsewhere and a
%                reference to the 2-D estimator.
% AT, 10/28/2026: Option opts.method = 'qp': L is the minimizer of ||P*L - Z||^2
%                over L of bounded degree, from QP_COPVAR (a QP in the
%                coefficients of L), one QP per column of Z; no inverse of P is
%                formed. The previous computation is the method 'inv', the
%                default, and is unchanged. Options opts.qp_deg, qp_basis,
%                qp_assembly, qp_solver, qp_lambda, qp_verbose. Lop is the
%                'opvar' of L for at most one spatial variable, otherwise
%                [] with a warning (as Kop of GETCONTROLLER_SOP).
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Default options; 'opvar' operators to containers
if nargin<3,    opts = struct();    end
if isa(P,'opvar'),  P = opvar2copvar(P);    end
if isa(Z,'opvar'),  Z = opvar2copvar(Z);    end
% The method: 'inv' (default) or 'qp'
method = 'inv';                                                             % AT, 10/28/2026
if isfield(opts,'method'),  method = opts.method;   end                     % AT, 10/28/2026
if ~any(strcmpi(method,{'inv','qp'}))                                       % AT, 10/28/2026
    error('getObserver_sop:method','opts.method must be ''inv'' or ''qp''.');     % AT, 10/28/2026
end                                                                         % AT, 10/28/2026
if strcmpi(method,'qp')                                                     % AT, 10/28/2026
    % L = argmin ||P*L - Z||^2 over L of bounded degree: for each column j of Z % AT, 10/28/2026
    % (input R^ny) a least squares problem for QP_COPVAR with A = P and        % AT, 10/28/2026
    % B = Z e_j, since it takes a B with one input. No inverse of P is formed.  % AT, 10/28/2026
    % The degree of L is bounded by opts.qp_deg; the true L is not a polynomial, % AT, 10/28/2026
    % which limits the accuracy.                                                % AT, 10/28/2026
    qo = struct();                                                          % AT, 10/28/2026
    for f = {'deg','basis','assembly','solver','lambda','verbose'}          % AT, 10/28/2026
        if isfield(opts,['qp_' f{1}]),  qo.(f{1}) = opts.(['qp_' f{1}]);    end % AT, 10/28/2026
    end                                                                     % AT, 10/28/2026
    ny = Z.dim_in(1);                                                       % AT, 10/28/2026
    pvar s_qp s_qp_dum                                                      % AT, 10/28/2026
    Lc = [];    info = struct('qp',{cell(1,ny)});                           % AT, 10/28/2026
    for j = 1:ny                                                            % AT, 10/28/2026
        % column selector e_j: R^1 -> R^ny                                  % AT, 10/28/2026
        Ej = opvar2copvar(mat2opvar(full(sparse(j,1,1,ny,1)),[ny 1;0 0],[s_qp,s_qp_dum],[0,1])); % AT, 10/28/2026
        [Xj,info.qp{j}] = qp_copvar(P,Z*Ej,qo);                             % AT, 10/28/2026
        info.qp{j}.H = [];   info.qp{j}.f = [];     % keep the diagnostics, not the n x n matrix  % AT, 10/28/2026
        Ejt = Ej';                                                          % AT, 10/28/2026
        if j==1,    Lc = Xj*Ejt;    else,   Lc = Lc + Xj*Ejt;   end         % AT, 10/28/2026
    end                                                                     % AT, 10/28/2026
    L = Lc;                                 % L = sum_j X_j e_j'             % AT, 10/28/2026
else                                                                        % AT, 10/28/2026
    % The inverse is one-dimensional                                        % AT, 10/28/2026
    if numel(P.vars)>1 || numel(Z.vars)>1                                   % AT, 10/28/2026
        error('getObserver_sop:ndInverse',['P and Z depend on %d spatial variables, and the inverse '... % AT, 10/28/2026
            '(method ''inv'') is one-dimensional: use opts.method = ''qp''.'],max(numel(P.vars),numel(Z.vars)));% AT, 10/28/2026
    end                                                                     % AT, 10/28/2026
    % L = P^{-1}*Z
    [Pinv,info] = inv(P,opts);
    L = Pinv*Z;
end                                                                         % AT, 10/28/2026

Lop = [];                                                                   % AT, 10/28/2026
if nargout>1                                                                % AT, 10/28/2026
    if numel(L.vars)<=1                                                     % AT, 10/28/2026
        Lop = copvar2opvar(L);                                              % AT, 10/28/2026
    else    % no 'opvar' form in N-D: Lop = []                              % AT, 10/28/2026
        warning('getObserver_sop:Lop','L depends on %d spatial variables and has no ''opvar'' form: Lop = [].',numel(L.vars));% AT, 10/28/2026
    end                                                                     % AT, 10/28/2026
end                                                                         % AT, 10/28/2026
end
