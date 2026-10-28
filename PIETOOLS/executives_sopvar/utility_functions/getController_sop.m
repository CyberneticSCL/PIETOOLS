function [K,Kop,info] = getController_sop(P,Z,opts)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% [K,KOP,INFO] = GETCONTROLLER_SOP(P,Z[,OPTS]) returns the state-feedback
% gain
%   K = Z*P^{-1},   u = K*x,   x in R^k x L2^m[s],
% of the solved H-infinity or H2 controller synthesis LPI, with P^{-1} the
% inverse of P by COPVAR/INV (the Gohberg-Krein inverse of the L2 block,
% its kernels fitted by polynomials, and the finite-dimensional Schur
% complement for an R^k space) and the product the composition of the
% containers. The coefficients of K are not truncated.
%
% INPUT
% - P:      'copvar' object (or 'opvar', converted by OPVAR2COPVAR), the
%           solved Lyapunov operator: the same spaces and dimensions in and
%           out, R^k x L2^m[s] in one spatial variable, at most one R^k and
%           one L2 space;
% - Z:      'copvar' object (or 'opvar'), the solved free operator
%           R^k x L2^m[s] -> R^nu;
% - opts:   (optional) 'struct', the options of the inverse, each field
%           optional: N (nodes of the RK4 grid, default 101), tol (target
%           relative RMS residual of the polynomial fit, 1e-8), deg0 (first
%           fit degree, 4), degmax (largest fit degree, 16), deg (one fixed
%           fit degree, in place of deg0, degmax, tol), tol_rank (relative
%           SVD cut of the kernel ranks, 1e-13).
%
% OUTPUT
% - K:      'copvar' object, R^k x L2^m[s] -> R^nu;
% - Kop:    'opvar' object, K converted by COPVAR2OPVAR ([] when fewer than
%           two outputs are requested);
% - info:   'struct' of the inverse, with fields N, h (grid step), r1, r2
%           (kernel ranks), n (= r1 + r2), rcondU22, condR0 (largest
%           condition number of R0 on the grid), d (fit degree), relrms
%           (1x3 relative RMS residuals of the fits), note, and condT
%           (condition number of the Schur complement) when P has both an
%           R^k and an L2 space; a 'struct' with no fields when P has no L2
%           space.
%
% See also GETOBSERVER_SOP, COPVAR/INV, GETCONTROLLER, CLOSEDLOOPPIE.
%
% For support, contact M. Peet, Arizona State University at mpeet@asu.edu

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PIETOOLS - getController_sop
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
%                the stock getController and the uses of Kop elsewhere.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Default options; 'opvar' operators to containers
if nargin<3,    opts = struct();    end
if isa(P,'opvar'),  P = opvar2copvar(P);    end
if isa(Z,'opvar'),  Z = opvar2copvar(Z);    end
% K = Z*P^{-1}
[Pinv,info] = inv(P,opts);
K = Z*Pinv;
% The 'opvar' form, when requested
Kop = [];
if nargout>1,   Kop = copvar2opvar(K);  end
end
