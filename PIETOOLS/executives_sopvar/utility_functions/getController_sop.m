function [K,Kop,info] = getController_sop(P,Z,opts)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% [K,KOP,INFO] = GETCONTROLLER_SOP(P,Z[,OPTS]) the state-feedback gain
%   K = Z P^{-1}
% of the dual H-infinity / H2 synthesis LPIs on the container path, the
% counterpart of the stock 'getController' (K = Z*inv(P,tol), then
% clean_opvar). P^{-1} is COPVAR/INV (Gohberg-Krein on the L2 block, the
% finite-dimensional Schur complement for an R^k space); the product is the
% container composition. Nothing is truncated: the stock routine zeroes
% every coefficient of K below tol = 1e-4, which on a monomial basis is not
% a bound on the operator.
%
% INPUTS
% - P:      the solved Lyapunov operator, 'copvar' (or 'opvar', converted by
%           OPVAR2COPVAR), square on R^k x L2^m[s];
% - Z:      the solved free operator, 'copvar' (or 'opvar'), R^k x L2^m -> R^nu;
% - opts:   (optional) the options of SOPVAR/INV.
% OUTPUTS
% - K:      'copvar', R^k x L2^m -> R^nu;
% - Kop:    the same gain as an 'opvar' (COPVAR2OPVAR), for 'closedLoopPIE',
%           'piess' and PIESIM, which take the legacy classes only;
% - info:   the INFO of COPVAR/INV.
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
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if nargin<3,    opts = struct();    end
if isa(P,'opvar'),  P = opvar2copvar(P);    end
if isa(Z,'opvar'),  Z = opvar2copvar(Z);    end
[Pinv,info] = inv(P,opts);
K = Z*Pinv;
Kop = [];
if nargout>1,   Kop = copvar2opvar(K);  end
end
