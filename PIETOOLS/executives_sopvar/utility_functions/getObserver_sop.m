function [L,Lop,info] = getObserver_sop(P,Z,opts)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% [L,LOP,INFO] = GETOBSERVER_SOP(P,Z[,OPTS]) the observer gain
%   L = P^{-1} Z
% of the H-infinity / H2 estimator LPIs on the container path, the
% counterpart of the stock 'getObserver' (inv(P,tol)*Z). P^{-1} is
% COPVAR/INV; the product is the container composition. 1-D only (the 2-D
% estimator keeps 'getObserver_2D').
%
% INPUTS
% - P:      the solved Lyapunov operator, 'copvar' (or 'opvar', converted by
%           OPVAR2COPVAR), square on R^k x L2^m[s];
% - Z:      the solved free operator, 'copvar' (or 'opvar'), R^ny -> R^k x L2^m;
% - opts:   (optional) the options of SOPVAR/INV.
% OUTPUTS
% - L:      'copvar', R^ny -> R^k x L2^m;
% - Lop:    the same gain as an 'opvar' (COPVAR2OPVAR), for 'closedLoopPIE'
%           and PIESIM;
% - info:   the INFO of COPVAR/INV.
%
% See also GETCONTROLLER_SOP, COPVAR/INV, GETOBSERVER, GETOBSERVER_2D.
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
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if nargin<3,    opts = struct();    end
if isa(P,'opvar'),  P = opvar2copvar(P);    end
if isa(Z,'opvar'),  Z = opvar2copvar(Z);    end
[Pinv,info] = inv(P,opts);
L = Pinv*Z;
Lop = [];
if nargout>1,   Lop = copvar2opvar(L);  end
end
