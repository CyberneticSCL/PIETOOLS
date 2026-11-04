function PIEcl = closedLoopPIE_OF_sop(PIE,Kc)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PIEcl = CLOSEDLOOPPIE_OF_SOP(PIE,KC) returns the closed loop of the PIE
%
%       T x_t = A x + B1 w + B2 u,
%       z = C1 x + D11 w + D12 u,       y = C2 x + D21 w + D22 u,
%
% with the dynamic output feedback controller
%
%       TK xK_t = AK xK + BK y,         u = CK xK + DK y,
%
% as the PIE in the state [x; xK],
%
%   [T 0; 0 TK] [x; xK]_t = [A + B2 DK2 C2, B2 CK2; BK2 C2, AK2] [x; xK]
%                           + [B1 + B2 DK2 D21; BK2 D21] w,
%   z = [C1 + D12 DK2 C2, D12 CK2] [x; xK] + (D11 + D12 DK2 D21) w,
%
% where, for D22 nonzero, Q = (I - D22 DK)^{-1}, AK2 = AK + BK Q D22 CK,
% BK2 = BK Q, CK2 = (I + DK Q D22) CK, DK2 = DK Q; only D22 = 0 is implemented, where
% these are AK, BK, CK, DK. The control input u and the observed output y
% are eliminated.
%
% INPUT
% - PIE:    'pie_struct' object of the plant, with 'copvar' or 'opvar'
%           operators, and Tw, Tu zero or absent (no algebraic constraint);
% - Kc:     'struct' with fields TK, AK, BK, CK, DK, the operators of the
%           controller as 'copvar' or 'opvar' objects, as returned by
%           PIETOOLS_HINF_OFCONTROL_SOP (Kc, or Kc.G);
%
% OUTPUT
% - PIEcl:  initialized 'pie_struct' object with 'copvar' operators, of the
%           state [x; xK], the exogenous input w and the regulated output z.
%           The maps to the PDE are removed, as the PIE no longer
%           corresponds to the original PDE;
%
% For support, contact M. Peet, Arizona State University at mpeet@asu.edu

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PIETOOLS - closedLoopPIE_OF_sop
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
% Initial coding: Claude 10/31/2026, in the structure of CLOSEDLOOPPIE for the
%                observer closed loop.
% 11/03/2026: Superseded by CLOSEDLOOPPIE(PIE,Kc,'dynamic'); moved to cleaning.

if ~isa(PIE,'pie_struct')
    error('closedLoopPIE_OF_sop:pie','The PIE should be a ''pie_struct'' object.');
end
if isfield(Kc,'G'),     Kc = Kc.G;  end
PIE = initialize(PIE);
if size(PIE.K,1)>0
    error('closedLoopPIE_OF_sop:K','PIEs with an algebraic constraint (K, Kw, Ku) are not supported.');
end

% The controller operators as containers
for f = {'TK','AK','BK','CK','DK'}
    if isa(Kc.(f{1}),'opvar') || isa(Kc.(f{1}),'opvar2d'),  Kc.(f{1}) = op2copvar_sop(Kc.(f{1}));   end
end
TK = Kc.TK;     AK = Kc.AK;     BK = Kc.BK;     CK = Kc.CK;     DK = Kc.DK;

% D22 = 0 (AK2 = AK, BK2 = BK, CK2 = CK, DK2 = DK)
if ~(PIE.Dyu==0)
    error('closedLoopPIE_OF_sop:D22','The closed loop is implemented for D22 = 0.');
end

% The closed-loop PIE in the state [x; xK]
PIEcl = pie_struct();
PIEcl.T  = blkdiag(PIE.T,TK);
PIEcl.A  = [PIE.A+PIE.Bu*DK*PIE.Cy,     PIE.Bu*CK;
            BK*PIE.Cy,                  AK];
PIEcl.Bw = [PIE.Bw+PIE.Bu*DK*PIE.Dyw;   BK*PIE.Dyw];
PIEcl.Cz = [PIE.Cz+PIE.Dzu*DK*PIE.Cy,   PIE.Dzu*CK];
PIEcl.Dzw = PIE.Dzw + PIE.Dzu*DK*PIE.Dyw;

% The closed loop no longer corresponds to the original PDE, so the maps to
% the PDE are removed, and the PIE is initialized
default_pie = pie_struct();
PIEcl.space_pde = default_pie.space_pde;
PIEcl.pie2pde = default_pie.pie2pde;
PIEcl.pde2pie = default_pie.pde2pie;
PIEcl = initialize(PIEcl);
end
