%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% HINF_GAIN_1D_SOP computes an upper bound gam on the H-infinity gain from
% w to z of the 1-D PDE
%
%       x_t = x_ss + (pi^2/2) x + s(1-s) w,   s in [0,1],   x(0) = x(1) = 0,
%       z   = int_0^1 x(s) ds,
%
% by an LPI on the container path. The PDE is converted to the PIE
% T v_t = A v + B w, z = C v + D w, its operators to 'copvar' objects, and
% the script minimizes gam subject to P = P' >= 0 and
%
%       [ -gam*I     D'       B'*P*T         ]
%       [  D        -gam*I    C              ]  <= 0,
%       [  T'*P*B    C'       T'*P*A + A'*P*T]
%
% where, in the code, the scalar -gam in a block matrix stands for -gam
% times the identity of its block. POSCOPVAR declares P on L2[s] at
% degree 2, and
% LPI_INEQ_SOP imposes the inequality with a certificate sized from P,
% printing its decision variables first. The script then computes the
% bound of PIETOOLS_Hinf_gain_sop(PIE,lpisettings('light')) and prints
% both bounds.
%
% INPUT
% none.
%
% OUTPUT
% the workspace variables
% - PDE:            'pde_struct' object, the initialized PDE;
% - PIE:            'pie_struct' object, its PIE;
% - T, A, B, C, D:  'copvar' objects, PIE.T, PIE.A, PIE.B1, PIE.C1 and
%                   PIE.D11;
% - prog:           'struct', the solved LPI program;
% - gam:            'dpvar' object, the decision variable gam;
% - P, K:           'cdopvar' objects, the positive operator and the block
%                   operator of the LPI;
% - gain:           double, the bound of the LPI;
% - settings:       'struct', lpisettings('light');
% - gain_exec:      double, the bound of PIETOOLS_Hinf_gain_sop;
% - s, t:           'polynomial' objects, the spatial and time variables;
% - x, w, z:        'pde_struct' objects, the state, input and output.
%
% For support, contact M. Peet, Arizona State University at mpeet@asu.edu

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PIETOOLS - hinf_gain_1d_sop
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
% Initial coding MMP, 10/21/2026. The end-to-end check that had this name
% is lpi_programming_sopvar/tests/e2e_hinf_gain_1d.m.
% MMP, 10/21/2026: lpi_ineq_sop reads P from the program (poscopvar records
%                it), so the options are gone. Deleted (same date):
%                "ineq_opts.like = P;" and the ineq_opts argument of
%                lpi_ineq_sop, with the help text on the option. Kept out
%                of the body as (was) lines, which echo would print.
% MMP, 10/21/2026: No solver options: lpisolve without them takes the first
%                solver installed, MOSEK and then SeDuMi (sossolve), as
%                the executive does from lpisettings (exec_sop_settings);
%                options belong to the settings files (maintainer).
%                Deleted (same date): "% MOSEK when it is installed, SeDuMi
%                otherwise", "sos_opts = struct('solver','sedumi');", "if
%                ~isempty(which('mosekopt')), sos_opts.solver = 'mosek';
%                end", the sos_opts argument of lpisolve and
%                "settings.sos_opts = sos_opts;".
% MMP, 10/28/2026: Header to the header protocol (what the script
%                computes, inputs, outputs), the licence block added,
%                step comments in the body; code unchanged. The deleted
%                header text was the reference to the stock executive of
%                the LPI and its folder. The help text of the 10/21/2026
%                entries is in the paragraph after the LPI, without
%                markers.

echo on
% The PDE and its PIE
pvar s t
x = pde_var('state',1,s,[0,1]);
w = pde_var('in',1);
z = pde_var('out',1);
PDE = [diff(x,t) == diff(x,s,2) + (pi^2/2)*x + s*(1-s)*w;
       z == int(x,s,[0,1]);
       subs(x,s,0) == 0;
       subs(x,s,1) == 0];
PDE = initialize(PDE);
PIE = convert(PDE);
% The PIE operators as containers
T = op2copvar_sop(PIE.T);   A = op2copvar_sop(PIE.A);
B = op2copvar_sop(PIE.B1);  C = op2copvar_sop(PIE.C1);  D = op2copvar_sop(PIE.D11);

% The LPI
prog = lpiprogram_sop(T);
[prog,gam] = lpidecvar(prog,'gam');
[prog,P] = poscopvar(prog,1,{{'s'}},[0,1],2);       % P >= 0 on L2[s], degree 2
K = [-gam,       D',     B'*P*T;
      D,         -gam,   C;
      T'*P*B,    C',     T'*P*A + A'*P*T];
prog = lpi_ineq_sop(prog,-K);                       % K <= 0                % MMP, 10/21/2026
prog = lpisetobj(prog,gam);                         % minimize gam
prog = lpisolve(prog);                              % MOSEK, else SeDuMi    % MMP, 10/21/2026

% The solved bound
gain = double(lpigetsol_sop(prog,gam));

% The same problem through the executive
settings = lpisettings('light');
[~,~,gain_exec] = PIETOOLS_Hinf_gain_sop(PIE,settings);
echo off

% Print the two bounds
fprintf('\n H-infinity gain bound: %.6f by the LPI above, %.6f by PIETOOLS_Hinf_gain_sop.\n',...
        gain,gain_exec);
