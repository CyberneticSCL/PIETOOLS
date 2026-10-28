%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% STABILITY_2D_SOP tests stability of the 2-D PDE
%
%       x_t = x_{s1 s1} + x_{s2 s2} + lam*x  on [0,1]^2,  x = 0 on the
%       boundary,  lam = pi^2,
%
% by an LPI on the container path. The PDE is stable since lam < 2 pi^2.
% Its PIE is T v_t = A v, and the PIE is stable if, for some eps > 0,
% there exists P = P' >= eps*I such that
%
%       T'*P*A + A'*P*T <= 0.
%
% The PIE operators are converted to 'copvar' objects by OP2COPVAR_SOP;
% POSCOPVAR declares a positive operator on L2[s1,s2] at degree 1, P is
% that operator plus 1e-3*I, and LPI_INEQ_SOP imposes the inequality with
% a certificate sized from the operator of POSCOPVAR, the sum of a plain
% term and one term per face of [0,1]^2, printing its decision variables
% first. The script prints that the PIE is stable when the solver reports
% pinf = 0, dinf = 0 and numerr = 0, else those three values.
% PIETOOLS_stability_sop(PIE,settings,'P') poses this LPI, with the
% degrees and margin of SETTINGS, in one call.
%
% INPUT
% none.
%
% OUTPUT
% the workspace variables
% - lam:       double, pi^2;
% - PDE:       'pde_struct' object, the initialized PDE;
% - PIE:       'pie_struct' object, its PIE;
% - T, A:      'copvar' objects, PIE.T and PIE.A;
% - dom:       double, [0,1;0,1], the domain;
% - prog:      'struct', the solved LPI program;
% - P:         'cdopvar' object, the Lyapunov operator (P >= 1e-3*I);
% - info:      'struct', prog.solinfo.info, the solver report (fields
%              pinf, dinf, numerr);
% - s1, s2, t: 'polynomial' objects, the spatial and time variables;
% - x:         'pde_struct' object, the state.
%
% For support, contact M. Peet, Arizona State University at mpeet@asu.edu

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PIETOOLS - stability_2d_sop
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
% is lpi_programming_sopvar/tests/e2e_stability_2d.m.
% MMP, 10/21/2026: lpi_ineq_sop reads P from the program (poscopvar records
%                it), so the options are gone. Deleted (same date):
%                "ineq_opts.like = P;" and the ineq_opts argument of
%                lpi_ineq_sop, with the help text on the option and on the
%                12036880 decision variables without it. Kept out of the
%                body as (was) lines, which echo would print.
% MMP, 10/21/2026: No solver options: lpisolve without them takes the first
%                solver installed, MOSEK and then SeDuMi (sossolve), which
%                the deleted lines chose by hand; options belong to the
%                settings files (maintainer). Deleted (same date): "%
%                MOSEK when it is installed, SeDuMi otherwise", "sos_opts =
%                struct('solver','sedumi');", "if ~isempty(which(
%                'mosekopt')), sos_opts.solver = 'mosek'; end" and the
%                sos_opts argument of lpisolve.
% MMP, 10/28/2026: Header to the header protocol (what the script
%                computes, inputs, outputs), the licence block added,
%                step comments in the body; code unchanged. The deleted
%                header text was the reference to stability_1d_sop. The
%                one-call executive is named PIETOOLS_stability_sop with
%                form 'P', which replaced PIETOOLS_stability_2D_sop on
%                10/28/2026. The help text of the 10/21/2026 entries is in
%                the paragraph after the LPI, without markers.

echo on
% The PDE and its PIE
pvar s1 s2 t
lam = pi^2;
x = pde_var('state',1,[s1;s2],[0,1;0,1]);
PDE = [diff(x,t) == diff(x,s1,2) + diff(x,s2,2) + lam*x;
       subs(x,s1,0) == 0;   subs(x,s1,1) == 0;
       subs(x,s2,0) == 0;   subs(x,s2,1) == 0];
PDE = initialize(PDE);
PIE = convert(PDE);
T = op2copvar_sop(PIE.T);                   % the PIE operators as containers
A = op2copvar_sop(PIE.A);                   % on L2[s1,s2]

% The LPI
dom = [0,1;0,1];
prog = lpiprogram_sop(T);
[prog,P] = poscopvar(prog,1,{{'s1','s2'}},dom,1);     % P >= 0 on L2[s1,s2]
P = P + 1e-3*eye_copvar_sop(1,{{'s1','s2'}},dom);     % P >= 1e-3*I
prog = lpi_ineq_sop(prog,-(T'*P*A + A'*P*T));        % T'*P*A + A'*P*T <= 0 % MMP, 10/21/2026
prog = lpisolve(prog);                      % MOSEK, else SeDuMi            % MMP, 10/21/2026
echo off

% Report the verdict of the solver
info = prog.solinfo.info;
if info.pinf==0 && info.dinf==0 && info.numerr==0
    fprintf('\n The PIE is stable: the LPI is feasible.\n');
else
    fprintf('\n The LPI is not solved (pinf %d, dinf %d, numerr %d).\n',info.pinf,info.dinf,info.numerr);
end
