%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% VOLTERRA_NORM_SOP computes an upper bound on the norm of the Volterra
% operator
%
%       (T x)(s) = int_0^s x(r) dr,   s in [0,1],   ||T|| = 2/pi,
%
% from the LPI
%
%       min gam   s.t.   gam*I - T'*T >= 0,
%
% on the container path: T is declared as an 'opvar' object and converted
% to a 'copvar' object by OP2COPVAR_SOP, and LPI_INEQ_SOP imposes the
% inequality with a certificate sized from gam*I - T'*T, so no degree is
% given, printing its decision variables first. The script prints the
% bound sqrt(gam) on ||T|| and 2/pi.
%
% INPUT
% none.
%
% OUTPUT
% the workspace variables
% - Top:            'opvar' object, the Volterra operator;
% - T:              'copvar' object, Top as a container;
% - prog:           'struct', the solved LPI program;
% - gam:            'dpvar' object, the decision variable gam;
% - operator_norm:  double, sqrt(gam) at the solution, the bound on ||T||.
%
% For support, contact M. Peet, Arizona State University at mpeet@asu.edu

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PIETOOLS - volterra_norm_sop
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
% Initial coding MMP, 09/29/2026
% MMP, 10/06/2026: cx_space_list -> library copvar_space_list.
% MMP, 10/21/2026: Rewritten to use the user functions only. Deleted: the
%                  legacy DEMO2 solve run alongside, the positive operator
%                  declared by hand (two poscopvar terms at degree 1) with
%                  gam - T'*T - W = 0 imposed by lpi_eq_cdopvar, and the
%                  residual check by the test helper cx_resid (lines of
%                  09/29 and 10/06/2026). lpi_ineq_sop now does what those
%                  lines did, and the example is a script, as the DEMOs are.
% MMP, 10/21/2026: No solver options: lpisolve without them takes the first
%                  solver installed, MOSEK and then SeDuMi (sossolve), which
%                  the deleted lines chose by hand; options belong to the
%                  settings files (maintainer). Deleted (same date): "%
%                  MOSEK when it is installed, SeDuMi otherwise", "sos_opts
%                  = struct('solver','sedumi');", "if ~isempty(which(
%                  'mosekopt')), sos_opts.solver = 'mosek'; end" and the
%                  sos_opts argument of lpisolve.
% MMP, 10/28/2026: Header to the header protocol (what the script
%                  computes, inputs, outputs), the licence block added,
%                  step comments in the body; code unchanged. The deleted
%                  header text was the comparison with DEMO2 and the list
%                  of the user and stock functions the script calls.

% % % MMP, 10/21/2026: start of the rewritten body (header note).
echo on
% The operator, as in DEMO2, and its container form on L2[0,1]
opvar Top;
Top.R.R1 = 1;   Top.I = [0,1];
T = op2copvar_sop(Top);

% The LPI
prog = lpiprogram_sop(T);                   % a program on the variables of T
[prog,gam] = lpidecvar(prog,'gam');         % a scalar decision variable
prog = lpi_ineq_sop(prog,gam - T'*T);       % gam*I - T'*T >= 0
prog = lpisetobj(prog,gam);                 % minimize gam
prog = lpisolve(prog);                      % MOSEK, else SeDuMi            % MMP, 10/21/2026

% The bound on ||T||
operator_norm = sqrt(double(lpigetsol_sop(prog,gam)));
echo off

% Print the bound and the exact norm
fprintf('\n Upper bound on the norm of the Volterra operator: %.6f (exact: 2/pi = %.6f).\n',...
        operator_norm,2/pi);
% % % MMP, 10/21/2026: end of the rewritten body.
