%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% VOLTERRA_NORM_SOP  An upper bound on the norm of the Volterra operator,
%
%       (T x)(s) = int_0^s x(r) dr,   s in [0,1],   ||T|| = 2/pi,
%
% from the LPI
%
%       min gam   s.t.   gam*I - T'*T >= 0,
%
% on the container path: DEMO2 (PIETOOLS_demos/DEMO2_volterra_operator_
% norm.m) with the operator converted to a container by op2copvar_sop and
% the LPI written with the container user functions lpiprogram_sop,
% lpi_ineq_sop and lpigetsol_sop; lpidecvar, lpisetobj and lpisolve are the
% stock functions. lpi_ineq_sop sizes the positive operator that certifies
% gam*I - T'*T >= 0 from that operator itself, so no degree is given.
%
% For support, contact M. Peet, Arizona State University at mpeet@asu.edu
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Initial coding MMP, 09/29/2026
% MMP, 10/06/2026: cx_space_list -> library copvar_space_list.
% MMP, 10/21/2026: Rewritten to use the user functions only. Deleted: the
%                  legacy DEMO2 solve run alongside, the positive operator
%                  declared by hand (two poscopvar terms at degree 1) with
%                  gam - T'*T - W = 0 imposed by lpi_eq_cdopvar, and the
%                  residual check by the test helper cx_resid (lines of
%                  09/29 and 10/06/2026). lpi_ineq_sop now does what those
%                  lines did, and the example is a script, as the DEMOs are.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% % % MMP, 10/21/2026: start of the rewritten body (header note).
echo on
% The operator, as in DEMO2, and its container form on L2[0,1]
opvar Top;
Top.R.R1 = 1;   Top.I = [0,1];
T = op2copvar_sop(Top);

% MOSEK when it is installed, SeDuMi otherwise
sos_opts = struct('solver','sedumi');
if ~isempty(which('mosekopt')),     sos_opts.solver = 'mosek';  end

% The LPI
prog = lpiprogram_sop(T);                   % a program on the variables of T
[prog,gam] = lpidecvar(prog,'gam');         % a scalar decision variable
prog = lpi_ineq_sop(prog,gam - T'*T);       % gam*I - T'*T >= 0
prog = lpisetobj(prog,gam);                 % minimize gam
prog = lpisolve(prog,sos_opts);

operator_norm = sqrt(double(lpigetsol_sop(prog,gam)));
echo off

fprintf('\n Upper bound on the norm of the Volterra operator: %.6f (exact: 2/pi = %.6f).\n',...
        operator_norm,2/pi);
% % % MMP, 10/21/2026: end of the rewritten body.
