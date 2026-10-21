%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% STABILITY_1D_SOP  Stability of a 1-D PDE by an LPI on the container path.
%
%       x_t = x_ss + lam*x,   s in [0,1],   x(0) = x(1) = 0,   lam = pi^2/2.
%
% The PDE is stable since lam < pi^2. Its PIE is T v_t = A v, and the PIE
% is stable if some P = P' >= eps*I satisfies
%
%       T'*P*A + A'*P*T <= 0,
%
% the LPI of PIETOOLS_stability (executives/). It is written here with the
% container user functions: the PIE operators become containers with
% op2copvar_sop, poscopvar declares P >= 0 on L2[s], and lpi_ineq_sop
% imposes the inequality, sizing the positive operator that certifies it.
% The inequality is linear in the decision operator P, so the option
% 'like' sizes that operator from P: sized from T'*P*A + A'*P*T itself, it
% covers every P of degree 2 (get_lift_degs), at 9225 decision variables
% instead of 2041. The executive PIETOOLS_stability_sop solves the same
% LPI in one call.
%
% For support, contact M. Peet, Arizona State University at mpeet@asu.edu
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Initial coding MMP, 10/21/2026. The end-to-end check that had this name
% is lpi_programming_sopvar/tests/e2e_stability_1d.m.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

echo on
% The PDE and its PIE
pvar s t
lam = pi^2/2;
x = pde_var('state',1,s,[0,1]);
PDE = [diff(x,t) == diff(x,s,2) + lam*x;
       subs(x,s,0) == 0;
       subs(x,s,1) == 0];
PDE = initialize(PDE);
PIE = convert(PDE);
T = op2copvar_sop(PIE.T);                   % the PIE operators as containers
A = op2copvar_sop(PIE.A);                   % on L2[s]

% MOSEK when it is installed, SeDuMi otherwise
sos_opts = struct('solver','sedumi');
if ~isempty(which('mosekopt')),     sos_opts.solver = 'mosek';  end

% The LPI
prog = lpiprogram_sop(T);
[prog,P] = poscopvar(prog,1,{{'s'}},[0,1],2);       % P >= 0 on L2[s], degree 2
P = P + 1e-3*eye_copvar_sop(1,{{'s'}},[0,1]);       % P >= 1e-3*I
ineq_opts.like = P;                         % size the certificate from P
prog = lpi_ineq_sop(prog,-(T'*P*A + A'*P*T),ineq_opts);    % T'*P*A + A'*P*T <= 0
prog = lpisolve(prog,sos_opts);

Psol = lpigetsol_sop(prog,P);               % the Lyapunov operator
echo off

info = prog.solinfo.info;
if info.pinf==0 && info.dinf==0 && info.numerr==0
    fprintf('\n The PIE is stable: the LPI is feasible.\n');
else
    fprintf('\n The LPI is not solved (pinf %d, dinf %d, numerr %d).\n',info.pinf,info.dinf,info.numerr);
end
