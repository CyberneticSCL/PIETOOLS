%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% STABILITY_2D_SOP  Stability of a 2-D PDE by an LPI on the container path.
%
%       x_t = x_{s1 s1} + x_{s2 s2} + lam*x  on [0,1]^2,  x = 0 on the
%       boundary,  lam = pi^2.
%
% The PDE is stable since lam < 2 pi^2. As in 1-D (stability_1d_sop), its
% PIE T v_t = A v is stable if some P = P' >= eps*I satisfies
%
%       T'*P*A + A'*P*T <= 0.
%
% The container functions are the same in two variables: poscopvar
% declares P >= 0 on L2[s1,s2], and lpi_ineq_sop sizes the positive
% operator that certifies the inequality, with the face terms of the box
% that two variables need. As in 1-D, the option 'like' sizes that
% operator from P: 197920 decision variables. Sized from T'*P*A + A'*P*T
% itself, it covers every P of degree 1 (get_lift_degs), at 12036880
% decision variables, unsolved after 900 s and 31 GB. The executive
% PIETOOLS_stability_2D_sop solves the same LPI in one call.
%
% For support, contact M. Peet, Arizona State University at mpeet@asu.edu
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Initial coding MMP, 10/21/2026. The end-to-end check that had this name
% is lpi_programming_sopvar/tests/e2e_stability_2d.m.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

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

% MOSEK when it is installed, SeDuMi otherwise
sos_opts = struct('solver','sedumi');
if ~isempty(which('mosekopt')),     sos_opts.solver = 'mosek';  end

% The LPI
dom = [0,1;0,1];
prog = lpiprogram_sop(T);
[prog,P] = poscopvar(prog,1,{{'s1','s2'}},dom,1);     % P >= 0 on L2[s1,s2]
P = P + 1e-3*eye_copvar_sop(1,{{'s1','s2'}},dom);     % P >= 1e-3*I
ineq_opts.like = P;                                   % size the certificate from P
prog = lpi_ineq_sop(prog,-(T'*P*A + A'*P*T),ineq_opts);    % T'*P*A + A'*P*T <= 0
prog = lpisolve(prog,sos_opts);
echo off

info = prog.solinfo.info;
if info.pinf==0 && info.dinf==0 && info.numerr==0
    fprintf('\n The PIE is stable: the LPI is feasible.\n');
else
    fprintf('\n The LPI is not solved (pinf %d, dinf %d, numerr %d).\n',info.pinf,info.dinf,info.numerr);
end
