%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% HINF_GAIN_1D_SOP  The H-infinity gain of a 1-D PDE by an LPI on the
% container path, for the system
%
%       x_t = x_ss + (pi^2/2) x + s(1-s) w,   s in [0,1],   x(0) = x(1) = 0,
%       z   = int_0^1 x(s) ds.
%
% The PIE is T v_t = A v + B w, z = C v + D w, and gam bounds the gain
% from w to z if some P = P' >= 0 satisfies
%
%       [ -gam*I     D'       B'*P*T         ]
%       [  D        -gam*I    C              ]  <= 0,
%       [  T'*P*B    C'       T'*P*A + A'*P*T]
%
% the LPI of PIETOOLS_Hinf_gain with the storage P (executives/). It is
% written here with the container user functions; a scalar decision
% variable placed in a block matrix is gam times the identity of that
% block. The inequality is linear in the decision operator P, so the
% option 'like' sizes the positive operator that certifies it from P:
% sized from the block operator itself, it covers every P of degree 2
% (get_lift_degs), at 27426 decision variables instead of 2274, for the
% same gain. The executive PIETOOLS_Hinf_gain_sop solves the same problem
% in one call and is run at the end for comparison.
%
% For support, contact M. Peet, Arizona State University at mpeet@asu.edu
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Initial coding MMP, 10/21/2026. The end-to-end check that had this name
% is lpi_programming_sopvar/tests/e2e_hinf_gain_1d.m.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

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
T = op2copvar_sop(PIE.T);   A = op2copvar_sop(PIE.A);
B = op2copvar_sop(PIE.B1);  C = op2copvar_sop(PIE.C1);  D = op2copvar_sop(PIE.D11);

% MOSEK when it is installed, SeDuMi otherwise
sos_opts = struct('solver','sedumi');
if ~isempty(which('mosekopt')),     sos_opts.solver = 'mosek';  end

% The LPI
prog = lpiprogram_sop(T);
[prog,gam] = lpidecvar(prog,'gam');
[prog,P] = poscopvar(prog,1,{{'s'}},[0,1],2);       % P >= 0 on L2[s], degree 2
K = [-gam,       D',     B'*P*T;
      D,         -gam,   C;
      T'*P*B,    C',     T'*P*A + A'*P*T];
ineq_opts.like = P;                                 % size the certificate from P
prog = lpi_ineq_sop(prog,-K,ineq_opts);             % K <= 0
prog = lpisetobj(prog,gam);                         % minimize gam
prog = lpisolve(prog,sos_opts);

gain = double(lpigetsol_sop(prog,gam));

% The same problem through the executive
settings = lpisettings('light');
settings.sos_opts = sos_opts;
[~,~,gain_exec] = PIETOOLS_Hinf_gain_sop(PIE,settings);
echo off

fprintf('\n H-infinity gain bound: %.6f by the LPI above, %.6f by PIETOOLS_Hinf_gain_sop.\n',...
        gain,gain_exec);
