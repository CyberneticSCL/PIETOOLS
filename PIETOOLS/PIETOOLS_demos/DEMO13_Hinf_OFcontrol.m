% DEMO13_Hinf_OFcontrol.m
%
% This document illustrates how PIETOOLS can be used to synthesize an
% H-infinity optimal dynamic output feedback controller for a PDE.
% Specifically, we consider the reaction-diffusion equation, with a
% disturbance w1 and a control input u acting on the whole domain,
% PDE:       d/dt v(t,s) = 12*v(t,s) + (d^2/ds^2) v(t,s) + w1(t) + u(t),   s in [0,1];
% BCs:                0 = v(t,0) = v(t,1);
% OUTPUTS:     z(t) = [int_0^1 v(t,s) ds ; u(t)],
%              y(t) = int_0^1 v(t,s) ds + w2(t),
% where w2 is sensor noise. The exogenous input is w = [w1; w2].
% The PDE is unstable (the first eigenvalue of v_t = 12v + v_ss is 12-pi^2).
% We convert the PDE to a PIE of the form
%   T \dot{x}(t) = A x(t) + B1 w(t) + B2 u(t),
%         z(t) = C1 x(t) + D11 w(t) + D12 u(t),
%         y(t) = C2 x(t) + D21 w(t),
% and look for a dynamic output feedback controller, also a PIE,
%   TK \dot{xK}(t) = AK xK(t) + BK y(t),
%            u(t) = CK xK(t) + DK y(t),
% that minimizes a bound gam on the L2-gain from w to z of the closed loop.
% The synthesis LPI is declared step by step: positive operators R, V with
% the equalities T'*Q = R, T'*T*U = V, the free operators Q, U, An, Bn, Cn,
% Dn, the KYP inequality in these variables and the coupling
% [R T'*T; T'*T V] >= 0; the slack degree is raised until the solver
% certifies a solution or the raise count of the settings is reached.
% The controller operators are recovered from the solution by a change of
% variables, with no operator inverse.
%    Then, we
%  - bound the H-infinity norm of the closed loop with the analysis
%    executive PIETOOLS_Hinf_gain, applied to the PIE of the closed loop,
%    and report the bound only when the solver certifies it and it is not
%    below the numerical gain of the closed loop,
% -  and simulate the closed-loop dynamics using PIESIM.
%  - The commented block at the end compares gam with the finite-dimensional
%    H-infinity optimal controllers (hinfsyn, Robust Control Toolbox) for
%    finite-difference discretizations of the PDE.
% The script prints gam, warns if the synthesis LPI is not certified, prints
% the closed-loop bound or the reason it is not certified, and plots the
% simulation. Workspace: PIE, the controller Kc (fields TK, AK, BK, CK, DK),
% the solved variables sol, gam, the closed-loop PIE PIEcl, gam_cl, info_cl
% and solution_cl.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PIETOOLS - DEMO13
%
% Copyright (C)2026  PIETOOLS Team
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
% For support, contact M. Peet, Arizona State University at mpeet@asu.edu
% Initial coding: DB 10/31/2026, based on formulation of Jungbae Chun.

% Clear the workspace
clear; clc; clear stateNameGenerator

%% Declare the PDE and convert to a PIE
pvar t s
v  = pde_var(s,[0,1]);
w1 = pde_var('in');         % disturbance
w2 = pde_var('in');         % sensor noise
u  = pde_var('control');
z  = pde_var('out',2);      % regulated output [average; u]
y  = pde_var('sense');      % observed output
PDE = [diff(v,t)==12*v+diff(v,s,2)+w1+u;
       subs(v,s,0)==0;
       subs(v,s,1)==0;
       z==[int(v,s,[0,1]); u];
       y==int(v,s,[0,1])+w2];
PDE = initialize(PDE);
PIE = convert(PDE,'pie');
T = PIE.T;  A = PIE.A;  Bw = PIE.B1;  Bu = PIE.B2;
Cz = PIE.C1;  Dzw = PIE.D11;  Dzu = PIE.D12;  Cy = PIE.C2;  Dyw = PIE.D21;
dom = PIE.dom;
%% Synthesize the controller
% The positivity margins eppos of the operators R and V of the LPI are taken
% small: the free operators Q and U rescale R and V.
settings = lpisettings('light');
% margin on R, V: eppos on R^n (entry 1) and on L2[s] (entry 2)
ep = [1e-4; 1e-6];
%% Identity 'opvar' on the output spaces of an operator.
% Used for the -gam*I blocks of the KYP operator, on the disturbance space w
% (the output space of Bw') and on the output space z (of Cz).
Iw = eye_opvar(Bw');        % identity on w
Iz = eye_opvar(Cz);         % identity on z
% ep*I on the state (weights per block: R^n -> ep(1), L2 -> ep(2)) used to
% impose coercivity of the Lyapunov variables; for each block of the state
% (the input of T), the variables spi{k} ({} = R^n, {s} = L2[s]) and the
% dimension dmi(k).
[spi,dmi] = opvar_space_list(T,'in');
w = cellfun(@(c) ep(1+(~isempty(c))),spi);
blk = arrayfun(@(k) w(k)*eye(dmi(k)),1:numel(spi),'UniformOutput',false);
% mat2opvar(M,DM,SP,DOM): turns the numeric matrix M (block structure DM)
%   into a constant 'opvar' on the spaces SP over DOM (blocks: multiplication by
%   the matrix on that space). Here M = blkdiag(ep*I) is the margin operator
%   added to R and V so that they are >= ep*I, not just >= 0, in the raise loop.
Ieps = mat2opvar(blkdiag(blk{:}),dmi,spi,dom);

%% raise the slack degree until a feasible solution is found
% The raise count and the weight raise of the settings hold one entry per
% number of spatial variables: take the entry for the dimension of the PIE
nvp = PIE.dim;
max_raise = settings.sop.max_raise;  max_raise = max_raise(min(nvp,end));
dw_Q      = settings.sop.dw_Q;          dw_Q = dw_Q(min(nvp,end));
% LPI program on the variables and domain of the PIE, with the settings;
% the field sopeq records the rows of each equality constraint
prog = lpiprogram(PIE.vars(:,1),PIE.vars(:,2),dom,settings);
prog.sopeq = {};
[prog,gamv] = lpidecvar(prog,'gam');           % gam >= 0, minimized
prog = lpi_ineq(prog,gamv);
prog = lpisetobj(prog,gamv);
% Lyapunov variables R, V >= eppos I on the state (the output space of
% T'), with the degrees of the storage row of the settings (POSLPIVAR
% with DEG omitted)
[prog,R] = poslpivar(prog,T');
R = R + Ieps;
[prog,V] = poslpivar(prog,T');
V = V + Ieps;
% get_lpivar_degs(R): reads the monomial degrees actually used in the
%   positive operator R (multiplier, integral kernels), returned in the
%   format lpivar expects, so the free variables below get the
%   same degrees as R (dQ from R for Q,Bn,Dn; dU from V for U,An,Cn).
dQ = get_lpivar_degs(R);
dU = get_lpivar_degs(V);
%% indefinite operator variables Q,U,An,Bn,Cn,Dn.
% Lyapunov variables:  Q_cl = [Q S; S S], for functional <x,Rx>==<x,T*Qx>.
%  Define U := (Q − S)^{-1}, such that V = T U. The congruence transform
%  U_cl = [I U; 0 −U] linearizes the closed-loop inequality.
% Each free operator takes the output and input spaces of an operator of
% the same shape: Q those of T (state -> state), U and An those of R
% (state -> state), Bn those of Cy' (y -> state), Cn those of Bu'
% (state -> u), Dn those of (Cy*Bu)' (y -> u).
[prog,Q]  = lpivar(prog,T,dQ);
[prog,U]  = lpivar(prog,R,dU);
% controller variables (AK,BK,CK,DK) defined as:
%  (assumes D22 = 0)
%   An := Q'*(A + Bu*DK*Cy - Bu*CK)*U + S'*(BK*Cy - AK)*U
%   Bn := Q'*Bu*DK + S'*BK
%   Cn := (DK*Cy - CK)*U
%   Dn := DK
% (AK,BK,CK,DK are recovered from the solution at the end).
[prog,An] = lpivar(prog,R,dU);
[prog,Bn] = lpivar(prog,Cy',dQ);
[prog,Cn] = lpivar(prog,Bu',dU);
[prog,Dn] = lpivar(prog,(Cy*Bu)',dQ);
%% declare the operator equality and inequality constraints.
Tc = T'*T;
prog = lpi_eq(prog,T'*Q - R);                  % T'Q = R
prog = lpi_eq(prog,Tc*U - V);                  % T'T U = V
% KYP operator blocks
P11 = Q'*A + Bn*Cy;           P11 = P11 + P11';
P12 = An + (A + Bu*Dn*Cy)'*T;
P13 = Q'*Bw + Bn*Dyw;
P14 = Cz' + Cy'*Dn'*Dzu';
P22 = T'*(A*U + Bu*Cn);       P22 = P22 + P22';
P23 = T'*(Bw + Bu*Dn*Dyw);
P24 = U'*Cz' + Cn'*Dzu';
P34 = Dzw' + Dyw'*Dn'*Dzu';
Km = [P11,  P12,   P13,        P14;
      P12', P22,   P23,        P24;
      P13', P23',  -(gamv*Iw), P34;
      P14', P24',  P34',       -(gamv*Iz)];
Cp = [R, Tc;
         Tc', V];                          % coupling [R T'T; T'T V] >= 0
% Slack options: degrees like R, weight raise dw_Q + dw, margin of the settings
o = struct('like',R,'dw',dw_Q+max_raise,'margin',settings.sop.margin,'cap',[]);
prog = lpi_ineq(prog,-Km,o);                   % Km <= 0
prog = lpi_ineq(prog,Cp,o);                    % Cp >= 0
% Solve with the solver options of the settings
prog = lpisolve(prog,settings.sos_opts);

% Solver status: primal and dual infeasibility, numerical error
si = prog.solinfo.info;   pinf = 0;  dinf = 0;  nerr = 0;      % missing flag = 0
if isfield(si,'pinf')
    pinf = double(si.pinf);
end      % 1: primal infeasible
if isfield(si,'dinf')
    dinf = double(si.dinf);
end      % 1: dual infeasible
if isfield(si,'numerr')
    nerr = double(si.numerr);
end      % 1: inaccurate, 2: failed
gam = double(lpigetsol(prog,gamv));
ok = pinf~=1 && dinf~=1 && nerr<=1;            % optimal or inaccurate
% gam is a bound only when the solver certified it: warn otherwise; print gam
if ~ok
    warning('LPI not certified (pinf=%d dinf=%d numerr=%d).',pinf,dinf,nerr);
end
fprintf('Hinf_OFcontrol_demo13: gam = %.8g (certified: %d, dw = %d)\n',gam,ok,dw_Q+max_raise);
%% controller structure with the solved variables   
% getsol_lpivar(prog,X) needs the operator variable X itself, not its name,
% so the variables are listed next to their names.
sol_vars  = {Q, U, An, Bn, Cn, Dn, R, V};
sol_names = {'Q','U','An','Bn','Cn','Dn','R','V'};
for k = 1:numel(sol_vars)
    sol.(sol_names{k}) = getsol_lpivar(prog,sol_vars{k});
end
% The controller from the change of variables, with no inverse
DK = sol.Dn;
CK = DK*Cy*sol.U - sol.Cn;
BK = sol.Bn - sol.Q'*Bu*DK;
AK = BK*Cy*sol.U + sol.Q'*(A*sol.U + Bu*sol.Cn) - sol.An;
TK = T'*(sol.Q*sol.U - T);
% Controller operators TK, AK, BK, CK, DK ('opvar' objects)
Kc = struct('TK',TK,'AK',AK,'BK',BK,'CK',CK,'DK',DK);


%% Certify the H-infinity norm of the closed loop
% The closed loop of the plant and the controller is a PIE with state [x; xK].
% Its L2-gain is bounded with the KYP LPI of PIETOOLS_Hinf_gain, in the form
% 'Q', for which the storage functional is <Qcl x, Tcl x>. The controller
% operators carry factors of T, so the storage is declared with higher
% degrees than the default.
PIEcl = closedLoopPIE(PIE,Kc,'dynamic');
settings_an = lpisettings('light');
% storage degrees in one variable (multiplier 4; integrals 6, 6 in total),
% as lengths of the storage row, entry 1 (one variable): d1 = d2 = 7,
% d3 = [5;7;7], theta + s <= 6
settings_an.lpi.pos.len_inc{1} = 7;     settings_an.lpi.pos.len_com{1} = 7;
settings_an.lpi.pos.len_sos{1} = [5;7;7];   settings_an.lpi.pos.cap{1} = 6;
settings_an.eppos = 1e-6;   settings_an.eppos2 = 1e-6;
% The witness (on by default) is the gain of the discretized closed loop, a
% lower bound on it: a value below it is no bound, whatever the status
[~,~,gam_cl,info_cl] = PIETOOLS_Hinf_gain(PIEcl,settings_an,'Q');
cert_cl = info_cl.certified && ~(info_cl.gap<0);
if cert_cl
    fprintf('The H-infinity norm of the closed loop is bounded by %.4f.\n',gam_cl);
else
    fprintf(['The analysis LPI certifies no bound on the closed loop (status %s, value %.4f;\n' ...
             'numerical closed-loop gain %.4f).\n'],info_cl.status,gam_cl,info_cl.gain_lb);
end

%% Simulate the closed-loop PIE with PIESIM
% The closed-loop PIE PIEcl has state [x; xK] (plant PIE state, controller state),
% inputs w = [w1; w2] (disturbance, sensor noise) and output z.
syms st sx real
% initial conditions, one entry per scalar state component, here [x; xK]:
% plant PIE state x = v_ss starts at sin(pi*sx), controller starts at rest
uinput.ic = [sin(pi*sx); 0*sx];
% disturbances, functions of time st: w1 pulse-like, w2 small sensor noise
uinput.w(1) = exp(-st)*sin(2*pi*st);    % w1
uinput.w(2) = 0.1*sin(5*st);            % w2
uinput.ifexact = false;
% simulation options: Chebyshev order N, final time tf, BDF scheme of order Norder
opts_sim.N = 8;         opts_sim.tf = 2;        opts_sim.dt = 0.003;
opts_sim.intScheme = 1; opts_sim.Norder = 2;
opts_sim.plot = 'yes';  opts_sim.ploteig = 'no';    opts_sim.ifexact = false;
[solution_cl,grid_cl] = PIESIM(PIEcl,opts_sim,uinput);

%% Compare with the finite-dimensional baseline 
% If Robust Control Toolbox available, discretize the PDE by finite differences on N+1 grid points, and synthesize
% an H-infinity optimal output feedback controller with hinfsyn.
% if exist('hinfsyn','file')
%     Nlist = [10 20 50 100];     gam_N = zeros(size(Nlist));
%     for k = 1:numel(Nlist)
%         N = Nlist(k);   h = 1/N;
%         LN = -2*eye(N-1) + diag(ones(N-2,1),1) + diag(ones(N-2,1),-1);
%         AN = 12*eye(N-1) + LN/h^2;
%         B1N = [ones(N-1,1), zeros(N-1,1)];      B2N = ones(N-1,1);
%         C1N = [h*ones(1,N-1); zeros(1,N-1)];    C2N = h*ones(1,N-1);
%         D11N = zeros(2,2);  D12N = [0; 1];      D21N = [0 1];   D22N = 0;
%         PN = ss(AN,[B1N B2N],[C1N; C2N],[D11N D12N; D21N D22N]);
%         % RelTol 1e-7: the default 1e-2 returns gamma up to 1% above the optimum
%         [~,~,gam_N(k)] = hinfsyn(PN,1,1,hinfsynOptions('RelTol',1e-7));
%     end
%     fprintf('\nBaseline: finite-difference discretization and hinfsyn\n');
%     for k = 1:numel(Nlist),     fprintf('  N = %3d:  gamma = %.4f\n',Nlist(k),gam_N(k));     end
%     fprintf('PIE:  gam (LPI) = %.4f\n',gam);
%     if cert_cl,     fprintf('      closed-loop norm bound = %.4f\n',gam_cl);    end
% end
