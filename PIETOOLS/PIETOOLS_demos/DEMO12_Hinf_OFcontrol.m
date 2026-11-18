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
% The controller operators are recovered from the solution of the LPI in PIETOOLS_Hinf_OFcontrol_sop.
%    Then, we
%  - bound the H-infinity norm of the closed loop with the analysis
%    executive PIETOOLS_Hinf_gain_sop, applied to the PIE of the closed loop,
%    and report the bound only when the solver certifies it and it is not
%    below the numerical gain of the closed loop,
%  - compare gam with the finite-dimensional H-infinity optimal controllers
%    (hinfsyn) for finite-difference discretizations of the PDE.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PIETOOLS - DEMO12
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
% Initial coding: DB 10/31/2026, based on formulation of Jungbae Chun.
% MMP 11/05/2026: A bound is printed only when certified: gam when
%                 INFO.certified; the closed-loop bound when INFO.certified
%                 and not below the numerical gain (the witness, no longer
%                 switched off). The closed-loop LPI stalls (relres 1.6e-5),
%                 MOSEK reports the stalled point as optimal, and its value
%                 5.6129 was below the closed-loop gain 5.6150. hinfsyn at
%                 RelTol 1e-7: the default 1e-2 returned gamma up to 1% above
%                 the optimum, which put it above the certified PIE bound.
% MMP 11/08/2026: Settings lpisettings_sop('light') for the synthesis and
%                 the closed-loop analysis (was the stock lpisettings('light'),
%                 whose folder leaves the path).
% MMP 11/12/2026: The storage degrees of the closed-loop analysis in the
%                 native settings, lpi.pos (maintainer, 11/12/2026: phase 2
%                 of the knob reduction): mult 4, int 6, joint 6, the former
%                 dd1 = {4,[6 6 6],[6 6 6]}, which the settings no longer
%                 accept. Same program.
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
% settings = lpisettings('light');                                          % MMP, 11/08/2026 (was)
settings = lpisettings_sop('light');                                        % MMP, 11/08/2026
% settings.eppos = 1e-7;      settings.eppos2 = 1e-7;         % DB, 11/17/2026 (was)
% margin on R, V: eppos on R^n (entry 1) and on L2[s] (entry 2) % DB, 11/17/2026 
ep = [1e-4; 1e-6];                    % DB, 11/17/2026 
%% For each 'copvar' P, get for each block row ('out') or column ('in')   % DB, 11/17/2026
% the variables sp{k} ({} = R^n,  {s} = L2[s]) and its dimensions dm(k).           % DB, 11/17/2026
% These lists describe the spaces of the PIE (state x in R^n x L2[s]; inputs w,u; outputs z,y) and      % DB, 11/17/2026
% are the size/space arguments lpivar_cdopvar needs for Q,U,An,Bn,Cn,Dn. % DB, 11/17/2026
[spi,dmi] = copvar_space_list(T,'in');                                    % DB, 11/17/2026 
[spo,dmo] = copvar_space_list(T,'out');                                % DB, 11/17/2026
[spu,dmu] = copvar_space_list(Bu,'in');                                 % DB, 11/17/2026
[spy,dmy] = copvar_space_list(Cy,'out');                                % DB, 11/17/2026
%% Identity 'copvar' on the spaces (SP,DM) over the domain DOM.                 % DB, 11/17/2026
% Used for the -gam*I blocks of the KYP operator, on the disturbance space w        % DB, 11/17/2026
% (columns of B1) and on the output space z (rows of C1).                                           % DB, 11/17/2026
[s,d] = copvar_space_list(Bw,'in');   Iw = eye_copvar_sop(d,s,dom);   % identity on w       % DB, 11/17/2026
[s,d] = copvar_space_list(Cz,'out');  Iz = eye_copvar_sop(d,s,dom);   % identity on z           % DB, 11/17/2026
% ep*I on the state (weights per block: R^n -> ep(1), L2 -> ep(2)) used to                      % DB, 11/17/2026
% impose coercivity of the Lyapunov variables.                                                                % DB, 11/17/2026
w = cellfun(@(c) ep(1+(~isempty(c))),spi);                                                                         % DB, 11/17/2026
blk = arrayfun(@(k) w(k)*eye(dmi(k)),1:numel(spi),'UniformOutput',false);                         % DB, 11/17/2026
% mat2copvar_sop(M,DM,SP,DOM): turns the numeric matrix M (block structure DM)      % DB, 11/17/2026
%   into a constant 'copvar' on the spaces SP over DOM (blocks: multiplication by           % DB, 11/17/2026
%   the matrix on that space). Here M = blkdiag(ep*I) is the margin operator                      % DB, 11/17/2026
%   added to R and V so that they are >= ep*I, not just >= 0, at lines 138 and 141.                                               % DB, 11/17/2026
Ieps = mat2copvar_sop(blkdiag(blk{:}),dmi,spi,dom);                                                             % DB, 11/17/2026

% % [prog,Kc,gam] = PIETOOLS_Hinf_OFcontrol_sop(PIE,settings);                % MMP, 11/05/2026 (was)
% [prog,Kc,gam,info] = PIETOOLS_Hinf_OFcontrol_sop(PIE,settings);
% % MMP, 11/05/2026  % DB 11/17/2026(was)
% size/space arguments of lpivar_cdopvar for each free operator (rows 'out' x columns 'in'):   % DB 11/17/2026
%   dimension struct with fields out/in, and space struct (variables of each block)   % DB 11/17/2026
dimQ  = struct('out',dmo,'in',dmi);   spQ  = struct('out',{spo},'in',{spi});   % Q  : state <- state   % DB 11/17/2026
dimU  = struct('out',dmi,'in',dmi);   spU  = struct('out',{spi},'in',{spi});   % U, An : state -> state   % DB 11/17/2026
dimBn = struct('out',dmi,'in',dmy);   spBn = struct('out',{spi},'in',{spy});   % Bn : y -> state   % DB 11/17/2026
dimCn = struct('out',dmu,'in',dmi);   spCn = struct('out',{spu},'in',{spi});   % Cn : state -> u   % DB 11/17/2026
dimDn = struct('out',dmu,'in',dmy);   spDn = struct('out',{spu},'in',{spy});   % Dn : y -> u   % DB 11/17/2026

%% raise the slack degree until a feasible solution is found   % DB 11/17/2026
for dw = 0:settings.sop.max_raise   % DB 11/17/2026
    prog = lpiprogram(PIE.vars(:,1),PIE.vars(:,2),dom,settings);   % DB 11/17/2026
    prog.sopeq = {};   % DB 11/17/2026
    [prog,gamv] = lpidecvar(prog,'gam');           % gam >= 0, minimized   % DB 11/17/2026
    prog = lpi_ineq_sop(prog,gamv);   % DB 11/17/2026
    prog = lpisetobj(prog,gamv);   % DB 11/17/2026
    % Lyapunov variables R, V >= eppos I   % DB 11/17/2026
    [prog,R] = poslpivar_settings_sop(prog,T,settings,'lf','in',dom);   % DB 11/17/2026
    R = R + Ieps;   % DB 11/17/2026
    [prog,V] = poslpivar_settings_sop(prog,T,settings,'lf','in',dom);  %  DB 11/17/2026 
    V = V + Ieps;         % DB 11/17/2026 
    % get_lpivar_degs_sop(R): reads the monomial degrees actually used in the   % DB 11/17/2026
    %   positive operator R (multiplier, integral kernels), returned in the   % DB 11/17/2026
    %   format lpivar_cdopvar expects, so the free variables below get the   % DB 11/17/2026
    %   same degrees as R (dQ from R for Q,Bn,Dn; dU from V for U,An,Cn).   % DB 11/17/2026
    dQ = get_lpivar_degs_sop(R);  % DB 11/17/2026
    dU = get_lpivar_degs_sop(V);   % DB 11/17/2026
    %% indefinite operator variables Q,U,An,Bn,Cn,Dn.   % DB 11/17/2026
    % Lyapunov variables:  Q_cl = [Q S   , for functional <x,Rx>==<x,T*Qx>.   % DB 11/17/2026
    %                                               S S]   % DB 11/17/2026
    %  Define U := (Q − S)^{-1}, such that V = T U. The congruence transform U_cl = [I U   % DB 11/17/2026
    %                                                                                                 0 −U]   % DB 11/17/2026
    % linearizes the closed-loop inequality.   % DB 11/17/2026
    [prog,Q]  = lpivar_cdopvar(prog,dimQ, spQ, dom,dQ);   % DB 11/17/2026
    [prog,U]  = lpivar_cdopvar(prog,dimU, spU, dom,dU);   % DB 11/17/2026
    % controller variables (AK,BK,CK,DK) defined as:   % DB 11/17/2026
    %  (assumes D22 = 0)   % DB 11/17/2026
    %   An := Q'*(A + Bu*DK*Cy - Bu*CK)*U + S'*(BK*Cy - AK)*U   % DB 11/17/2026
    %   Bn := Q'*Bu*DK + S'*BK   % DB 11/17/2026
    %   Cn := (DK*Cy - CK)*U   % DB 11/17/2026
    %   Dn := DK   % DB 11/17/2026
    % (AK,BK,CK,DK are recovered from the solution at the end).   % DB 11/17/2026
    [prog,An] = lpivar_cdopvar(prog,dimU, spU, dom,dU);   % DB 11/17/2026
    [prog,Bn] = lpivar_cdopvar(prog,dimBn,spBn,dom,dQ);   % DB 11/17/2026
    [prog,Cn] = lpivar_cdopvar(prog,dimCn,spCn,dom,dU);   % DB 11/17/2026
    [prog,Dn] = lpivar_cdopvar(prog,dimDn,spDn,dom,dQ);   % DB 11/17/2026
    %% declare the operator equality and inequality constraints.   % DB 11/17/2026
    Tc = T'*T;   % DB 11/17/2026
    prog = lpi_eq_sop(prog,T'*Q - R);              % T'Q = R   % DB 11/17/2026
    prog = lpi_eq_sop(prog,Tc*U - V);              % T'T U = V   % DB 11/17/2026
    % KYP operator blocks   % DB 11/17/2026
    P11 = Q'*A + Bn*Cy;           P11 = P11 + P11';   % DB 11/17/2026
    P12 = An + (A + Bu*Dn*Cy)'*T;   % DB 11/17/2026
    P13 = Q'*Bw + Bn*Dyw;   % DB 11/17/2026
    P14 = Cz' + Cy'*Dn'*Dzu';   % DB 11/17/2026
    P22 = T'*(A*U + Bu*Cn);       P22 = P22 + P22';   % DB 11/17/2026
    P23 = T'*(Bw + Bu*Dn*Dyw);   % DB 11/17/2026
    P24 = U'*Cz' + Cn'*Dzu';   % DB 11/17/2026
    P34 = Dzw' + Dyw'*Dn'*Dzu';   % DB 11/17/2026
    Km = [P11,  P12,   P13,        P14;   % DB 11/17/2026
          P12', P22,   P23,        P24;   % DB 11/17/2026
          P13', P23',  -(gamv*Iw), P34;   % DB 11/17/2026
          P14', P24',  P34',       -(gamv*Iz)];   % DB 11/17/2026
    Cp = [R, Tc;   % DB 11/17/2026
             Tc', V];                          % coupling [R T'T; T'T V] >= 0   % DB 11/17/2026
    o = struct('like',R,'dw',settings.sop.dw_Q+dw,'margin',settings.sop.margin,'cap',[]);   % DB 11/17/2026
    prog = lpi_ineq_sop(prog,-Km,o);               % Km <= 0   % DB 11/17/2026
    prog = lpi_ineq_sop(prog,Cp,o);                % Cp >= 0   % DB 11/17/2026
    prog = lpisolve(prog,settings.sos_opts);   % DB 11/17/2026

    % status (as lpi_classify_sop)   % DB 11/17/2026
    si = prog.solinfo.info;   pinf = 0;  dinf = 0;  nerr = 0;      % missing flag = 0   % DB 11/17/2026
    if isfield(si,'pinf')   % DB 11/17/2026
        pinf = double(si.pinf);   % DB 11/17/2026
    end      % 1: primal infeasible   % DB 11/17/2026
    if isfield(si,'dinf')   % DB 11/17/2026
        dinf = double(si.dinf);   % DB 11/17/2026
    end      % 1: dual infeasible   % DB 11/17/2026
    if isfield(si,'numerr')   % DB 11/17/2026
        nerr = double(si.numerr);   % DB 11/17/2026
    end      % 1: inaccurate, 2: failed   % DB 11/17/2026
    gam = double(lpigetsol_sop(prog,gamv));   % DB 11/17/2026
    ok = pinf~=1 && dinf~=1 && nerr<=1;            % optimal or inaccurate   % DB 11/17/2026
    if ok || dw==settings.sop.max_raise || ~(nerr==2 || pinf==1)   % DB 11/17/2026
        break   % DB 11/17/2026
    end   % DB 11/17/2026
end   % DB 11/17/2026
% gam is a bound only when the solver certified it
% % MMP, 11/05/2026      % DB 11/17/2026(was)
% if ~info.certified                                                          % MMP, 11/05/2026   % DB 11/17/2026(was)
%     error('The synthesis LPI certifies no bound (status %s).',info.status); % MMP, 11/05/2026  % DB 11/17/2026(was)
% end                                                                         % MMP, 11/05/2026      % DB 11/17/2026(was)
% fprintf('\nThe LPI certifies a closed-loop H-infinity bound gam = %.4f.\n',gam);  % DB 11/17/2026 (was)
if ~ok         % DB 11/17/2026
    warning('LPI not certified (pinf=%d dinf=%d numerr=%d).',pinf,dinf,nerr);     % DB 11/17/2026
end         % DB 11/17/2026
fprintf('Hinf_OFcontrol_demo12: gam = %.8g (certified: %d, dw = %d)\n',gam,ok,dw);     % DB 11/17/2026
%% controller structure with the solved variables   
% getsol_lpivar_sop(prog,X) needs the operator variable X itself, not its name,   % DB 11/17/2026
% so the variables are listed next to their names.   % DB 11/17/2026
sol_vars  = {Q, U, An, Bn, Cn, Dn, R, V};   % DB 11/17/2026
sol_names = {'Q','U','An','Bn','Cn','Dn','R','V'};   % DB 11/17/2026
for k = 1:numel(sol_vars)   % DB 11/17/2026
    sol.(sol_names{k}) = getsol_lpivar_sop(prog,sol_vars{k});   % DB 11/17/2026
end   % DB 11/17/2026
% The controller from the change of variables, with no inverse
DK = sol.Dn;
CK = DK*Cy*sol.U - sol.Cn;
BK = sol.Bn - sol.Q'*Bu*DK;
AK = BK*Cy*sol.U + sol.Q'*(A*sol.U + Bu*sol.Cn) - sol.An;   % DB 11/17/2026 (Am -> A)
TK = T'*(sol.Q*sol.U - T);   % DB 11/17/2026 (Tm -> T)
%   Kc = struct('TK',tl.toop(TK),'AK',tl.toop(AK),'BK',tl.toop(BK),'CK',tl.toop(CK),'DK',tl.toop(DK)); % MMP, 11/05/2026 (was)
% The containers: tl.toop needs the opvar classes in 1-D and 2-D        % MMP, 11/05/2026
Kc = struct('TK',TK,'AK',AK,'BK',BK,'CK',CK,'DK',DK);                   % MMP, 11/05/2026


%% Certify the H-infinity norm of the closed loop
% The closed loop of the plant and the controller is a PIE with state [x; xK].
% Its L2-gain is bounded with the KYP LPI of PIETOOLS_Hinf_gain_sop, in the form
% 'Q', for which the storage functional is <Qcl x, Tcl x>. The controller
% operators carry factors of T, so the storage is declared with higher
% degrees than the default.
PIEcl = closedLoopPIE(PIE,Kc,'dynamic');
% settings_an = lpisettings('light');                                       % MMP, 11/08/2026 (was)
settings_an = lpisettings_sop('light');                                     % MMP, 11/08/2026
% settings_an.dd1 = {4,[6 6 6],[6 6 6]};                                    % MMP, 11/12/2026 (was)
% storage degrees in one variable: multiplier 4; integrals 6, 6 in total    % MMP, 11/12/2026
settings_an.lpi.pos.mult(1) = 4;    settings_an.lpi.pos.int(1) = 6;         % MMP, 11/12/2026
settings_an.lpi.pos.joint(1) = 6;                                           % MMP, 11/12/2026
settings_an.eppos = 1e-6;   settings_an.eppos2 = 1e-6;
% settings_an.sop.witness = false;                                          % MMP, 11/05/2026 (was)
% [~,~,gam_cl] = PIETOOLS_Hinf_gain_sop(PIEcl,settings_an,'Q');             % MMP, 11/05/2026 (was)
% fprintf('The H-infinity norm of the closed loop is bounded by %.4f.\n',gam_cl); % MMP, 11/05/2026 (was)
% The witness (on by default) is the gain of the discretized closed loop, a % MMP, 11/05/2026
% lower bound on it: a value below it is no bound, whatever the status      % MMP, 11/05/2026
[~,~,gam_cl,info_cl] = PIETOOLS_Hinf_gain_sop(PIEcl,settings_an,'Q');       % MMP, 11/05/2026
cert_cl = info_cl.certified && ~(info_cl.gap<0);                            % MMP, 11/05/2026
if cert_cl                                                                  % MMP, 11/05/2026
    fprintf('The H-infinity norm of the closed loop is bounded by %.4f.\n',gam_cl); % MMP, 11/05/2026
else                                                                        % MMP, 11/05/2026
    fprintf(['The analysis LPI certifies no bound on the closed loop (status %s, value %.4f;\n' ...
             'numerical closed-loop gain %.4f).\n'],info_cl.status,gam_cl,info_cl.gain_lb); % MMP, 11/05/2026
end                                                                         % MMP, 11/05/2026

%% Simulate the closed-loop PIE with PIESIM   % DB 11/17/2026
% The closed-loop PIE PIEcl has state [x; xK] (plant PIE state, controller state),   % DB 11/17/2026
% inputs w = [w1; w2] (disturbance, sensor noise) and output z.   % DB 11/17/2026
syms st sx real   % DB 11/17/2026
% initial conditions, one entry per scalar state component, here [x; xK]:   % DB 11/17/2026
% plant PIE state x = v_ss starts at sin(pi*sx), controller starts at rest   % DB 11/17/2026
uinput.ic = [sin(pi*sx); 0*sx];   % DB 11/17/2026
% disturbances, functions of time st: w1 pulse-like, w2 small sensor noise   % DB 11/17/2026
uinput.w(1) = exp(-st)*sin(2*pi*st);    % w1   % DB 11/17/2026
uinput.w(2) = 0.1*sin(5*st);            % w2   % DB 11/17/2026
uinput.ifexact = false;   % DB 11/17/2026
% simulation options: Chebyshev order N, final time tf, BDF scheme of order Norder   % DB 11/17/2026
opts_sim.N = 8;         opts_sim.tf = 2;        opts_sim.dt = 0.003;   % DB 11/17/2026
opts_sim.intScheme = 1; opts_sim.Norder = 2;   % DB 11/17/2026
opts_sim.plot = 'yes';  opts_sim.ploteig = 'no';    opts_sim.ifexact = false;   % DB 11/17/2026
[solution_cl,grid_cl] = PIESIM(PIEcl,opts_sim,uinput);   % DB 11/17/2026

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
% %       [~,~,gam_N(k)] = hinfsyn(PN,1,1);                                   % MMP, 11/05/2026 (was)
%         % RelTol 1e-7: the default 1e-2 returns gamma up to 1% above the optimum % MMP, 11/05/2026
%         [~,~,gam_N(k)] = hinfsyn(PN,1,1,hinfsynOptions('RelTol',1e-7));     % MMP, 11/05/2026
%     end
%     fprintf('\nBaseline: finite-difference discretization and hinfsyn\n');
%     for k = 1:numel(Nlist),     fprintf('  N = %3d:  gamma = %.4f\n',Nlist(k),gam_N(k));     end
% %   fprintf('PIE:  gam (LPI) = %.4f, closed-loop norm bound = %.4f\n',gam,gam_cl); % MMP, 11/05/2026 (was)
%     fprintf('PIE:  gam (LPI) = %.4f\n',gam);                                % MMP, 11/05/2026
%     if cert_cl,     fprintf('      closed-loop norm bound = %.4f\n',gam_cl);    end % MMP, 11/05/2026
% end
