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

%% Synthesize the controller
% The positivity margins eppos of the operators R and V of the LPI are taken
% small: the free operators Q and U rescale R and V.
settings = lpisettings('light');
settings.eppos = 1e-7;      settings.eppos2 = 1e-7;
% [prog,Kc,gam] = PIETOOLS_Hinf_OFcontrol_sop(PIE,settings);                % MMP, 11/05/2026 (was)
[prog,Kc,gam,info] = PIETOOLS_Hinf_OFcontrol_sop(PIE,settings);             % MMP, 11/05/2026
% gam is a bound only when the solver certified it                          % MMP, 11/05/2026
if ~info.certified                                                          % MMP, 11/05/2026
    error('The synthesis LPI certifies no bound (status %s).',info.status); % MMP, 11/05/2026
end                                                                         % MMP, 11/05/2026
fprintf('\nThe LPI certifies a closed-loop H-infinity bound gam = %.4f.\n',gam);

%% Certify the H-infinity norm of the closed loop
% The closed loop of the plant and the controller is a PIE with state [x; xK].
% Its L2-gain is bounded with the KYP LPI of PIETOOLS_Hinf_gain_sop, in the form
% 'Q', for which the storage functional is <Qcl x, Tcl x>. The controller
% operators carry factors of T, so the storage is declared with higher
% degrees than the default.
PIEcl = closedLoopPIE(PIE,Kc,'dynamic');
settings_an = lpisettings('light');
settings_an.dd1 = {4,[6 6 6],[6 6 6]};
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

%% Compare with the finite-dimensional baseline 
% If Robust Control Toolbox available, discretize the PDE by finite differences on N+1 grid points, and synthesize
% an H-infinity optimal output feedback controller with hinfsyn.
if exist('hinfsyn','file')
    Nlist = [10 20 50 100];     gam_N = zeros(size(Nlist));
    for k = 1:numel(Nlist)
        N = Nlist(k);   h = 1/N;
        LN = -2*eye(N-1) + diag(ones(N-2,1),1) + diag(ones(N-2,1),-1);
        AN = 12*eye(N-1) + LN/h^2;
        B1N = [ones(N-1,1), zeros(N-1,1)];      B2N = ones(N-1,1);
        C1N = [h*ones(1,N-1); zeros(1,N-1)];    C2N = h*ones(1,N-1);
        D11N = zeros(2,2);  D12N = [0; 1];      D21N = [0 1];   D22N = 0;
        PN = ss(AN,[B1N B2N],[C1N; C2N],[D11N D12N; D21N D22N]);
%       [~,~,gam_N(k)] = hinfsyn(PN,1,1);                                   % MMP, 11/05/2026 (was)
        % RelTol 1e-7: the default 1e-2 returns gamma up to 1% above the optimum % MMP, 11/05/2026
        [~,~,gam_N(k)] = hinfsyn(PN,1,1,hinfsynOptions('RelTol',1e-7));     % MMP, 11/05/2026
    end
    fprintf('\nBaseline: finite-difference discretization and hinfsyn\n');
    for k = 1:numel(Nlist),     fprintf('  N = %3d:  gamma = %.4f\n',Nlist(k),gam_N(k));     end
%   fprintf('PIE:  gam (LPI) = %.4f, closed-loop norm bound = %.4f\n',gam,gam_cl); % MMP, 11/05/2026 (was)
    fprintf('PIE:  gam (LPI) = %.4f\n',gam);                                % MMP, 11/05/2026
    if cert_cl,     fprintf('      closed-loop norm bound = %.4f\n',gam_cl);    end % MMP, 11/05/2026
end
