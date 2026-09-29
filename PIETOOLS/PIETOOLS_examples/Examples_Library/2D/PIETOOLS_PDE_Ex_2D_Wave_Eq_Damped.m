function [PDE_t] = PIETOOLS_PDE_Ex_2D_Wave_Eq_Damped(GUI,params)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% PIETOOLS PDE Examples
% INPUT
% - GUI:        Binary index {0,1} indicating whether or not a GUI
%               implementation of the example should be produced.
% - params:     Optional parameters for the example, that should be
%               specified as a cell of strings e.g. {'kap=2;'}.
%
% OUTPUT
% - PDE_t:      PDE structure defining the example system in the term-based
%               format.
%
% %---------------------------------------------------------------------% %
% % % 2D wave equation with stabilizing feedback (Jagt & Peet, 2025, Sec. 7.2.2,
% % % Eqns. (18)-(19)):
% % % PDE         phi_{tt} = phi_{s1s1} + phi_{s2s2} - 2*kap*phi_{t} - kap^2*phi,
% % % With BCs    phi(s1=0) = 0;   phi_{s1}(s1=1) = 0;
% % %             phi(s2=0) = 0;   phi_{s2}(s2=1) = 0,       (s1,s2) in [0,1]^2.
% % %
% % % Use the state u = [u1; u2] = [phi; phi_{t}], both in S_2^{(2,2)}, then
% % %
% % % PDE         u_{t} = [0, 1; -kap^2, -2*kap]*u + [0, 0; 1, 0]*(u_{s1s1} + u_{s2s2})
% % % With BCs    u(s1=0) = 0;   u_{s1}(s1=1) = 0;
% % %             u(s2=0) = 0;   u_{s2}(s2=1) = 0,
% % %
% % % i.e. u2 = phi_t satisfies the same BCs and regularity as u1 = phi, as in
% % % (19) of the paper.
% %
% % Parameter kap (>= 0) can be set. For kap > 0 the system is exponentially
% % PIE-to-PDE stable with maximal rate kap (paper, Appx. C.2): phi =
% % exp(-kap*t)*psi with psi an undamped wave; kap = 0 is the undamped wave
% % itself, which is not exponentially stable. It is NOT exponentially stable
% % in the classical L2 sense for (phi,phi_t), so the stock PIETOOLS_stability_2D
% % test (P >= eppos*I) is not the right one; the paper certifies the rate
% % with the LPI of its Cor. 35, which is what
% % PIETOOLS_demos/multivariate_PIE_2508_14840/jp_wave_plate_rates.m runs.
% %
% % Reference values (paper, Table 2; d = 1, eps = 0.1, MOSEK): the maximal
% % certified rate for kap = 1..7 is 0.9999, 1.9998, 2.9874, 3.9973, 4.9980,
% % 5.9906, 6.9645.
% %---------------------------------------------------------------------% %
%
% CC, 09/28/2026: new example, recreating the damped 2D wave equation of
% D. S. Jagt and M. M. Peet, "A State-Space Representation of Coupled Linear
% Multivariate PDEs and Stability Analysis using SDP", arXiv:2508.14840 (v4).
% No executive flag is set (unlike the stability examples in this folder),
% because the default stability test cannot certify this system (see above).

% Initialize variables
pvar s1 s2

kap = 1;
npars = length(params);
if npars~=0
    %%% Specify potential parameters
    for j=1:npars
        eval(params{j});
    end
end

% % % Construct the PDE.
clear stateNameGenerator
u = pde_var(2,[s1;s2],[0,1;0,1],[2;2]);   % u = [phi; phi_t], order (2,2) in both components
PDE_t = [diff(u,'t')==[0,1;-kap^2,-2*kap]*u + [0,0;1,0]*(diff(u,s1,2)+diff(u,s2,2));
         subs(u,s1,0)==0;     subs(diff(u,s1),s1,1)==0;
         subs(u,s2,0)==0;     subs(diff(u,s2),s2,1)==0];

if GUI~=0
    disp('No GUI representation available for this system.')
end

end
% D. S. Jagt and M. M. Peet, "A State-Space Representation of Coupled Linear
% Multivariate PDEs and Stability Analysis using SDP," arXiv:2508.14840, 2025.
