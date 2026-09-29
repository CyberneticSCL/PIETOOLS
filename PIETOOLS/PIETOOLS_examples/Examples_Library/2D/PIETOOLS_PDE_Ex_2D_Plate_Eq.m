function [PDE_t] = PIETOOLS_PDE_Ex_2D_Plate_Eq(GUI,params)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% PIETOOLS PDE Examples
% INPUT
% - GUI:        Binary index {0,1} indicating whether or not a GUI
%               implementation of the example should be produced.
% - params:     Optional parameters for the example, that should be
%               specified as a cell of strings e.g. {'alp0=0.2;'}.
%
% OUTPUT
% - PDE_t:      PDE structure defining the example system in the term-based
%               format.
%
% %---------------------------------------------------------------------% %
% % % Clamped Kirchhoff plate with structural damping (Jagt & Peet, 2025,
% % % Sec. 7.2.3, Eqns. (20)-(21)):
% % % PDE         w_{tt} = -w_{s1s1s1s1} - 2*w_{s1s1s2s2} - w_{s2s2s2s2}
% % %                       + alp0*(w_{t s1s1} + w_{t s2s2}),
% % % With BCs    w(s1=0) = w_{s1}(s1=0) = w(s1=1) = w_{s1}(s1=1) = 0;
% % %             w(s2=0) = w_{s2}(s2=0) = w(s2=1) = w_{s2}(s2=1) = 0,
% % %             (s1,s2) in [0,1]^2.
% % %
% % % Use the state u = [u1; u2] = [w; w_{t}], both in S_2^{(4,4)}, then
% % %
% % % PDE         u_{t} = [0, 1; 0, 0]*u - [0, 0; 1, 0]*(u_{s1s1s1s1} + 2*u_{s1s1s2s2} + u_{s2s2s2s2})
% % %                     + [0, 0; 0, alp0]*(u_{s1s1} + u_{s2s2})
% % % With BCs    clamped on all four edges for the whole vector u, as in (21)
% % %             of the paper (so w_t is clamped too).
% %
% % Parameter alp0 (> 0) can be set. The paper certifies an exponential
% % PIE-to-PDE decay rate of at least 3.6328 at alp0 = 0.2 (Cor. 35, d = 0,
% % eps = 0.1, MOSEK); the paper gives no analytic rate. As for the damped wave, the
% % stock PIETOOLS_stability_2D test is not the right one; the paper's LPI is
% % run by PIETOOLS_demos/multivariate_PIE_2508_14840/jp_wave_plate_rates.m.
% %---------------------------------------------------------------------% %
%
% CC, 09/28/2026: new example, recreating the clamped damped plate of
% D. S. Jagt and M. M. Peet, "A State-Space Representation of Coupled Linear
% Multivariate PDEs and Stability Analysis using SDP", arXiv:2508.14840 (v4).
% It is the first 2D fourth-order PDE in the example library. No executive
% flag is set (see the damped wave example for why).

% Initialize variables
pvar s1 s2

alp0 = 0.2;
npars = length(params);
if npars~=0
    %%% Specify potential parameters
    for j=1:npars
        eval(params{j});
    end
end

% % % Construct the PDE.
clear stateNameGenerator
u = pde_var(2,[s1;s2],[0,1;0,1],[4;4]);   % u = [w; w_t], order (4,4) in both components
PDE_t = [diff(u,'t')==[0,1;0,0]*u - [0,0;1,0]*(diff(u,s1,4)+2*diff(u,[s1;s2],[2;2])+diff(u,s2,4)) ...
                      + [0,0;0,alp0]*(diff(u,s1,2)+diff(u,s2,2));
         subs(u,s1,0)==0;     subs(diff(u,s1),s1,0)==0;
         subs(u,s1,1)==0;     subs(diff(u,s1),s1,1)==0;
         subs(u,s2,0)==0;     subs(diff(u,s2),s2,0)==0;
         subs(u,s2,1)==0;     subs(diff(u,s2),s2,1)==0];

if GUI~=0
    disp('No GUI representation available for this system.')
end

end
% D. S. Jagt and M. M. Peet, "A State-Space Representation of Coupled Linear
% Multivariate PDEs and Stability Analysis using SDP," arXiv:2508.14840, 2025.
