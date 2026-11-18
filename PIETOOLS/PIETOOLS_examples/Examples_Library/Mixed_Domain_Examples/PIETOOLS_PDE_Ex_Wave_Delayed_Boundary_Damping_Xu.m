function [PDE_t] = PIETOOLS_PDE_Ex_Wave_Delayed_Boundary_Damping_Xu(GUI,params)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% PIETOOLS PDE Examples
% INPUT
% - GUI:        Binary index {0,1} indicating whether or not a GUI
%               implementation of the example should be produced.
% - params:     Optional parameters for the example, that should be
%               specified as a cell of strings e.g. {'mu=0.75;','tau=1;'}.
%
% OUTPUT
% - PDE_t:      PDE structure defining the example system in the term-based
%               format.
%
% %---------------------------------------------------------------------% %
% % Wave equation with partially delayed boundary damping (Xu, Yung, Li
% % 2006):
% % PDE         w_{tt}(t,s1) = w_{s1s1}(t,s1);                s1 in [0,1]
% % With BCs    w(t,s1=0) = 0;
% %             w_{s1}(t,1) = -k*mu*w_{t}(t,1) - k*(1-mu)*w_{t}(t-tau,1);
% %
% % Let e = w_{s1}, v = w_{t} (Riemann form), and y(t,s2) = v(t-tau*s2,1),
% % s2 in [0,1]:
% % PDE         e_{t} = v_{s1};   v_{t} = e_{s1};   y_{t} = -(1/tau)*y_{s2};
% % With BCs    v(t,s1=0) = 0;
% %             e(t,s1=1) = -k*mu*v(t,1) - k*(1-mu)*y(t,s2=1);
% %             y(t,s2=0) = v(t,s1=1);
% %
% % Parameters mu, k and tau can be set; default mu = 0.6, k = 1, tau = 2.
% %
% % Threshold [abstract]: for mu > 1/2 some gain k gives exponential
% % stability for every delay; for mu = 1/2 at most asymptotic stability;
% % for mu < 1/2 instability.
% % Spectrum [derived, exact]: w = sinh(lam*s1) gives
% % cosh(lam) + k*sinh(lam)*(mu + (1-mu)*exp(-lam*tau)) = 0; for tau = 2,
% % X = exp(2*lam) solves a quadratic. For k = 1, tau = 2 the spectral
% % abscissa is 0.0957 (mu = 0.4), 0 (mu = 0.5), -0.105968 (mu = 0.6),
% % -0.301728 (mu = 0.75): exponentially stable iff mu > mu* = 1/2.
% %
% % Certified (PIETOOLS_stability_sop, 'P' form, lpisettings_sop('light'),
% % MOSEK, 11/17/2026), k = 1, tau = 2:
% %   - non-increase of a weighted norm of Tx (epneg = 0): certified at
% %     mu = 0.5, 0.500391, 0.6, 0.75 (stock light, 47,304 decision
% %     variables); infeasible at mu = 0.40, 0.45, 0.49. The proposed light
% %     row of the settings proposal of 11/15/2026 gives the same verdicts
% %     (50,580 variables at mu = 0.6).
% %   - decay rate (epneg > 0), mu = 0.6, stock light: certified for
% %     epneg = 0.10, infeasible for epneg = 0.11 (exact rate 0.105968).
% % Library example 44 treats the same paper in second-order form.
% %                                                                         % MMP, 11/17/2026
% % No second output st: the light tier certifies this example              % MMP, 11/17/2026
% % (record above). Its storage is the proposed light row since             % MMP, 11/17/2026
% % 11/17/2026; its slack of 11/17/2026 (wcross 0, 2-D faces at offset      % MMP, 11/17/2026
% % 1) was not rerun.                                                       % MMP, 11/17/2026
% %---------------------------------------------------------------------% %
%
% MMP, 11/17/2026: Initial coding;
% MMP, 11/17/2026: Header: the light tier certifies the example, so
%                   no second output st (application-specific settings);

% Determine the location of this example file <-- DO NOT MOVE THE FILE
loc = mfilename('fullpath');
root = fileparts(loc);

% Initialize variables
pvar s1 s2

%%% Executive Function:
evalin('base','stability = 1;');

mu = 0.6;   k = 1;  tau = 2;
npars = length(params);
if npars~=0
    %%% Specify potential parameters
    for j=1:npars
        eval(params{j});
    end
end


% % % Construct the PDE.
%%% pde_var input
clear stateNameGenerator
e = pde_var(s1,[0,1],1);
v = pde_var(s1,[0,1],1);
y = pde_var(s2,[0,1],1);
PDE_t = [diff(e,'t')==diff(v,s1);
         diff(v,'t')==diff(e,s1);
         diff(y,'t')==-(1/tau)*diff(y,s2);
         subs(v,s1,0)==0;
         subs(e,s1,1)==-k*mu*subs(v,s1,1)-k*(1-mu)*subs(y,s2,1);
         subs(y,s2,0)==subs(v,s1,1)];


if GUI~=0
    disp('No GUI representation available for this system.')
end

end
% @article{XU_2006,
%      author = {Xu, Gen Qi and Yung, Siu Pang and Li, Leong Kwan},
%      title = {Stabilization of wave systems with input delay in the boundary control},
%      journal = {ESAIM: Control, Optimisation and Calculus of Variations},
%      pages = {770--785},
%      volume = {12},
%      number = {4},
%      year = {2006},
%      doi = {10.1051/cocv:2006021}
% }
