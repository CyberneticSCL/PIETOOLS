function [PDE_t,st] = PIETOOLS_PDE_Ex_Wave_Delayed_Boundary_Damping_Xu(GUI,params)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% PIETOOLS PDE Examples
% [PDE_t,ST] = PIETOOLS_PDE_EX_WAVE_DELAYED_BOUNDARY_DAMPING_XU(GUI,PARAMS)
% returns the wave equation with partially delayed boundary damping below,
% in Riemann form with the delayed boundary velocity as a transport state
% y(s2), as a 'pde_struct' object, and the recommended settings ST of its
% stability test. It sets stability = 1 and stability_form = 'Q'
% (generalized storage T'Q) in the base workspace; a nonzero GUI prints
% that no GUI representation exists.
%
% INPUT
% - GUI:        Binary index {0,1} indicating whether or not a GUI
%               implementation of the example should be produced.
% - params:     Optional parameters for the example, that should be
%               specified as a cell of strings e.g. {'mu=0.75;','tau=1;'}.
%
% OUTPUT
% - PDE_t:      PDE structure defining the example system in the term-based
%               format.
% - st:         recommended settings (lpisettings structure) for the
%               test this file names.
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
% %---------------------------------------------------------------------% %

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PIETOOLS - PIETOOLS_PDE_Ex_Wave_Delayed_Boundary_Damping_Xu
%
% Copyright (C) 2026 PIETOOLS Team
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
% Initial coding MMP - 11/17/2026

% Determine the location of this example file <-- DO NOT MOVE THE FILE
loc = mfilename('fullpath');
root = fileparts(loc);

% Initialize variables
pvar s1 s2

%%% Executive Function:
evalin('base','stability = 1;');
evalin('base','stability_form = ''Q'';');  % form of the stability test

% Default parameters, overridden by the assignment strings in params
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


% Example-specific settings: the light tier, with the slack lift degree at
% its minimum Dmin, the 2-D slack with the plain term only (psatz 'none'),
% Q of degree 0, and up to 2 raises in one or two variables, each raising
% the degree of Q (raise_Q 'deg')
st = lpisettings('light');
st.lpi.ineq.dD = 0;
st.lpi.ineq.psatz{2} = 'none';
st.sop.deg_Q{2} = 0;
st.sop.max_raise_Q = [2 2 0];
st.sop.raise_Q = 'deg';

% No GUI representation of this system
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
