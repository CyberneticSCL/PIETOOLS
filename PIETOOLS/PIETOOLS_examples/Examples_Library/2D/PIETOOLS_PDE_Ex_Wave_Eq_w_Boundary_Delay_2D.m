function [PDE_t,st] = PIETOOLS_PDE_Ex_Wave_Eq_w_Boundary_Delay_2D(GUI,params) % MMP, 11/21/2026
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% PIETOOLS PDE Examples
% [PDE_t,ST] = PIETOOLS_PDE_EX_WAVE_EQ_W_BOUNDARY_DELAY_2D(GUI,PARAMS)
% returns the wave equation below with a delayed damping feedback at
% s1 = 1, the delayed boundary values written as transport states x4, x5
% in a second variable s2, as a 'pde_struct' object, and the recommended
% settings ST of its stability test. It sets stability = 1 and
% stability_form = 'Q' (generalized storage T'Q) in the base workspace; a
% nonzero GUI prints that no GUI representation exists.
%
% INPUT
% - GUI:        Binary index {0,1} indicating whether or not a GUI
%               implementation of the example should be produced.
% - params:     Optional parameters for the example, that should be
%               specified as a cell of strings e.g. {'lam=1;','dom=[0,1]'}.
%
% OUTPUT
% - PDE_t:      PDE structure defining the example system in the term-based 
%               format.
% - st:         example-specific settings (lpisettings_sop structure)       % MMP, 11/21/2026
%               for the tests this file names                               % MMP, 11/21/2026
%
% %---------------------------------------------------------------------% %
% % Wave Equation with Delay in BCs (Xu 2006):
% % PDE         x_{tt}  = x_{s1s1}(t,s1);       s1 in [0,1]   
% % With BCs    x(t,0) = 0;   
% %             x_{s1}(t,1) = -k*(1-mu)*x_{t}(t,1) - k*mu*x_{t}(t-tau,1);
% %
% % Introduce:
% %     x1(t)    = - k*(1-mu)*x(t,1) - k*mu*x(t-tau,1);
% %     x2(t,s1) = x(t,s1);
% %     x3(t,s1) = x_{t}(t,s1);
% %     x4(t,s2) = x(t-tau*s2,1);
% %     x5(t,s2) = x_{t}(t-tau*s2,1);
% % Then:
% % ODE:    x1_{t}(t)    = x2_{s1}(t,1);
% % PDE:    x2_{t}(t,s1) = x3(t,s1);                    s1 in (0,1)
% %         x3_{t}(t,s1) = x2_{s1s1}(t,s1);                
% %         x4_{t}(t,s2) = -(1/tau) x4_{s2}(t,s2);      s2 in (0,1)
% %         x5_{t}(t,s2) = -(1/tau) x5_{s2}(t,s2);
% % BCs:    x2(t,0) = 0;
% %         x3(t,0) = 0;
% %         x4(t,0) = x2(t,1);
% %         x5(t,0) = x3(t,1);
% %         x1(t)   = - k*(1-mu)*x2(t,1) - k*mu*x4(t,1);
% %         x2_{s1}(t,1) = - k*(1-mu)*x3(t,1) - k*mu*x5(t,1);
% %
% % Parameters tau, k, and mu can be set; a and b (default 0, 1) set the
% % domain of s1, but the BCs take x2, x3 at s1 = 1, so b = 1 is assumed.
% % Stable whenever mu < 1/2
% %---------------------------------------------------------------------% %

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PIETOOLS - PIETOOLS_PDE_Ex_Wave_Eq_w_Boundary_Delay_2D
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
% Initial coding DJ - 11/30/2022
% MMP, 11/21/2026: Changes since PIETOOLS 2025: second output st, the
%   recommended settings of the stability test; stability_form = 'Q' set in
%   the base workspace; help revised (OUTPUT lists no PDE_b, which is not
%   returned; the BCs assume b = 1), license added.
% function [PDE_t] = PIETOOLS_PDE_Ex_Wave_Eq_w_Boundary_Delay_2D(GUI,params) % MMP, 11/21/2026 (was)

% Determine the location of this example file <-- DO NOT MOVE THE FILE
loc = mfilename('fullpath');
root = fileparts(loc);

% Initialize variables
pvar s1 s2

%%% Executive Function:
evalin('base','stability = 1;');
evalin('base','stability_form = ''Q'';');  % form of the stability test     % MMP, 11/21/2026
%evalin('base','stability_dual = 1;');

% Default parameters, overridden by the assignment strings in params
a = 0;      b = 1;
k = 1;
tau = 1;    mu = 0.5-1e-1;
npars = length(params);
if npars~=0
    %%% Specify potential parameters
    for j=1:npars
        eval(params{j});
    end
end


% % % Construct the PDE.
%%% pde_var input format
clear stateNameGenerator
x1 = pde_var();
x2 = pde_var(s1,[a,b]);
x3 = pde_var(s1,[a,b],2);       % x3 must be 2nd order differentiable wrt s1
x4 = pde_var(s2,[0,1]);
x5 = pde_var(s2,[0,1]);
PDE_t = [diff(x1,'t')==subs(diff(x2,s1),s1,b);
         diff(x2,'t')==x3;
         diff(x3,'t')==diff(x2,s1,2);
         diff(x4,'t')==-(1/tau)*diff(x4,s2);
         diff(x5,'t')==-(1/tau)*diff(x5,s2);
         subs([x2;x3],s1,a)==0;
         x1==-k*(1-mu)*subs(x2,s1,1)-k*mu*subs(x4,s2,1);
         subs(diff(x2,s1),s1,b)==-k*(1-mu)*subs(x3,s1,1)-k*mu*subs(x5,s2,1);
         subs(x4,s2,0)==subs(x2,s1,b);
         subs(x5,s2,0)==subs(x3,s1,b)];


% %%% Term-based input format
% % Initialize the state variables and input.
% pde_struct PDE_t;
% PDE_t.x{1}.vars = [];
% PDE_t.x{2}.vars = s1;       PDE_t.x{2}.dom = [a,b];
% PDE_t.x{3}.vars = s1;       PDE_t.x{3}.dom = [a,b];
% PDE_t.x{3}.diff = 2;    % x3 must be 2nd order differentiable wrt s1
% PDE_t.x{4}.vars = s2;       PDE_t.x{4}.dom = [0,1];
% PDE_t.x{5}.vars = s2;       PDE_t.x{5}.dom = [0,1];
% 
% % PDE: x1_{t}(t) = x2_{s1}(t,1)
% PDE_t.x{1}.term{1}.x = 2;
% PDE_t.x{1}.term{1}.D = 1;
% PDE_t.x{1}.term{1}.loc = b;
% 
% % PDE: x2_{t}(t,s1) = x3(t,s1)
% PDE_t.x{2}.term{1}.x = 3;
% 
% % PDE: x3_{t}(t,s1) = x2_{s1s1}(t,s1);
% PDE_t.x{3}.term{1}.x = 2;
% PDE_t.x{3}.term{1}.D = 2;
% 
% % PDE: x4_{t}(t,s2) = -(1/tau)*x4_{s2}(t,s2);
% PDE_t.x{4}.term{1}.x = 4;
% PDE_t.x{4}.term{1}.D = 1;
% PDE_t.x{4}.term{1}.C = -1/tau;
% 
% % PDE: x5_{t}(t,s2) = -(1/tau)*x5_{s2}(t,s2);
% PDE_t.x{5}.term{1}.x = 5;
% PDE_t.x{5}.term{1}.D = 1;
% PDE_t.x{5}.term{1}.C = -1/tau;
% 
% 
% % BC1: x2(t,s1=a) = 0;              BC4: x3(t,s1=a) = 0;
% PDE_t.BC{1}.term{1}.x = 2;          PDE_t.BC{4}.term{1}.x = 3;
% PDE_t.BC{1}.term{1}.loc = a;        PDE_t.BC{4}.term{1}.loc = a;
% 
% % BC2: x1(t) = -k*(1-mu)*x2(t,s1=b) - k*mu*x4(t,s2=1)
% PDE_t.BC{2}.term{1}.x = 1;  PDE_t.BC{2}.term{2}.x = 2;        PDE_t.BC{2}.term{3}.x = 4;
%                             PDE_t.BC{2}.term{2}.C = k*(1-mu); PDE_t.BC{2}.term{3}.C = k*mu;
%                             PDE_t.BC{2}.term{2}.loc = b;      PDE_t.BC{2}.term{3}.loc = 1;
% 
% % BC3: x2_{s1}(t,s1=b) = -k*(1-mu)*x3(t,s1=b) - k*mu*x5(t,s2=1)
% PDE_t.BC{3}.term{1}.x = 2;    PDE_t.BC{3}.term{2}.x = 3;        PDE_t.BC{3}.term{3}.x = 5;
% PDE_t.BC{3}.term{1}.D = 1;    PDE_t.BC{3}.term{2}.C = k*(1-mu); PDE_t.BC{3}.term{3}.C = k*mu;
% PDE_t.BC{3}.term{1}.loc = b;  PDE_t.BC{3}.term{2}.loc = b;      PDE_t.BC{3}.term{3}.loc = 1;
% 
% % BC5: (Sewing condition) x4(t,s2=0) = x2(t,s1=b);
% PDE_t.BC{5}.term{1}.x = 4;      PDE_t.BC{5}.term{2}.x = 2;
% PDE_t.BC{5}.term{1}.loc = 0;    PDE_t.BC{5}.term{2}.loc = b;
%                                 PDE_t.BC{5}.term{2}.C = -1;
% 
% % BC6: (Sewing condition) x5(t,s2=0) = x3(t,s1=b);
% PDE_t.BC{6}.term{1}.x = 5;      PDE_t.BC{6}.term{2}.x = 3;
% PDE_t.BC{6}.term{1}.loc = 0;    PDE_t.BC{6}.term{2}.loc = b;
%                                 PDE_t.BC{6}.term{2}.C = -1;

                                
% Example-specific settings: the light tier                                 % MMP, 11/21/2026
st = lpisettings_sop('light');                                              % MMP, 11/21/2026

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
%      publisher = {EDP-Sciences},
%      volume = {12},
%      number = {4},
%      year = {2006},
%      doi = {10.1051/cocv:2006021},
%      zbl = {1105.35016},
%      mrnumber = {2266817},
%      language = {en},
%      url = {http://www.numdam.org/articles/10.1051/cocv:2006021/}
% }