function [PDE_t,PDE_b,st] = PIETOOLS_PDE_Ex_Heat_Eq_with_ODE_Das(GUI,params) % MMP, 11/19/2026
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% PIETOOLS PDE Examples
% [PDE_t,PDE_b,ST] = PIETOOLS_PDE_EX_HEAT_EQ_WITH_ODE_DAS(GUI,PARAMS)
% returns the heat equation coupled with a 4-state ODE below (Das et al.,
% Example 2), in the term-based and in the batch format, and the
% recommended settings ST of its stability tests. It sets stability = 1,
% stability_dual = 1 and stability_form = 'Q' (generalized storage T'Q;
% dual form 'Qd') in the base workspace. With GUI = 1 it also loads the
% associated GUI save file into a new PIETOOLS_PDE_GUI window.
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
% - PDE_b:      PDE structure defining the example system in the batch
%               format.
% - st:         example-specific settings (lpisettings_sop structure)       % MMP, 11/19/2026
%               for the stability tests                                     % MMP, 11/19/2026
%
% %---------------------------------------------------------------------% %
% % Heat Equation coupled with ODE from Das (Example 2):
% % % ODE        xo_{t} = A * xo + Bxr * x_{s}(s=a)                         % MMP, 11/19/2026 (was)
% % ODE        xo_{t} = A * xo - Bxr * x_{s}(s=a)                           % MMP, 11/19/2026
% % PDE        x_{t} = lam * x + x_{ss} + Bpv * xo
% % with BCs   x(s=a) = 0
% %            x(s=b) = 0
% %
% %
% % Parameters lam, a and b (domain) can be set.
% % Stable for lam<pi^2                                                     % MMP, 11/19/2026 (was)
% % Stable for lam < 9.7946. The paper (Das et al. 2018, arXiv 1812.07186,  % MMP, 11/19/2026
% % Sec. 6.2) has +Bxr, which is stable only for lam < -5.7799 against its  % MMP, 11/19/2026
% % claim lam < pi^2; the sign of Bxr is taken as a typo.                   % MMP, 11/19/2026
% %---------------------------------------------------------------------% %

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PIETOOLS - PIETOOLS_PDE_Ex_Heat_Eq_with_ODE_Das
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
% Initial coding DJ - 08/20/2022
% MMP, 11/19/2026: the flux input of the ODE takes the sign -Bxr, in both
%                   formats (maintainer: "probably a typo in the paper. so
%                   reverse the sign"). With +Bxr the system is stable only
%                   for lam < -5.7799 (Chebyshev, 56 and 72 nodes), against
%                   the lam < pi^2 of the paper; with -Bxr, lam < 9.7946.
% MMP, 11/19/2026: output st, the example-specific settings of the two
%                   stability tests, and their form ('Q'), from the search
%                   of 11/19/2026 (maintainer: no default settings);
% function [PDE_t,PDE_b] = PIETOOLS_PDE_Ex_Heat_Eq_with_ODE_Das(GUI,params) % MMP, 11/19/2026 (was)
% MMP, 11/19/2026: Release comment pass (comment-only).
% MMP, 11/21/2026: The GUI file PIETOOLS_PDE_Ex_Heat_Eq_with_ODE_Das_GUI.mat
%                   now stores -Bxr, as PDE_t and PDE_b: the entry of
%                   11/19/2026 changed the code only, so GUI = 1 loaded the
%                   +Bxr system, unstable at the default lam.

% Determine the location of this example file <-- DO NOT MOVE THE FILE
loc = mfilename('fullpath');
root = fileparts(loc);

% Initialize variables
pvar s

%%% Executive Function:
evalin('base','stability = 1;');
evalin('base','stability_dual = 1;')
evalin('base','stability_form = ''Q'';');  % form of the stability tests    % MMP, 11/19/2026

% Default parameters, overridden by the assignment strings in params
lam = pi^2-1;  a = 0;  b = 1;
npars = length(params);
if npars~=0
    %%% Specify potential parameters
    for j=1:npars
        eval(params{j});
    end
end


% Construct the PDE.
%%% Batch input format
PDE_b.n0 = 0;   PDE_b.n1 = 0;   PDE_b.n2 = 2;   PDE_b.nx = 4;
PDE_b.dom = [a,b];

% ODE
% ODE dynamic matrix
PDE_b.A = [-1.2142,  1.9649,  0.2232,  0.5616;
    -1.8042, -0.7260, -0.3479,  5.4355;
    -0.2898,  0.7381, -1.7606,  0.8294;
    -0.9417, -5.3399, -1.0704, -0.7590];
% PDE state to ODE dynamics
%PDE_b.E0 = [-1.5368 0;0 0.8871;1.0656 0;1.1882 0] * [zeros(PDE_b.n2) zeros(PDE_b.n2) eye(PDE_b.n2) zeros(PDE_b.n2)]; % MMP, 11/19/2026 (was)
PDE_b.E0 = -[-1.5368 0;0 0.8871;1.0656 0;1.1882 0] * [zeros(PDE_b.n2) zeros(PDE_b.n2) eye(PDE_b.n2) zeros(PDE_b.n2)]; % MMP, 11/19/2026

% PDE
% PDE dynamics
PDE_b.A0 = lam*eye(PDE_b.n2);   PDE_b.A2 = eye(PDE_b.n2);
% ODE state to PDE dynamics
PDE_b.E = [-2.5575 0 1.0368 0;-1.8067 0.4630 1.3621 0];

% BCs
PDE_b.B = [eye(PDE_b.n2), zeros(PDE_b.n2), zeros(PDE_b.n2),zeros(PDE_b.n2);
    zeros(PDE_b.n2), eye(PDE_b.n2), zeros(PDE_b.n2), zeros(PDE_b.n2)];


%%% Term-based input format
% Initialize ODE state component.
PDE_t.x{1}.vars = [];
% Initialize PDE state component.
PDE_t.x{2}.vars = s;   PDE_t.x{2}.dom = [a,b];

% ODE: xo_{t} = A * xo;
PDE_t.x{1}.term{1}.x = 1;
PDE_t.x{1}.term{1}.C = [-1.2142,  1.9649,  0.2232,  0.5616;
    -1.8042, -0.7260, -0.3479,  5.4355;
    -0.2898,  0.7381, -1.7606,  0.8294;
    -0.9417, -5.3399, -1.0704, -0.7590];

% % ODE: xo_{t} = ... + Bxr * x_{s}(s=0)                                    % MMP, 11/19/2026 (was)
% ODE: xo_{t} = ... - Bxr * x_{s}(s=0)                                      % MMP, 11/19/2026
PDE_t.x{1}.term{2}.x = 2;
PDE_t.x{1}.term{2}.D = 1;
PDE_t.x{1}.term{2}.loc = 0;
%PDE_t.x{1}.term{2}.C = [-1.5368 0;0 0.8871;1.0656 0;1.1882 0];             % MMP, 11/19/2026 (was)
PDE_t.x{1}.term{2}.C = -[-1.5368 0;0 0.8871;1.0656 0;1.1882 0];             % MMP, 11/19/2026

% PDE: x_{t} = lam * x
PDE_t.x{2}.term{1}.x = 2;
PDE_t.x{2}.term{1}.C = lam*eye(2);

% PDE: x_{t} = ... + x_{ss}
PDE_t.x{2}.term{2}.x = 2;
PDE_t.x{2}.term{2}.D = 2;

% PDE: x_{t} = ... + Bpv * v
PDE_t.x{2}.term{3}.x = 1;
PDE_t.x{2}.term{3}.C = [-2.5575 0 1.0368 0;-1.8067 0.4630 1.3621 0];

% BCs: 0 = x(0)
PDE_t.BC{1}.term{1}.x = 2;
PDE_t.BC{1}.term{1}.loc = a;
PDE_t.BC{1}.term{1}.C = eye(2);

% BCs: 0 = x(1)
PDE_t.BC{2}.term{1}.x = 2;
PDE_t.BC{2}.term{1}.loc = b;
PDE_t.BC{2}.term{1}.C = eye(2);
% command-line format
% PDE_c=sys();
% xo=state('ode');x=state('pde');
% % eq_ODEPDE=[diff(xo,t) == A*xo+Bxr*subs(diff(x,s),s,a)                   % MMP, 11/19/2026 (was)
% eq_ODEPDE=[diff(xo,t) == A*xo-Bxr*subs(diff(x,s),s,a)                     % MMP, 11/19/2026
%                         diff(x,t) == lam*x+diff(x,s,2)+Bpv*xo];
% eq_BCs=[subs(x,s,a)==0;subs(x,s,b)==0];
% PDE_c=addequation(PDE_c,[eq_ODEPDE;eq_BCs]);
% Example-specific settings                                                 % MMP, 11/19/2026
% (stripped tier; Q of degree 2 in the multiplier and 1 in the integral
% kernels)
st = lpisettings_sop('stripped');                                           % MMP, 11/19/2026
st.sop.deg_Q{1} = [2 1 1];                                                  % MMP, 11/19/2026

if GUI
    %%% Associated GUI save file
    app = PIETOOLS_PDE_GUI;
    load(fullfile(root,'PIETOOLS_PDE_Ex_Heat_Eq_with_ODE_Das_GUI.mat'));
    logval = app.loadData(data);
    if logval
        disp("Failed to load data object. Incorrect structure");
    end
end

end
% @article{das2018representation,
% author = {Das, Amritam and Shivakumar, Sachin and Weiland, Siep and Peet, Matthew},
% year = {2018},
% pages = {},
% title = {Representation and Stability Analysis of PDE-ODE Coupled Systems}
% }