function [PDE_t,PDE_b,st] = PIETOOLS_PDE_Ex_Heat_Eq_with_ODE_Tang(GUI,params)  % MMP, 11/19/2026
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% PIETOOLS PDE Examples
% [PDE_t,PDE_b,ST] = PIETOOLS_PDE_EX_HEAT_EQ_WITH_ODE_TANG(GUI,PARAMS)
% returns the ODE-heat equation system below, coupled through both
% boundaries, in the term-based and in the batch format, and the
% recommended settings ST of its stability test. It sets stability = 1
% and stability_form = 'Q' (generalized storage T'Q) in the base
% workspace. With GUI = 1 it also loads the associated GUI save file into
% a new PIETOOLS_PDE_GUI window.
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
%               for the tests this file names                               % MMP, 11/19/2026
%
% %---------------------------------------------------------------------% %
% % Coupled ODE-PDE system from Tang (see reference below):
% % ODE             xo_{t} = xo + x_{s}(s=0)
% % PDE             x_{t} = x_{ss}
% % With BCs        x(s=0) = -xo
% %                 x(s=1) = k * xo
% %
% % Parameter k can be set.
% %---------------------------------------------------------------------% %

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PIETOOLS - PIETOOLS_PDE_Ex_Heat_Eq_with_ODE_Tang
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
% MMP, 11/19/2026: Output st, the example-specific settings, and the form
%                   of its tests named in the workspace (stability_form,
%                   Hinf_gain_form, H2_norm_form), from the search of
%                   11/19/2026 (maintainer: no default settings);
% function [PDE_t,PDE_b] = PIETOOLS_PDE_Ex_Heat_Eq_with_ODE_Tang(GUI,params)  % MMP, 11/19/2026 (was)
% MMP, 11/19/2026: Release comment pass (comment-only).

% Determine the location of this example file <-- DO NOT MOVE THE FILE
loc = mfilename('fullpath');
root = fileparts(loc);

% Initialize variables
pvar s

%%% Executive Function:
evalin('base','stability = 1;');
evalin('base','stability_form = ''Q'';');  % form of the stability test     % MMP, 11/19/2026

% Default parameter, overridden by the assignment strings in params
k = -2;
npars = length(params);
if npars~=0
    %%% Specify potential parameters
    for j=1:npars
        eval(params{j});
    end
end


% Construct the PDE.
%%% Batch input format
PDE_b.nx = 1;   PDE_b.nw = 0;   PDE_b.ny = 0;   PDE_b.nz = 0;   PDE_b.nu = 0;
PDE_b.n0 = 0;   PDE_b.n1 = 0;   PDE_b.n2 = 1;
PDE_b.dom = [0,1];

PDE_b.A = 1;   PDE_b.E0 = [0 0 1 0];
PDE_b.A2 = 1;
PDE_b.Bx = [-1; k];

PDE_b.B = [1 0 0 0;
           0 1 0 0];


%%% Term-based input format
% Initialize ODE state component.
PDE_t.x{1}.vars = [];
% Initialize 1D PDE state component.
PDE_t.x{2}.vars = s;   PDE_t.x{2}.dom = [0,1];

% ODE: xo_{t} = xo
PDE_t.x{1}.term{1}.x = 1;

% ODE: xo_{t} = ... + x_{s}(s=0)
PDE_t.x{1}.term{2}.x = 2;
PDE_t.x{1}.term{2}.D = 1;
PDE_t.x{1}.term{2}.loc = 0;

% PDE: x_{t} = x_{ss}
PDE_t.x{2}.term{1}.x = 2;
PDE_t.x{2}.term{1}.D = 2;

% BC 1: 0 = x(0)
PDE_t.BC{1}.term{1}.x = 2;
PDE_t.BC{1}.term{1}.loc = 0;

% BC 1: 0 = ... xo
PDE_t.BC{1}.term{2}.x = 1;
PDE_t.BC{1}.term{2}.C = 1;

% BC 2: 0 = -x(1)
PDE_t.BC{2}.term{1}.x = 2;
PDE_t.BC{2}.term{1}.loc = 1;
PDE_t.BC{2}.term{1}.C = -1;

% BC 2: 0 = ... k * xo
PDE_t.BC{2}.term{2}.x = 1;
PDE_t.BC{2}.term{2}.C = k;


if GUI
    %%% Associated GUI save file
    app = PIETOOLS_PDE_GUI;
    load(fullfile(root,'PIETOOLS_PDE_Ex_Heat_Eq_with_ODE_Tang_GUI.mat'));
    logval = app.loadData(data);
    if logval
        disp("Failed to load data object. Incorrect structure");
    end
end

% Example-specific settings                                                 % MMP, 11/19/2026
% (stripped tier; slack lift degree at its minimum Dmin; slack integral
% bases capped at total degree theta+s <= Dmax-1; Q of degree 2 in the
% multiplier and 1 in the integral kernels)
st = lpisettings_sop('stripped');                                           % MMP, 11/19/2026
st.lpi.ineq.dD = 0;                                                         % MMP, 11/19/2026
st.lpi.ineq.cap_inc = 'auto';                                               % MMP, 11/19/2026
st.sop.deg_Q{1} = [2 1 1];                                                  % MMP, 11/19/2026

end
% @article{tang2011state,
% title = {State and output feedback boundary control for a coupled PDE–ODE system},
% journal = {Systems & Control Letters},
% volume = {60},
% number = {8},
% pages = {540-545},
% year = {2011},
% issn = {0167-6911},
% doi = {https://doi.org/10.1016/j.sysconle.2011.04.011},
% url = {https://www.sciencedirect.com/science/article/pii/S0167691111000922},
% author = {Shuxia Tang and Chengkang Xie},
% keywords = {Coupled system, PDEs, Boundary control, Output feedback},
% abstract = {This note is devoted to stabilizing a coupled PDE–ODE system with interaction at the interface. First, a state feedback boundary controller is designed, and the system is transformed into an exponentially stable PDE–ODE cascade with an invertible integral transformation, where PDE backstepping is employed. Moreover, the solution to the resulting closed-loop system is derived explicitly. Second, an observer is proposed, which is proved to exhibit good performance in estimating the original coupled system, and then an output feedback boundary controller is obtained. For both the state and output feedback boundary controllers, exponential stability analyses in the sense of the corresponding norms for the resulting closed-loop systems are provided. The boundary controller and observer for a scalar coupled PDE–ODE system as well as the solutions to the closed-loop systems are given explicitly.}
% }