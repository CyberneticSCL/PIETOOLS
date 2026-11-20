function [PDE_t,PDE_b,st] = PIETOOLS_PDE_Ex_Reaction_Diffusion_Ahmadi_2(GUI,params)  % MMP, 11/19/2026
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% PIETOOLS PDE Examples
% [PDE_t,PDE_b,ST] = PIETOOLS_PDE_EX_REACTION_DIFFUSION_AHMADI_2(GUI,
% PARAMS) returns the 3-state reaction-diffusion equation below with an
% s-dependent coupling Cm, in the term-based and in the batch format, and
% the recommended settings ST of its stability test. It sets stability =
% 1, stability_dual = 0 and stability_form = 'P' (storage <T x, P T x>,
% P >= eppos I) in the base workspace. With GUI = 1 it also loads the
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
%               for the tests this file names                               % MMP, 11/19/2026
%
% %---------------------------------------------------------------------% %
% % Adapted Example B (Diffusion-Reaction equation) from Ahmadi:
% % PDE        x_{t} = Cm*x + (1/R)*x_{ss} 
% % with BCs   x(s=0) = 0
% %            x_{s}(s=1) = 0
% %
% % Parameters Cm and R can be set.
% % Stable for R<21.
% %---------------------------------------------------------------------% %

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PIETOOLS - PIETOOLS_PDE_Ex_Reaction_Diffusion_Ahmadi_2
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
% function [PDE_t,PDE_b] = PIETOOLS_PDE_Ex_Reaction_Diffusion_Ahmadi_2(GUI,params)  % MMP, 11/19/2026 (was)
% MMP, 11/19/2026: Release comment pass (comment-only).

% Determine the location of this example file <-- DO NOT MOVE THE FILE
loc = mfilename('fullpath');
root = fileparts(loc);

% Initialize variables
pvar s

%%% Executive Function:
evalin('base','stability = 1;');
evalin('base','stability_form = ''P'';');  % form of the stability test     % MMP, 11/19/2026
evalin('base','stability_dual = 0;')

% Specify the parameters
Cm = [0 0 0; s 0 0; s -s^2 0];      R = 30;
npars = length(params);
if npars~=0
    %%% Specify potential parameters
    for j=1:npars
        eval(params{j});
    end
end


% Construct the PDE.
ne = size(Cm,2);
on = eye(ne);      ze = zeros(ne);


%%% Batch input format
PDE_b.n0 = 0;   PDE_b.n1 = 0;   PDE_b.n2 = ne;
PDE_b.dom = [0,1];

PDE_b.A0 = Cm;
PDE_b.A2 = (1/R)*on;

PDE_b.B = [on ze ze ze;
    ze on ze ze];


%%% Term-based input format
PDE_t.x{1}.vars = s;   PDE_t.x{1}.dom = [0,1];

% PDE: x_{t} = Cm * x
PDE_t.x{1}.term{1}.C = Cm;

% PDE: x_{t} = ... + (1/R) * x_{ss}
PDE_t.x{1}.term{2}.D = 2;
PDE_t.x{1}.term{2}.C = (1/R)*eye(ne);

% BCs: 0 = x(0)
PDE_t.BC{1}.term{1}.loc = 0;

% BCs: 0 = x_{s}(1)
PDE_t.BC{2}.term{1}.loc = 1;
PDE_t.BC{2}.term{1}.D = 1;


if GUI
    %%% Associated GUI save file
    app = PIETOOLS_PDE_GUI;
    load(fullfile(root,'PIETOOLS_PDE_Ex_Reaction_Diffusion_Ahmadi_2_GUI.mat'));
    logval = app.loadData(data);
    if logval
        disp("Failed to load data object. Incorrect structure");
    end
end

% Example-specific settings                                                 % MMP, 11/19/2026
% (light tier; slack lift degree at its minimum Dmin)
st = lpisettings_sop('light');                                              % MMP, 11/19/2026
st.lpi.ineq.dD = 0;                                                         % MMP, 11/19/2026

end
% @article{valmorbida2015stability,
%   title={Stability analysis for a class of partial differential equations via semidefinite programming},
%   author={Valmorbida, Giorgio and Ahmadi, Mohamadreza and Papachristodoulou, Antonis},
%   journal={IEEE Transactions on Automatic Control},
%   volume={61},
%   number={6},
%   pages={1649--1654},
%   year={2015},
%   publisher={IEEE}
% }