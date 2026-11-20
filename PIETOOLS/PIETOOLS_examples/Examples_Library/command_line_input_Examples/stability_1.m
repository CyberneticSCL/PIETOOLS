% STABILITY_1 declares the transport equation x_t = v*x_s on [0,1] with
% x(0) = 0 (v = -1 by default), displays it, converts it to a PIE and
% displays the PIE. Workspace: PDE ('pde_struct') and PIE ('pie_struct').
% If a cell array params of assignment strings (e.g. {'v=-2;'}) exists in
% the workspace, its entries are evaluated after the default parameters;
% if a variable GUI exists and is true, the associated GUI save file is
% loaded into a new PIETOOLS_PDE_GUI window.
%
% 	PDE: x_{t} = v*x_{s}                            | v=-1;   (Exponentially PDE stable ) 
%   BCs: x(s=0) = 0                                 |                (Finite Energy PDE stable ) 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PIETOOLS - stability_1
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
% Initial coding SS - 11/16/2022
% MMP, 11/19/2026: Release comment pass (comment-only).

%clc; clear; 
clear stateNameGenerator
pvar s t; % define independent variables
%% Define dependent variables and system variable
% 	PDE: x_{t} = x_{s}                                 
%   BCs: x(s=0) = 0 
x = pde_var(s,[0,1]);
%% Specify the parameters
v=-1;
if exist('params','var')
    npars = length(params);
    %%% Specify potential parameters
    for j=1:npars
        eval(params{j});
    end
end
%% Define equations
eq_PDE = diff(x,t)==v*diff(x,s); % 	PDE: x_{t} = -x_{s}
eq_BC = subs(x,s,0)==0;        %   BCs: x(s=0) = 0 
%% initialize pde system;
PDE = initialize([eq_PDE;eq_BC]);
%% display pde to verify and convert to pie
display(PDE);
PIE = convert(PDE,'pie');
display(PIE);
%% Load the GUI save file of the example, if requested
if exist('GUI','var')
    if GUI
            loc = mfilename('fullpath');
        if isempty(loc)
            loc = which('stability_1');
        end
        script_dir = fileparts(loc);
        examples_dir = fileparts(script_dir);
        gui_dir = fullfile(examples_dir,'Stability_Examples');
        %%% Associated GUI save file
        app = PIETOOLS_PDE_GUI;
        load(fullfile(gui_dir,'PIETOOLS_PDE_Example_Transport_Eq_GUI.mat'));
        logval = app.loadData(data);
        if logval
            disp("Failed to load data object. Incorrect structure");
        end
    end
end


% %% Define dependent variables and system variable
% % 	PDE: x_{t} = -x_{s}                                 
% %   BCs: x(s=0) = 0 
% x = state('pde');
% pde = sys();
% %% Define equations
% eq_PDE = diff(x,t)==-diff(x,s); % 	PDE: x_{t} = x_{s}
% eq_BC = subs(x,s,0)==0;        %   BCs: x(s=0) = 0 
% %% addequations to pde system; set control inputs/observed inputs, if any
% pde = addequation(pde,[eq_PDE;eq_BC]);
% %% display pde to verify and convert to pie
% display(pde);
% pie = convert(pde,'pie');
% display(pie);
