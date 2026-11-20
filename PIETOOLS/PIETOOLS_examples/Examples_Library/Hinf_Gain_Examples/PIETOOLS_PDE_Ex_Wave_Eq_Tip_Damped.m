function [PDE_t,PDE_b,st] = PIETOOLS_PDE_Ex_Wave_Eq_Tip_Damped(GUI,params)  % MMP, 11/19/2026
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% PIETOOLS PDE Examples
% [PDE_t,PDE_b,ST] = PIETOOLS_PDE_EX_WAVE_EQ_TIP_DAMPED(GUI,PARAMS) returns
% the tip-damped wave equation below, with disturbance w and integral
% output z, in the term-based and in the batch format, and the recommended
% settings ST of its H-infinity gain test. It sets Hinf_gain = 1,
% Hinf_gain_dual = 1 and Hinf_gain_form = 'P' (storage <T x, R T x>,
% R >= eppos I; dual form 'Pd') in the base workspace. With GUI = 1 it
% also loads the associated GUI save file into a new PIETOOLS_PDE_GUI
% window.
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
% % Example tip damped wave equation:
% % PDE         phi_{tt} = phi_{ss} + w(t)
% % With BCs    phi(s=0) = 0
% %             phi_{s}(s=1) = -k*phi_{t}(s=1)
% % And output  z(t) = int(phi_{t}(t,s),s,0,1)
% %
% % We use states x1 = phi_{s}, x2 = phi_{t}
% % Then x2(0) = 0,  x1(1) + k*x2(1) = 0.
% %
% % Parameter k can be set.
% % gamma = 0.8419 for k = 0.5.
% %---------------------------------------------------------------------% %

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PIETOOLS - PIETOOLS_PDE_Ex_Wave_Eq_Tip_Damped
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
% function [PDE_t,PDE_b] = PIETOOLS_PDE_Ex_Wave_Eq_Tip_Damped(GUI,params)   % MMP, 11/19/2026 (was)
% MMP, 11/19/2026: Release comment pass (comment-only).

% Determine the location of this example file <-- DO NOT MOVE THE FILE
loc = mfilename('fullpath');
root = fileparts(loc);

% Initialize variables
pvar s

%%% Executive Function:
evalin('base','Hinf_gain = 1;')
evalin('base','Hinf_gain_dual = 1;')
evalin('base','Hinf_gain_form = ''P'';');  % form of the gain test          % MMP, 11/19/2026

% Default parameter, overridden by the assignment strings in params
k = 0.5;
npars = length(params);
if npars~=0
    %%% Specify potential parameters
    for j=1:npars
        eval(params{j});
    end
end


% % % Construct the PDE.
%%% Batch input format
PDE_b.n0 = 0;   PDE_b.n1 = 2;   PDE_b.n2 = 0;   PDE_b.nw = 1;   PDE_b.nz = 1;
PDE_b.dom = [0,1];
np = PDE_b.n0 + PDE_b.n1 + PDE_b.n2;

PDE_b.A1= [0 1; 1 0];
PDE_b.B21 = [0; 1];
PDE_b.Ca1 = [0 1];

PDE_b.B = [0 1 0 0; 0 0 1 k];


%%% Term-based input format
% Initialize 1D PDE state components.
PDE_t.x{1}.vars = s;   PDE_t.x{1}.dom = [0,1];
PDE_t.x{2}.vars = s;   PDE_t.x{2}.dom = [0,1];
% Initialize finite-dimensional inputs and outputs.
PDE_t.w{1}.vars = [];  PDE_t.z{1}.vars = [];

% PDE: x1_{t} = x2_{s}
PDE_t.x{1}.term{1}.x = 2;
PDE_t.x{1}.term{1}.D = 1;

% PDE: x2_{t} = x1_{s}
PDE_t.x{2}.term{1}.x = 1;
PDE_t.x{2}.term{1}.D = 1;

% PDE: x2_{t} = ... + w
PDE_t.x{2}.term{2}.w = 1;

% Output: z = int_{0}^{1} x2(s) ds
PDE_t.z{1}.term{1}.x = 2;
PDE_t.z{1}.term{1}.I{1} = PDE_t.x{2}.dom;

% BC 1: 0 = x2(0)
PDE_t.BC{1}.term{1}.x = 2;
PDE_t.BC{1}.term{1}.loc = 0;

% BC 2: 0 = x1(1)...
PDE_t.BC{2}.term{1}.x = 1;
PDE_t.BC{2}.term{1}.loc = 1;
%...+ k*x2(1)
PDE_t.BC{2}.term{2}.x=2;
PDE_t.BC{2}.term{2}.loc=1;
PDE_t.BC{2}.term{2}.C=k;

if GUI
    %%% Associated GUI save file
    app = PIETOOLS_PDE_GUI;
    load(fullfile(root,'PIETOOLS_PDE_Ex_Wave_Eq_Tip_Damped_GUI.mat'));
    logval = app.loadData(data);
    if logval
        disp("Failed to load data object. Incorrect structure");
    end
end

% Example-specific settings                                                 % MMP, 11/19/2026
% (stripped tier; slack lift degree at its minimum Dmin)
st = lpisettings_sop('stripped');                                           % MMP, 11/19/2026
st.lpi.ineq.dD = 0;                                                         % MMP, 11/19/2026

end