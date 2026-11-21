function [PDE_t,PDE_b,st] = PIETOOLS_PDE_Ex_Heat_Eq_with_Distributed_Disturbance(GUI,params) % MMP, 11/21/2026
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% PIETOOLS PDE Examples
% [PDE_t,PDE_b,ST] = PIETOOLS_PDE_EX_HEAT_EQ_WITH_DISTRIBUTED_DISTURBANCE(GUI,
% PARAMS) returns the heat equation below, with a distributed disturbance
% s*w and an integral output z, in the term-based and in the batch format,
% and the recommended settings ST of its H-infinity gain test. It sets
% Hinf_gain = 1 and Hinf_gain_form = 'P' (storage <T x, R T x>,
% R >= eppos I) in the base workspace. With GUI = 1 it also loads the
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
% - st:         example-specific settings (lpisettings_sop structure)       % MMP, 11/21/2026
%               for the tests this file names                               % MMP, 11/21/2026
%
% %---------------------------------------------------------------------% %
% % Heat equation with distributed disturbance:
% % PDE         x_{t} = x_{ss} + s*w(t)
% % With BCs    x(s=0) = 0
% %             x_{s}(s=1) = 0
% % And Output  z(t) = int(x(t,s),s,0,1)
% %
% % gamma = 0.2083 (= 5/24, attained at frequency 0).                       % MMP, 11/21/2026
% %---------------------------------------------------------------------% %

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PIETOOLS - PIETOOLS_PDE_Ex_Heat_Eq_with_Distributed_Disturbance
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
% Initial coding DJ - 08/18/2022
% MMP, 11/21/2026: Changes since PIETOOLS 2025: third output st, the
%   recommended settings of the H-infinity gain test; Hinf_gain_form = 'P'
%   set in the base workspace; help revised: the stated gain 0.3333, that
%   of x_t = x_ss + w, corrected to gamma = 5/24 = 0.2083, the gain at
%   frequency 0 with the disturbance s*w (steady state x = s/2 - s^3/6);
%   license added.
% function [PDE_t,PDE_b] = PIETOOLS_PDE_Ex_Heat_Eq_with_Distributed_Disturbance(GUI,params) % MMP, 11/21/2026 (was)

% Determine the location of this example file <-- DO NOT MOVE THE FILE
loc = mfilename('fullpath');
root = fileparts(loc);

% Initialize variables
pvar s

%%% Executive Function:
evalin('base','Hinf_gain = 1;')
evalin('base','Hinf_gain_form = ''P'';');  % form of the gain test          % MMP, 11/21/2026

% Default state size, overridden by the assignment strings in params
ne = 1;
npars = length(params);
if npars~=0
    %%% Specify potential parameters
    for j=1:npars
        eval(params{j});
    end
end


% % % Construct the PDE.
on = eye(ne);     ze = zeros(ne);

%%% Batch input format
PDE_b.n0 = 0;   PDE_b.n1 = 0;   PDE_b.n2 = ne;   PDE_b.nw = ne;   PDE_b.nz = ne;
PDE_b.dom = [0,1];

PDE_b.A2 = on;    PDE_b.B21 = s*on;
PDE_b.Ca1 = on;

PDE_b.B = [on ze ze ze;
           ze ze ze on];


%%% Term-based input format
% Initialize 1D PDE state component.
PDE_t.x{1}.vars = s;   PDE_t.x{1}.dom = [0,1];
% Initialize finite-dimensional inputs and outputs.
PDE_t.w{1}.vars = [];  PDE_t.z{1}.vars = [];

% PDE: x_{t} = x_{ss}
PDE_t.x{1}.term{1}.x = 1;
PDE_t.x{1}.term{1}.D = 2;

% PDE: x_{t} = ... + s * w,
PDE_t.x{1}.term{2}.w = 1;
PDE_t.x{1}.term{2}.C = s*eye(ne);

% Output: z = int_{0}^{1} x(s) ds
PDE_t.z{1}.term{1}.x = 1;
PDE_t.z{1}.term{1}.I{1} = PDE_t.x{1}.dom;

% BC 1: 0 = x(0)
PDE_t.BC{1}.term{1}.x = 1;
PDE_t.BC{1}.term{1}.loc = 0;

% BC 2: 0 = x_{s}(1)
PDE_t.BC{2}.term{1}.x = 1;
PDE_t.BC{2}.term{1}.D = 1;
PDE_t.BC{2}.term{1}.loc = 1;


if GUI
    %%% Associated GUI save file
    app = PIETOOLS_PDE_GUI;
    load(fullfile(root,'PIETOOLS_PDE_Ex_Heat_Eq_with_Distributed_Disturbance_GUI.mat'));
    logval = app.loadData(data);
    if logval
        disp("Failed to load data object. Incorrect structure");
    end
end

% Example-specific settings                                                 % MMP, 11/21/2026
% (stripped tier; 1-D storage without the upper-integral basis operator,
% s-length 0, and theta-lengths 2 (multiplier), 1 (lower integral); slack
% lift degree at its minimum Dmin)
st = lpisettings_sop('stripped');                                           % MMP, 11/21/2026
st.lpi.pos.len_com{1} = 0;                                                  % MMP, 11/21/2026
st.lpi.pos.len_sos{1} = [2;1;1];                                            % MMP, 11/21/2026
st.lpi.ineq.dD = 0;                                                         % MMP, 11/21/2026

end