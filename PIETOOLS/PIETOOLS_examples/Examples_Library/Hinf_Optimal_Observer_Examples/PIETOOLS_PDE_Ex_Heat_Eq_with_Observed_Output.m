function [PDE_t,PDE_b,st] = PIETOOLS_PDE_Ex_Heat_Eq_with_Observed_Output(GUI,params) % MMP, 11/21/2026
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% PIETOOLS PDE Examples
% [PDE_t,PDE_b,ST] = PIETOOLS_PDE_EX_HEAT_EQ_WITH_OBSERVED_OUTPUT(GUI,
% PARAMS) returns the heat equation below, with disturbance w, regulated
% output z and measured output y, in the term-based and in the batch
% format, and the recommended settings ST of its H-infinity optimal
% estimator synthesis. It sets Hinf_estimator = 1 in the base workspace.
% With GUI = 1 it also loads the associated GUI save file into a new
% PIETOOLS_PDE_GUI window.
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
% % Example Heat Equation:
% % PDE             x_{t} = x_{ss} + w(t)
% % With BCs        x(s=0) = 0
% %                        x(s=1) = 0
% % And Outputs     z(t) = int(x(t,s),s,0,1) + w(t)
% %                 y(t) = x_{s}(s=1)
% %
% % Answer: 1.0045.
% %---------------------------------------------------------------------% %

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PIETOOLS - PIETOOLS_PDE_Ex_Heat_Eq_with_Observed_Output
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
% MMP, 11/21/2026: Changes since PIETOOLS 2025: third output st, the
%   recommended settings of the H-infinity optimal estimator synthesis;
%   help revised, license added.
% function [PDE_t,PDE_b] = PIETOOLS_PDE_Ex_Heat_Eq_with_Observed_Output(GUI,params) % MMP, 11/21/2026 (was)

% Determine the location of this example file <-- DO NOT MOVE THE FILE
loc = mfilename('fullpath');
root = fileparts(loc);

% Initialize variables
pvar s

%%% Executive Function:
evalin('base','Hinf_estimator = 1;')

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
PDE_b.nw = ne;   PDE_b.ny = ne;   PDE_b.nz = ne;   PDE_b.nx = 0;
PDE_b.n0 = 0;    PDE_b.n1 = 0;    PDE_b.n2 = ne;
PDE_b.dom = [0,1];

PDE_b.A2 = on; PDE_b.Ca1 = on; PDE_b.B21 = on;   PDE_b.D11 = on;
PDE_b.C20 = [ze ze ze on];  

PDE_b.B = [on ze ze ze;
           ze on ze ze];


%%% Term-based input format
% Initialize 1D PDE state component.
PDE_t.x{1}.vars = s;   PDE_t.x{1}.dom = [0,1];
% Initialize finite-dimensional inputs and outputs.
PDE_t.z{1}.vars = [];  PDE_t.w{1}.vars = [];
PDE_t.y{1}.vars = [];

% PDE: x_{t} = x_{ss}
PDE_t.x{1}.term{1}.x = 1;
PDE_t.x{1}.term{1}.D = 2;
PDE_t.x{1}.term{1}.C = eye(ne);

% PDE: x_{t} = ... + w
PDE_t.x{1}.term{2}.w = 1;
PDE_t.x{1}.term{2}.C = eye(ne);

% Output: z = int_{0}^{1} x(s) ds
PDE_t.z{1}.term{1}.x = 1;
PDE_t.z{1}.term{1}.C = eye(ne);
PDE_t.z{1}.term{1}.I{1} = PDE_t.x{1}.dom;

% Output: z = ... + w
PDE_t.z{1}.term{2}.w  =1;

% Output: y = x_{s}(s=1)
PDE_t.y{1}.term{1}.x = 1;
PDE_t.y{1}.term{1}.D = 1;
PDE_t.y{1}.term{1}.loc = 1;

% BC 1: 0 = x(0)
PDE_t.BC{1}.term{1}.x = 1;
PDE_t.BC{1}.term{1}.loc = 0;

% BC 2: 0 = x(1)
PDE_t.BC{2}.term{1}.x = 1;
PDE_t.BC{2}.term{1}.loc = 1;


if GUI
    %%% Associated GUI save file
    app = PIETOOLS_PDE_GUI;
    load(fullfile(root,'PIETOOLS_PDE_Ex_Heat_Eq_with_Observed_Output_GUI.mat'));
    logval = app.loadData(data);
    if logval
        disp("Failed to load data object. Incorrect structure");
    end
end

% Example-specific settings                                                 % MMP, 11/21/2026
% (light tier; 1-D storage with theta-lengths 2, 1, 1 (multiplier, upper,
% lower integral); slack lift degree at its minimum Dmin; free operators
% of degree 1 in the multiplier and 2 in the integral kernels; slack weight
% raise of the synthesis LPI 0 in 1-D)
st = lpisettings_sop('light');                                              % MMP, 11/21/2026
st.lpi.pos.len_sos{1} = [2;1;1];                                            % MMP, 11/21/2026
st.lpi.ineq.dD = 0;                                                         % MMP, 11/21/2026
st.lpi.var.deg{1} = [1 2 2];                                                % MMP, 11/21/2026
st.sop.dw_P = [0 2 1];                                                      % MMP, 11/21/2026

end