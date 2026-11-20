function [PDE_t,PDE_b,st] = PIETOOLS_PDE_Ex_Transport_Eq_with_Observed_Output(GUI,params)  % MMP, 11/19/2026
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% PIETOOLS PDE Examples
% [PDE_t,PDE_b,ST] = PIETOOLS_PDE_EX_TRANSPORT_EQ_WITH_OBSERVED_OUTPUT(GUI,
% PARAMS) returns the transport equation below, with disturbance w,
% regulated output z and measured output y, in the term-based and in the
% batch format, and the recommended settings ST of its H-infinity optimal
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
% - st:         example-specific settings (lpisettings_sop structure)       % MMP, 11/19/2026
%               for the tests this file names                               % MMP, 11/19/2026
%
% %---------------------------------------------------------------------% %
% % Example 1D transport:
% % PDE             x_{t} = x_{s} + w(t)
% % with BC         x(s=1) = 0
% % And Outputs     z(t) = int(x(t,s),s,0,1) + w(t)
% %                 y(t) = int(x(t,s),s,0,1)
% %
% % Answer: 1.0012.
% %---------------------------------------------------------------------% %

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PIETOOLS - PIETOOLS_PDE_Ex_Transport_Eq_with_Observed_Output
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
% MMP, 11/19/2026: Output st, the example-specific settings, and the form
%                   of its tests named in the workspace (stability_form,
%                   Hinf_gain_form, H2_norm_form), from the search of
%                   11/19/2026 (maintainer: no default settings);
% function [PDE_t,PDE_b] = PIETOOLS_PDE_Ex_Transport_Eq_with_Observed_Output(GUI,params)  % MMP, 11/19/2026 (was)
% MMP, 11/19/2026: Release comment pass (comment-only).

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
PDE_b.n0 = 0;    PDE_b.n1 = ne;   PDE_b.n2 = 0;
PDE_b.dom = [0,1];

PDE_b.A1 = on;
PDE_b.Ca1 = on;   PDE_b.Ca2 = on;
PDE_b.B21 = on;   PDE_b.D11 = on;

PDE_b.B = [ze on];


%%% Term-based input format
% Initialize 1D PDE state component.
PDE_t.x{1}.vars = s;   PDE_t.x{1}.dom = [0,1];
% Initialize finite-dimensional inputs and outputs.
PDE_t.z{1}.vars = [];  PDE_t.w{1}.vars = [];
PDE_t.y{1}.vars = [];

% PDE: x_{t} = x_{s}
PDE_t.x{1}.term{1}.x = 1;
PDE_t.x{1}.term{1}.D = 1;

% PDE: x_{t} = ... + w;
PDE_t.x{1}.term{2}.w = 1;

% Output: z = int_{0}^{1} x(s) ds
PDE_t.z{1}.term{1}.x = 1;
PDE_t.z{1}.term{1}.I{1} = PDE_t.x{1}.dom;

% Output: z = ... + w
PDE_t.z{1}.term{2}.w  = 1;
PDE_t.z{1}.term{2}.C  = 1;

% Output: y = int_{0}^{1} x(s) ds
PDE_t.y{1}.term{1}.x = 1;
PDE_t.y{1}.term{1}.I{1} = PDE_t.x{1}.dom;

% BC 1: 0 = x(1)
PDE_t.BC{1}.term{1}.x = 1;
PDE_t.BC{1}.term{1}.loc = 1;


if GUI
    %%% Associated GUI save file
    app = PIETOOLS_PDE_GUI;
    load(fullfile(root,'PIETOOLS_PDE_Ex_Transport_Eq_with_Observed_Output_GUI.mat'));
    logval = app.loadData(data);
    if logval
        disp("Failed to load data object. Incorrect structure");
    end
end

% Example-specific settings                                                 % MMP, 11/19/2026
% (extreme tier; 1-D storage from the multiplier only, theta-length 1;
% free operators of degree 0 in the multiplier and 0, 2 in the integral
% kernels (output, integration variable); no slack weight raise of the
% synthesis LPI)
st = lpisettings_sop('extreme');                                            % MMP, 11/19/2026
st.lpi.pos.len_sos{1} = [1;0;0];                                            % MMP, 11/19/2026
st.lpi.var.deg{1} = [0 2 0];                                                % MMP, 11/19/2026
st.sop.dw_P = 0;                                                            % MMP, 11/19/2026

end