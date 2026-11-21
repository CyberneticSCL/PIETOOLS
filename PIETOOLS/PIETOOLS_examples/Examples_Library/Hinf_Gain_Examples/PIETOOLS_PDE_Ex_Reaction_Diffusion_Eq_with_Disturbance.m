function [PDE_t,PDE_b,st] = PIETOOLS_PDE_Ex_Reaction_Diffusion_Eq_with_Disturbance(GUI,params) % MMP, 11/21/2026
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% PIETOOLS PDE Examples
% [PDE_t,PDE_b,ST] = PIETOOLS_PDE_EX_REACTION_DIFFUSION_EQ_WITH_DISTURBANCE(
% GUI,PARAMS) returns the reaction-diffusion PDE below, with disturbance w
% and integral output z, in the term-based and in the batch format, and
% the recommended settings ST of its H-infinity gain test. It sets
% Hinf_gain = 1 and Hinf_gain_form = 'Q' (generalized storage T'Q) in the
% base workspace. With GUI = 1 it also loads the associated GUI save file
% into a new PIETOOLS_PDE_GUI window.
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
% % Diffusion-Reaction equation from Shivakumar (Example 3):
% % PDE         xi_{t} = lam*xi + sum(xk_{ss},k=1,i) + w(t)
% % With BCs    x(s=0) = 0
% %             x(s=1) = 0
% % And Output  z(t) = int(x(t,s),s,0,1)
% %
% % Parameter lam and the state size ne can be set.
% % gamma = 8.2141 for lam = (1-1e-2)*pi^2 and ne = 1, attained at
% % frequency 0: gamma = (2*tan(sqrt(lam)/2)/sqrt(lam) - 1)/lam.
% %---------------------------------------------------------------------% %

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PIETOOLS - PIETOOLS_PDE_Ex_Reaction_Diffusion_Eq_with_Disturbance
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
%   recommended settings of the H-infinity gain test; Hinf_gain_form = 'Q'
%   set in the base workspace; the SeDuMi warning states that the settings
%   st do not certify the gain with SeDuMi and the light tier does; help
%   revised (stated gain 8.1069 corrected to 8.2141, with its formula),
%   license added.
% function [PDE_t,PDE_b] = PIETOOLS_PDE_Ex_Reaction_Diffusion_Eq_with_Disturbance(GUI,params) % MMP, 11/21/2026 (was)

% Determine the location of this example file <-- DO NOT MOVE THE FILE
loc = mfilename('fullpath');
root = fileparts(loc);

% Initialize variables
pvar s

%%% Executive Function:
evalin('base','Hinf_gain = 1;')
evalin('base','Hinf_gain_form = ''Q'';');  % form of the gain test          % MMP, 11/21/2026
%%% evalin('base','Hinf_gain_dual = 1;')

% Default parameters, overridden by the assignment strings in params
lam = (1-1e-2)*pi^2;    ne = 1;
npars = length(params);
if npars~=0
    %%% Specify potential parameters
    for j=1:npars
        eval(params{j});
    end
end


% % % Construct the PDE.
on = eye(ne);     ze = zeros(ne);     lt = tril(ones(ne));

%%% Batch input format
PDE_b.n0 = 0;   PDE_b.n1 = 0;   PDE_b.n2 = ne;   PDE_b.nw = 1;   PDE_b.nz=1;
PDE_b.dom = [0,1];
PDE_b.A0 = lam*on;   PDE_b.A2 = lt;
PDE_b.B21 = ones(ne,1);
PDE_b.Ca1 = ones(1,ne);

PDE_b.B = [on ze ze ze;
           ze on ze ze];


%%% Term-based input format
% Initialize 1D PDE state component.
PDE_t.x{1}.vars = s;   PDE_t.x{1}.dom = [0,1];
PDE_t.x{1}.size = ne;
% Initialize finite-dimensional inputs and outputs.
PDE_t.w{1}.vars = [];  PDE_t.z{1}.vars = [];

% PDE: x_{t} = lam * xi
PDE_t.x{1}.term{1}.x = 1;
PDE_t.x{1}.term{1}.C = lam*eye(ne);

% PDE: x_{t} = ... + sum_{k=1}^{i}(xk_{ss})
PDE_t.x{1}.term{2}.x = 1;
PDE_t.x{1}.term{2}.D = 2;
PDE_t.x{1}.term{2}.C = lt; % lower-triangular matrix

% PDE: x_{t} = ... + w
PDE_t.x{1}.term{3}.w = 1;

% Output: z = int_{0}^{1} x(s) ds
PDE_t.z{1}.term{1}.x = 1;
PDE_t.z{1}.term{1}.I{1} = PDE_t.x{1}.dom;

% BC 1: 0 = x(0)
PDE_t.BC{1}.term{1}.x = 1;
PDE_t.BC{1}.term{1}.loc = 0;

% BC 2: 0 = x(1)
PDE_t.BC{2}.term{1}.x = 1;
PDE_t.BC{2}.term{1}.loc = 1;


if GUI
    %%% Associated GUI save file
    app = PIETOOLS_PDE_GUI;
    load(fullfile(root,'PIETOOLS_PDE_Ex_Reaction_Diffusion_Eq_with_Disturbance_GUI.mat'));
    logval = app.loadData(data);
    if logval
        disp("Failed to load data object. Incorrect structure");
    end
end
% disp('Warning: With sedumi solver, this example requires at least settings veryheavy.') % MMP, 11/21/2026 (was)
disp('Warning: with SeDuMi, the settings st do not certify the gain; lpisettings_sop(''light'') does.') % MMP, 11/21/2026
% Example-specific settings                                                 % MMP, 11/21/2026
% (light tier; slack lift degree at its minimum Dmin; slack integral bases
% capped at total degree theta+s <= Dmax-1)
st = lpisettings_sop('light');                                              % MMP, 11/21/2026
st.lpi.ineq.dD = 0;                                                         % MMP, 11/21/2026
st.lpi.ineq.cap_inc = 'auto';                                               % MMP, 11/21/2026

end
% @inproceedings{shivakumar2019computing,
%   title={Computing input-ouput properties of coupled linear pde systems},
%   author={Shivakumar, Sachin and Peet, Matthew M},
%   booktitle={2019 American Control Conference (ACC)},
%   pages={606--613},
%   year={2019},
%   organization={IEEE}
% }