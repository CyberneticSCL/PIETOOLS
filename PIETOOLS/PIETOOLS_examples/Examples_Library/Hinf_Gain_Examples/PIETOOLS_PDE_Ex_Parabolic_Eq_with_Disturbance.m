function [PDE_t,PDE_b,st] = PIETOOLS_PDE_Ex_Parabolic_Eq_with_Disturbance(GUI,params) % MMP, 11/21/2026
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% PIETOOLS PDE Examples
% [PDE_t,PDE_b,ST] = PIETOOLS_PDE_EX_PARABOLIC_EQ_WITH_DISTURBANCE(GUI,
% PARAMS) returns the parabolic PDE below, with disturbance w and boundary
% output z = x(t,1), in the term-based and in the batch format, and the
% recommended settings ST of its H-infinity gain test. It sets Hinf_gain =
% 1 and Hinf_gain_form = 'P' (storage <T x, R T x>, R >= eppos I) in the
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
% - st:         example-specific settings (lpisettings structure)       % MMP, 11/21/2026
%               for the tests this file names                               % MMP, 11/21/2026
%
% %---------------------------------------------------------------------% %
% % Parabolic PDE example from Shivakumar (2019) (Example 1):
% % PDE         x_{t} = A0(s)*x + A1(s)*x_{s} + A2(s)*x_{ss} + w(t)
% % With BCs    x(s=0) = 0
% %             x_{s}(s=1) = 0
% % And Output  z(t) = x(t,1)
% % 
% % where lam = 4.6, A0 = -0.5*s^3 + 1.3*s^2 - 1.5*s + 0.7 + lam,
% % A1 = 3*s^2 - 2*s, and A2 = s^3 - s^2 + 2.
% %
% % Parameters lam, A0, A1 and A2 can be set.
% % gamma = 23.706 for lam = 4.6, attained at frequency 0.
% %---------------------------------------------------------------------% %

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PIETOOLS - PIETOOLS_PDE_Ex_Parabolic_Eq_with_Disturbance
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
%   recommended settings of the H-infinity gain test; Hinf_gain_form = 'P'
%   set in the base workspace; the warning that SeDuMi requires the
%   veryheavy settings removed; help revised (stated gain 15.147 corrected
%   to 23.706, attained at frequency 0), license added.
% MMP, 11/28/2026: _sop dropped from the names of the executives, LPI
%   and settings functions (2025 names where they exist); names only.
% function [PDE_t,PDE_b] = PIETOOLS_PDE_Ex_Parabolic_Eq_with_Disturbance(GUI,params) % MMP, 11/21/2026 (was)

% Determine the location of this example file <-- DO NOT MOVE THE FILE
loc = mfilename('fullpath');
root = fileparts(loc);

% Initialize variables
pvar s

%%% Executive Function:
evalin('base','Hinf_gain = 1;')
evalin('base','Hinf_gain_form = ''P'';');  % form of the gain test          % MMP, 11/21/2026
%%% evalin('base','Hinf_gain_dual = 1;')

% Default parameters, overridden by the assignment strings in params
lam = 4.6; 
A0 = -0.5*s^3 + 1.3*s^2 - 1.5*s + 0.7 + lam;
A1 = 3*s^2 - 2*s; 
A2 = s^3 - s^2 + 2;
npars = length(params);
if npars~=0
    %%% Specify potential parameters
    for j=1:npars
        eval(params{j});
    end
end


% % % Construct the PDE.
%%% Batch input format
PDE_b.n0 = 0;   PDE_b.n1 = 0;   PDE_b.n2 = 1;   PDE_b.nw = 1;   PDE_b.nz = 1;
PDE_b.dom = [0,1];

PDE_b.A0 = A0;  PDE_b.A1 = A1;   PDE_b.A2 = A2;
PDE_b.B21 = 1;
PDE_b.C10 = [0,1,0,0];

PDE_b.B = [1 0 0 0;
           0 0 0 1];


%%% Term-based input format
% Initialize 1D PDE state component.
PDE_t.x{1}.vars = s;   PDE_t.x{1}.dom = [0,1];
% Initialize finite-dimensional inputs and outputs.
PDE_t.w{1}.vars = [];  PDE_t.z{1}.vars = [];

% PDE: x_{t} = A0 * x + A1 * x_{s} + A2 * x_{ss}
PDE_t.x{1}.term{1}.x = 1;
PDE_t.x{1}.term{1}.D = [0; 1; 2];
PDE_t.x{1}.term{1}.C = [A0, A1, A2];

% PDE: x_{t} = ... + w;
PDE_t.x{1}.term{2}.w = 1;

% Output: z = x(s=1)
PDE_t.z{1}.term{1}.x = 1;
PDE_t.z{1}.term{1}.loc = 1;

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
    load(fullfile(root,'PIETOOLS_PDE_Ex_Parabolic_Eq_with_Disturbance_GUI.mat'));
    logval = app.loadData(data);
    if logval
        disp("Failed to load data object. Incorrect structure");
    end
end
% disp('Warning: With sedumi solver, this example requires at least settings veryheavy.') % MMP, 11/21/2026 (was)
% Example-specific settings                                                 % MMP, 11/21/2026
% (stripped tier; 1-D storage with the basis operators a (multiplier) and
% b (upper integral) only, pattern 'separable')
st = lpisettings('stripped');                                               % MMP, 11/21/2026
st.lpi.pos.patterns{1} = 'separable';                                       % MMP, 11/21/2026

end
% @inproceedings{shivakumar2019computing,
%   title={Computing input-ouput properties of coupled linear pde systems},
%   author={Shivakumar, Sachin and Peet, Matthew M},
%   booktitle={2019 American Control Conference (ACC)},
%   pages={606--613},
%   year={2019},
%   organization={IEEE}
% }