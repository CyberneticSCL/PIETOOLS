function [PDE_t,PDE_b,st] = PIETOOLS_PDE_Ex_Wave_Eq_Datko_Boundary_Damped_1(GUI,params) % DB, MMP, 11/21/2026
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% PIETOOLS PDE Examples
% [PDE_t,PDE_b,ST] = PIETOOLS_PDE_EX_WAVE_EQ_DATKO_BOUNDARY_DAMPED_1(GUI,
% PARAMS) returns the damped wave equation below (Datko, Test 7.5c) in the
% states x1 = u_t, x2 = u, in the term-based and in the batch format, and
% the recommended settings ST of its dual stability test. It sets
% stability_dual = 1 and stability_form = 'P' (the dual form 'Pd') in the
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
% - st:         example-specific settings (lpisettings structure)       % DB, MMP, 11/21/2026
%               for the tests this file names                               % DB, MMP, 11/21/2026
%
% %---------------------------------------------------------------------% %
% % Test 7.5c - Datko Boundary-Damped Wave equation (see reference below):
% % PDE        u_{tt} + 2*ad*u_{t} = -ad^2*u + u_{ss} 
% % with BCs   u(s=0) = 0
% %            u_{s}(s=1) = -k*u_{t}(s=1)
% %
% % We use states x1 = u_{t}, x2 = u.
% % Then x1(0) = 0, x2(0) = 0,  k*x1(1) + x2_{s}(1) = 0.
% % Parameters k and ad can be set.
% %---------------------------------------------------------------------% %

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PIETOOLS - PIETOOLS_PDE_Ex_Wave_Eq_Datko_Boundary_Damped_1
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
% DB, MMP, 11/21/2026: Changes since PIETOOLS 2025: the term-based format
%   does not set the field PDE_t.PDE (an assignment of an older input
%   format, commented out); new output st, the recommended settings of the
%   dual stability test; stability_form = 'P' set in the base workspace;
%   help and comments revised.
% MMP, 11/28/2026: _sop dropped from the names of the executives, LPI
%   and settings functions (2025 names where they exist); names only.
% function [PDE_t,PDE_b] = PIETOOLS_PDE_Ex_Wave_Eq_Datko_Boundary_Damped_1(GUI,params) % DB, MMP, 11/21/2026 (was)

% Determine the location of this example file <-- DO NOT MOVE THE FILE
loc = mfilename('fullpath');
root = fileparts(loc);

% Initialize variables
pvar s

%%% Executive Function:
evalin('base','stability_dual = 1;');
evalin('base','stability_form = ''P'';');  % form of the stability test     % DB, MMP, 11/21/2026

% Default parameters, overridden by the assignment strings in params
k = 1;    ad = 1;
npars = length(params);
if npars~=0
    %%% Specify potential parameters
    for j=1:npars
        eval(params{j});
    end
end


% % % Construct the PDE.
%%% Batch input format
PDE_b.n0 = 0;   PDE_b.n1 = 1;   PDE_b.n2 = 1;
PDE_b.dom = [0,1];

PDE_b.A0 = [-2*ad, -ad^2; 1, 0];
PDE_b.A2 = [1; 0];

PDE_b.B = [0 0 1 0 0 0;
    1 0 0 0 0 0;
    0 k 0 0 0 1];


%%% Term-based input format
PDE_t.vars = s;        PDE_t.dom = [0,1];

% PDE: x1_{t} = -2*ad * x1 - ad^2 * x2
PDE_t.x{1}.term{1}.x = [1;2];
PDE_t.x{1}.term{1}.C = [-2*ad, -ad^2];

% PDE: x1_{t} = ... + x2_{ss}
%PDE_t.PDE.A{4}.coeff = 1; PDE_t.PDE.A{4}.Lstate = 1; PDE_t.PDE.A{4}.Rstate = 2; PDE_t.PDE.A{4}.D = 2;
PDE_t.x{1}.term{2}.x = 2;
PDE_t.x{1}.term{2}.D = 2;

% PDE: x2_{t} = x1
PDE_t.x{2}.term{1}.x = 1;

% BC 1: 0 = x1(0)
PDE_t.BC{1}.term{1}.x = 1;
PDE_t.BC{1}.term{1}.loc = 0;

% BC 2: 0 = x2(0)
PDE_t.BC{2}.term{1}.x = 2;
PDE_t.BC{2}.term{1}.loc = 0;

% BC 3: 0 = k * x1(1) + x2_{s}(1)
PDE_t.BC{3}.term{1}.x = [1;2];
PDE_t.BC{3}.term{1}.loc = 1;
PDE_t.BC{3}.term{1}.D = [0;1];
PDE_t.BC{3}.term{1}.C = [k, 1];


if GUI
    %%% Associated GUI save file
    app = PIETOOLS_PDE_GUI;
    load(fullfile(root,'PIETOOLS_PDE_Ex_Wave_Eq_Datko_Boundary_Damped_1_GUI.mat'));
    logval = app.loadData(data);
    if logval
        disp("Failed to load data object. Incorrect structure");
    end
end

% Example-specific settings                                                 % DB, MMP, 11/21/2026
% (stripped tier; 1-D storage without the upper-integral basis operator,
% s-length 0, and theta-lengths 1 (multiplier), 1 (lower integral); slack
% lift degree at its minimum Dmin)
st = lpisettings('stripped');                                               % DB, MMP, 11/21/2026
st.lpi.pos.len_com{1} = 0;                                                  % DB, MMP, 11/21/2026
st.lpi.pos.len_sos{1} = [1;2;1];                                            % DB, MMP, 11/21/2026
st.lpi.ineq.dD = 0;                                                         % DB, MMP, 11/21/2026

end
% @article{datko1986example,
%   title={An example on the effect of time delays in boundary feedback stabilization of wave equations},
%   author={Datko, Richard and Lagnese, John and Polis, MP},
%   journal={SIAM Journal on Control and Optimization},
%   volume={24},
%   number={1},
%   pages={152--156},
%   year={1986},
%   publisher={SIAM}
% }