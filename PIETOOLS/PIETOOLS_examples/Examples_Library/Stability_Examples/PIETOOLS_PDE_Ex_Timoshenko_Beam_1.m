function [PDE_t,PDE_b,st] = PIETOOLS_PDE_Ex_Timoshenko_Beam_1(GUI,params)   % DB, MMP, 11/21/2026
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% PIETOOLS PDE Examples
% [PDE_t,PDE_b,ST] = PIETOOLS_PDE_EX_TIMOSHENKO_BEAM_1(GUI,PARAMS) returns
% the Timoshenko beam below in hyperbolic first-order form, in the
% term-based and in the batch format, and the recommended settings ST of
% its stability test. It sets stability = 1 and stability_form = 'P'
% (storage <T x, P T x>, P >= eppos I) in the base workspace. With GUI = 1
% it also loads the associated GUI save file into a new PIETOOLS_PDE_GUI
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
% - st:         example-specific settings (lpisettings structure)       % DB, MMP, 11/21/2026
%               for the tests this file names                               % DB, MMP, 11/21/2026
%
% %---------------------------------------------------------------------% %
% % Timoschenko beam equation (hyperbolic) from Peet (Example 8.1.0.2):
% % PDE   r*aa * w_{tt} = k*aa*g * (-phi_{s} + w_{ss})
%       r*II * phi_{tt} = E*II * phi_{ss} + k*aa*g * (w_{s} - phi)-cphi_{t}
% % with BCs   phi(s=0) = 0
% %            w(s=0) = 0 
% %            phi_{s}(s=1) = 0 
% %            w_{s}(s=1) - phi(s=1) = 0 
%
% % We use states x1 = w_{t}, x2 = k*aa*g * (w_{s}-phi), x3 = phi_{t}, x4 = E*II * phi_{s}.
% % Then, our system becomes:
% % PDE        x1_{t} = (1/r/aa) * x2_{s}
% %            x2_{t} = k*aa*g * x1_{s} - k*aa*g * x3
% %            x3_{t} = (1/r/II) * x2 + (1/r/II) * x4_{s} -c x_3
% %            x4_{t} = E*II * x3_{s}
% % with BCs   x1(0) = 0, x3(0) = 0, x4(1) = 0, x2(1) = 0
% % 
% % Parameters k, aa, II, g, E, r and the damping c (default 0) can be set.
% %---------------------------------------------------------------------% %

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PIETOOLS - PIETOOLS_PDE_Ex_Timoshenko_Beam_1
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
% DB, MMP, 11/21/2026: Changes since PIETOOLS 2025: new parameter c
%   (default 0), the damping term -c*x3 = -c*phi_{t} in the dynamics of x3,
%   in both formats; new output st, the recommended settings of the
%   stability test; stability_form = 'P' set in the base workspace; help
%   and comments revised.
% MMP, 11/28/2026: _sop dropped from the names of the executives, LPI
%   and settings functions (2025 names where they exist); names only.
% function [PDE_t,PDE_b] = PIETOOLS_PDE_Ex_Timoshenko_Beam_1(GUI,params)    % DB, MMP, 11/21/2026 (was)

% Determine the location of this example file <-- DO NOT MOVE THE FILE
loc = mfilename('fullpath');
root = fileparts(loc);

% Initialize variables
pvar s

%%% Executive Function:
evalin('base','stability = 1;');
evalin('base','stability_form = ''P'';');  % form of the stability test     % DB, MMP, 11/21/2026

% Specify the parameters
k = 1;    aa = 1;    II = 1;    g = 1;    E = 1;   r = 1;  c=0;%5-30kPa
npars = length(params);
if npars~=0
    %%% Specify potential parameters
    for j=1:npars
        eval(params{j});
    end
end


% Construct the PDE.
%%% Batch input format
PDE_b.n0 = 0;   PDE_b.n1 = 4;   PDE_b.n2 = 0;
PDE_b.dom = [0,1];

PDE_b.A0 = [0 0 0 0;
    0 0 -k*aa*g 0;
    0 1/r/II -c 0;
    0 0 0 0];
PDE_b.A1 = [0 1/r/aa 0 0;
    k*aa*g 0 0 0;
    0 0 0 1/r/II;
    0 0 E*II 0];

PDE_b.B = [ 1 0 0 0 0 0 0 0;
    0 0 1 0 0 0 0 0;
    0 0 0 0 0 0 0 1;
    0 0 0 0 0 1 0 0];


%%% Term-based input format
PDE_t.x{1}.vars = s;    PDE_t.x{1}.dom = [0,1];

% PDE: x2_{t} = -k*aa*g * x3,   x3_{t} = (1/r/II) * x2-c x_3
PDE_t.x{1}.term{1}.C = [0, 0,      0,       0;
    0, 0,      -k*aa*g, 0;
    0, 1/r/II, -c,       0;
    0, 0,      0,       0];

% PDE: x1_{t} = (1/r/aa) * x2_{s},        x2_{t} = ... + k*aa*g * x1_{s}
%      x3_{t} = ... + (1/r/II) * x4_{s}   x4_{t} = E*II * x3_{s}
PDE_t.x{1}.term{2}.D = 1;
PDE_t.x{1}.term{2}.C = [0,      1/r/aa, 0,    0;
    k*aa*g, 0,      0,    0;
    0,      0,      0,    1/r/II;
    0,      0,      E*II, 0];

% BCs: 0 = x1(0),   0 = x3(0)
PDE_t.BC{1}.term{1}.loc = 0;
PDE_t.BC{1}.term{1}.C = [1,0,0,0;0,0,1,0;0,0,0,0;0,0,0,0];

% BCs: 0 = x4(1),   0 = x2(1)
PDE_t.BC{1}.term{2}.loc = 1;
PDE_t.BC{1}.term{2}.C = [0,0,0,0;0,0,0,0;0,0,0,1;0,1,0,0];


if GUI
    %%% Associated GUI save file
    app = PIETOOLS_PDE_GUI;
    load(fullfile(root,'PIETOOLS_PDE_Ex_Timoshenko_Beam_1_GUI.mat'));
    logval = app.loadData(data);
    if logval
        disp("Failed to load data object. Incorrect structure");
    end
end

% Example-specific settings                                                 % DB, MMP, 11/21/2026
% (stripped tier; 1-D storage from the multiplier only, theta-length 1,
% both integral s-lengths 0; slack lift degree at its minimum Dmin)
st = lpisettings('stripped');                                               % DB, MMP, 11/21/2026
st.lpi.pos.len_inc{1} = 0;                                                  % DB, MMP, 11/21/2026
st.lpi.pos.len_com{1} = 0;                                                  % DB, MMP, 11/21/2026
st.lpi.pos.len_sos{1} = [1;2;2];                                            % DB, MMP, 11/21/2026
st.lpi.ineq.dD = 0;                                                         % DB, MMP, 11/21/2026

end
% @article{peet2019discussion,
%   title={Discussion paper: A new mathematical framework for representation and analysis of coupled pdes},
%   author={Peet, Matthew M and Shivakumar, Sachin and Das, Amritam and Weiland, Seip},
%   journal={IFAC-PapersOnLine},
%   volume={52},
%   number={2},
%   pages={132--137},
%   year={2019},
%   publisher={Elsevier}
% }