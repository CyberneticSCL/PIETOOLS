function [PDE_t,PDE_b,st] = PIETOOLS_PDE_Ex_Timoshenko_Beam_2(GUI,params)   % MMP, 11/19/2026
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% PIETOOLS PDE Examples
% [PDE_t,PDE_b,ST] = PIETOOLS_PDE_EX_TIMOSHENKO_BEAM_2(GUI,PARAMS) returns
% the Timoshenko beam below in hyperbolic/diffusive form, all coefficients
% 1, in the term-based and in the batch format, and the recommended
% settings ST of its dual stability test. It sets stability_dual = 1 and
% stability_form = 'Q' (the dual form 'Qd') in the base workspace. With
% GUI = 1 it also loads the associated GUI save file into a new
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
%               for the dual stability test                                 % MMP, 11/19/2026
%
% %---------------------------------------------------------------------% %
% % Timoschenko Beam equation Example (hyperbolic/diffusive) from Peet:
% % PDE        r*aa * w_{tt} = k*aa*g * (-phi_{s} + w_{ss})
% %     r*II * phi_{tt} = E*II * phi_{ss} + k*aa*g * (w_{s} - phi) -c \phi_t
% % with BCs   phi(s=0) = 0
% %            w(s=0) = 0
% %            phi_{s}(s=1) = 0
% %            w_{s}(s=1) - phi(s=1) = 0
% %
% % We now use states x1 = w_{t}, x2 = w_{s}, x3 = phi_{t}, x4 = phi.
% % We also assume all coefficients are 1. Then, our system becomes:
% % PDE        x1_{t} = x2_{s} - x4_{s}
% %            x2_{t} = x1_{s}
% %            x3_{t} = x2 - x4-c x_3                                       % MMP, 11/19/2026 (was)
% %            x3_{t} = x2 - x4 + x4_{ss} - c x_3                           % MMP, 11/19/2026
% %            x4_{t} = x3
% % with BCs   x4(0) = 0, x4_{s}(1) = 0, x3(0) = 0, x1(0) = 0, x2(1)-x4(1) = 0
% %
% % Unstable.                                                               % MMP, 11/19/2026 (was)
% % With x4_{ss}: at c = 0 the conservative beam, stable, not               % MMP, 11/19/2026
% % exponentially stable (Peet, Automatica 2021, Sec. 11.1); stable         % MMP, 11/19/2026
% % for c > 0 (Peet et al. 2019, Sec. 8.1).                                 % MMP, 11/19/2026
% %---------------------------------------------------------------------% %

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PIETOOLS - PIETOOLS_PDE_Ex_Timoshenko_Beam_2
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
% MMP, 11/19/2026: the bending term x4_{ss} (E*II*phi_{ss} above) restored
%                   in x3_t, in both formats: the term-based format never
%                   had it and 72de7c6b (12/07/2022) set the batch A2 to
%                   zero, so the coded system grew linearly at c = 0.
%                   Source: Peet et al. 2019 (below), Sec. 8.1, A2 =
%                   [0;0;1;0] (maintainer, 11/19/2026: "Apply all four");
% MMP, 11/19/2026: output st, the example-specific settings of the dual
%                   stability test, and its form ('Q'), from the search of
%                   11/19/2026 (maintainer: no default settings);
% function [PDE_t,PDE_b] = PIETOOLS_PDE_Ex_Timoshenko_Beam_2(GUI,params)    % MMP, 11/19/2026 (was)
% MMP, 11/19/2026: Release comment pass (comment-only).

% Determine the location of this example file <-- DO NOT MOVE THE FILE
loc = mfilename('fullpath');
root = fileparts(loc);

% Initialize variables
pvar s

%%% Executive Function:
evalin('base','stability_dual= 1;');
evalin('base','stability_form = ''Q'';');  % form of the stability test     % MMP, 11/19/2026


% Specify the parameters
c=0;
npars = length(params);
if npars~=0
    %%% Specify potential parameters
    for j=1:npars
        eval(params{j});
    end
end

% Construct the PDE.
% Batch input format
PDE_b.n0 = 0;   PDE_b.n1 = 3;   PDE_b.n2=1;
PDE_b.dom = [0,1];

PDE_b.A0 = [0 0 0 0;
    0 0 0 0;
    0 1 -c -1;
    0 0 1 0];
PDE_b.A1 = [0 1 0 -1;
    1 0 0 0;
    0 0 0 0;
    0 0 0 0];
% PDE_b.A2 = [0;0;0;0];                                                     % MMP, 11/19/2026 (was)
PDE_b.A2 = [0;0;1;0];   % x4_{ss} in x3_t, the bending term                 % MMP, 11/19/2026

PDE_b.B = [1 0 0 0 0 0 0 0 0 0;
           0 0 1 0 0 0 0 0 0 0;
           0 0 0 0 0 0 1 0 0 0;
           0 0 0 0 0 0 0 0 0 1;
           0 0 0 0 1 0 0 -1 0 0];


%%% Term-based input format
PDE_t.vars = s;    PDE_t.dom = [0,1];

% PDE 1: x1_{t} = x2_{s} - x4_{s}
PDE_t.x{1}.term{1}.x = [2;4];
PDE_t.x{1}.term{1}.D = 1;
PDE_t.x{1}.term{1}.C = [1,-1];

% PDE 2: x2_{t} = x1_{s}
PDE_t.x{2}.term{1}.x = 1;
PDE_t.x{2}.term{1}.D = 1;

% PDE 3: x3_{t} = x2 - x4-c x_3                                             % MMP, 11/19/2026 (was)
% PDE 3: x3_{t} = x2 - x4 + x4_{ss} - c x_3                                 % MMP, 11/19/2026
PDE_t.x{3}.term{1}.x = [2;4;3];
PDE_t.x{3}.term{1}.C = [1,-1,-c];
PDE_t.x{3}.term{2}.x = 4;                                                   % MMP, 11/19/2026
PDE_t.x{3}.term{2}.D = 2;                                                   % MMP, 11/19/2026

% PDE 4: x4_{t} = x3
PDE_t.x{4}.term{1}.x = 3;

% BC 1: x4(0) = 0
PDE_t.BC{1}.term{1}.x = 4;
PDE_t.BC{1}.term{1}.loc = 0;

% BC 2: x4_{s}(1) = 0
PDE_t.BC{2}.term{1}.x = 4;
PDE_t.BC{2}.term{1}.D = 1;
PDE_t.BC{2}.term{1}.loc = 1;

% BC 3: x3(0) = 0
PDE_t.BC{3}.term{1}.x = 3;
PDE_t.BC{3}.term{1}.loc = 0;

% BC 4: x1(0) = 0
PDE_t.BC{4}.term{1}.x = 1;
PDE_t.BC{4}.term{1}.loc = 0;

% BC 5: x2(1)-x4(1) = 0
PDE_t.BC{5}.term{1}.x = [2;4];
PDE_t.BC{5}.term{1}.loc = 1;
PDE_t.BC{5}.term{1}.C = [1,-1];

% % % Alternatively, group together x1 x2 and x3
%  PDE_t.vars = s;    PDE_t.dom = [0,1];
%
%  % PDE: x3_{t} = x2
%  PDE_t.x{1}.term{1}.x = 1;
%  PDE_t.x{1}.term{1}.C = [0 0 0; 0 0 0; 0 1 0];
%
%  % PDE: x1_{t} = x2_{s},   x2_{t} = x1_{s}
%  PDE_t.x{1}.term{2}.x = 1;
%  PDE_t.x{1}.term{2}.D = 1;
%  PDE_t.x{1}.term{2}.C = [0 1 0; 1 0 0; 0 0 0];
%
%  % PDE: x3_{t} = ... - x4
%  PDE_t.x{1}.term{3}.x = 2;
%  PDE_t.x{1}.term{3}.C = [0; 0; -1];
%
%  % PDE: x1_{t} = ... - x4_{s}
%  PDE_t.x{1}.term{4}.x = 2;
%  PDE_t.x{1}.term{4}.D = 1;
%  PDE_t.x{1}.term{4}.C = [-1; 0; 0];
%
%  % PDE: x3_{t} = ... + x4_{ss}
%  PDE_t.x{1}.term{5}.x = 2;
%  PDE_t.x{1}.term{5}.D = 2;
%  PDE_t.x{1}.term{5}.C = [0; 0; 1];
%
%  % PDE: x4_{t} = x3
%  PDE_t.x{2}.term{1}.x = 1;
%  PDE_t.x{2}.term{1}.C = [0 0 1];
%
%  % BC 1: 0 = x1(0)
%  %       0 = x3(0)
%  PDE_t.BC{1}.term{1}.x = 1;
%  PDE_t.BC{1}.term{1}.loc = 0;
%  PDE_t.BC{1}.term{1}.C = [1,0,0;0,0,1];
%
%  % BC 2: 0 = x4(0)
%  PDE_t.BC{2}.term{1}.x = 2;
%  PDE_t.BC{2}.term{1}.loc = 0;
%
%  % BC 3: 0 = x4_{s}(1)
%  PDE_t.BC{3}.term{1}.x = 2;
%  PDE_t.BC{3}.term{1}.D = 1;
%  PDE_t.BC{3}.term{1}.loc = 1;
%
%  % BC 4: 0 = x2(1)
%  PDE_t.BC{4}.term{1}.x = 1;
%  PDE_t.BC{4}.term{1}.loc = 1;
%  PDE_t.BC{4}.term{1}.C = [0,1,0];
%
%  % BC 4: 0 = ... - x4(1)
%  PDE_t.BC{4}.term{2}.x = 2;
%  PDE_t.BC{4}.term{2}.loc = 1;
%  PDE_t.BC{4}.term{2}.C = -1;


% Example-specific settings                                                 % MMP, 11/19/2026
% (veryheavy tier; slack lift degree at its minimum Dmin)
st = lpisettings_sop('veryheavy');                                          % MMP, 11/19/2026
st.lpi.ineq.dD = 0;                                                         % MMP, 11/19/2026

if GUI
    %%% Associated GUI save file
    app = PIETOOLS_PDE_GUI;
    load(fullfile(root,'PIETOOLS_PDE_Ex_Timoshenko_Beam_2_GUI.mat'));
    logval = app.loadData(data);
    if logval
        disp("Failed to load data object. Incorrect structure");
    end
end

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