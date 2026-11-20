function [PDE_t,PDE_b,st] = PIETOOLS_PDE_Ex_Wave_Eq_Boundary_Damped(GUI,params)  % MMP, 11/19/2026
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% PIETOOLS PDE Examples
% [PDE_t,PDE_b,ST] = PIETOOLS_PDE_EX_WAVE_EQ_BOUNDARY_DAMPED(GUI,PARAMS)
% returns the boundary-damped wave equation below in the states
% x1 = u_s, x2 = u_t, in the term-based and in the batch format, and the
% recommended settings ST of its stability test. It sets stability = 1
% and stability_form = 'P' (storage <T x, P T x>, P >= eppos I) in the
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
% - st:         example-specific settings (lpisettings_sop structure)       % MMP, 11/19/2026
%               for the tests this file names                               % MMP, 11/19/2026
%
% %---------------------------------------------------------------------% %
% % Boundary-Damped Wave equation (Hyperbolic) [8] (Example 8.2)
% % PDE        u_{tt} = u_{ss}
% % with BCs   u(s=0) = 0
% %            u_{s}(s=1) = -k*u_{t}(s=1)
% %
% % We use states x1 = u_{s}, x2 = u_{t}.
% % Then x2(0) = 0,  x1(1) + k*x2(1) = 0.
% % Parameter k can be set.
% %---------------------------------------------------------------------% %

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PIETOOLS - PIETOOLS_PDE_Ex_Wave_Eq_Boundary_Damped
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
% function [PDE_t,PDE_b] = PIETOOLS_PDE_Ex_Wave_Eq_Boundary_Damped(GUI,params)  % MMP, 11/19/2026 (was)
% MMP, 11/19/2026: Release comment pass (comment-only).

% Determine the location of this example file <-- DO NOT MOVE THE FILE
loc = mfilename('fullpath');
root = fileparts(loc);

% Initialize variables
pvar s

%%% Executive Function:
evalin('base','stability = 1;');
evalin('base','stability_form = ''P'';');  % form of the stability test     % MMP, 11/19/2026

% Default parameter, overridden by the assignment strings in params
k = 1;
npars = length(params);
if npars~=0
    %%% Specify potential parameters
    for j=1:npars
        eval(params{j});
    end
end


% % % Construct the PDE.
%%% Batch input format
PDE_b.n0 = 0;   PDE_b.n1 = 2;   PDE_b.n2 = 0;
PDE_b.dom = [0,1];

PDE_b.A0 = [0 0; 0 0];
PDE_b.A1 = [0 1; 1 0];

PDE_b.B = [0 1 0 0;
    0 0 1 k];


%%% Term-based input format
PDE_t.vars = s;    PDE_t.dom = [0,1];

% PDE: x1_{t} = x2_{s},   x2_{t} = x1_{s}
PDE_t.x{1}.term{1}.x = 2;
PDE_t.x{1}.term{1}.D = 1;

PDE_t.x{2}.term{1}.x = 1;
PDE_t.x{2}.term{1}.D = 1;

% BC1: 0 = x2(0)
PDE_t.BC{1}.term{1}.x = 2;
PDE_t.BC{1}.term{1}.loc = 0;

% BC2: 0 = x1(1) + k*x2(1)
PDE_t.BC{2}.term{1}.x = 1;
PDE_t.BC{2}.term{1}.loc = 1;
PDE_t.BC{2}.term{2}.x = 2;
PDE_t.BC{2}.term{2}.loc = 1;
PDE_t.BC{2}.term{2}.C = k;


if GUI
    %%% Associated GUI save file
    app = PIETOOLS_PDE_GUI;
    load(fullfile(root,'PIETOOLS_PDE_Ex_Wave_Eq_Boundary_Damped_GUI.mat'));
    logval = app.loadData(data);
    if logval
        disp("Failed to load data object. Incorrect structure");
    end
end

% Example-specific settings                                                 % MMP, 11/19/2026
% (light tier; 1-D storage from the multiplier only, theta-length 1, both
% integral s-lengths 0; slack lift degree at its minimum Dmin)
st = lpisettings_sop('light');                                              % MMP, 11/19/2026
st.lpi.pos.len_inc{1} = 0;                                                  % MMP, 11/19/2026
st.lpi.pos.len_com{1} = 0;                                                  % MMP, 11/19/2026
st.lpi.pos.len_sos{1} = [1;2;2];                                            % MMP, 11/19/2026
st.lpi.ineq.dD = 0;                                                         % MMP, 11/19/2026

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