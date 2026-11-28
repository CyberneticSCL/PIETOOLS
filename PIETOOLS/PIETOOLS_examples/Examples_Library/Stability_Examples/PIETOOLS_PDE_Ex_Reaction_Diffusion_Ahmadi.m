function [PDE_t,PDE_b,st] = PIETOOLS_PDE_Ex_Reaction_Diffusion_Ahmadi(GUI,params) % DB, MMP, 11/22/2026
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% PIETOOLS PDE Examples
% [PDE_t,PDE_b,ST] = PIETOOLS_PDE_EX_REACTION_DIFFUSION_AHMADI(GUI,PARAMS)
% returns the 2-state reaction-diffusion equation below, in the term-based
% and in the batch format, and the recommended settings ST of its
% stability tests. It sets stability = 1, stability_dual = 1 and
% stability_form = 'P' (storage <T x, P T x>, P >= eppos I; dual form
% 'Pd') in the base workspace. With GUI = 1 it also loads the associated
% GUI save file into a new PIETOOLS_PDE_GUI window.
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
% - st:         example-specific settings (lpisettings structure)       % DB, MMP, 11/22/2026
%               for the tests this file names                               % DB, MMP, 11/22/2026
%
% %---------------------------------------------------------------------% %
% % Example Diffusion-reaction from Ahmadi (see reference below):
% % PDE        x_{t} = Cm*x + (1/R)*x_{ss} 
% % with BCs   x(s=0) = 0
% %            x(s=1) = 0
% %
% % Parameters Cm and R can be set.
% % Stable for R<2.7.
% %---------------------------------------------------------------------% %

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PIETOOLS - PIETOOLS_PDE_Ex_Reaction_Diffusion_Ahmadi
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
% DB, MMP, 11/22/2026: Changes since PIETOOLS 2025: the batch format
%   imposes x(1) = 0 in place of x_{s}(1) = 0, as the term-based format
%   and the help; new output st, the recommended settings of the stability
%   tests; stability_form = 'P' set in the base workspace; with GUI = 1 the
%   load check tests logval (it tested the undefined logvals); help and
%   comments revised.
% MMP, 11/28/2026: _sop dropped from the names of the executives, LPI
%   and settings functions (2025 names where they exist); names only.
% function [PDE_t,PDE_b] = PIETOOLS_PDE_Ex_Reaction_Diffusion_Ahmadi(GUI,params) % DB, MMP, 11/22/2026 (was)

% Determine the location of this example file <-- DO NOT MOVE THE FILE
loc = mfilename('fullpath');
root = fileparts(loc);

% Initialize variables
pvar s

%%% Executive Function:
evalin('base','stability = 1;');
evalin('base','stability_dual = 1;')
evalin('base','stability_form = ''P'';');  % form of the stability test     % DB, MMP, 11/22/2026

% Specify the parameters
Cm = [1, 1.5; 5, 0.2];      R = 2.9; % [2.93 2.94]
npars = length(params);
if npars~=0
    %%% Specify potential parameters
    for j=1:npars
        eval(params{j});
    end
end


% Construct the PDE.
ne = size(Cm,2);
on = eye(ne);      ze = zeros(ne);

%%% Batch input format
PDE_b.n0 = 0;   PDE_b.n1 = 0;   PDE_b.n2 = ne;
PDE_b.dom = [0,1];

PDE_b.A0 = Cm;
PDE_b.A2 = (1/R)*on;

PDE_b.B = [on ze ze ze;
                    ze  on ze ze];


%%% Term-based input format
PDE_t.x{1}.vars = s;   PDE_t.x{1}.dom = [0,1];

% PDE: x_{t} = Cm * x
PDE_t.x{1}.term{1}.C = Cm;

% PDE: x_{t} = ... + (1/R) * x_{ss}
PDE_t.x{1}.term{2}.D = 2;
PDE_t.x{1}.term{2}.C = (1/R)*eye(ne);

% BCs: 0 = x(0)
PDE_t.BC{1}.term{1}.loc = 0;

% BCs: 0 = x(1)
PDE_t.BC{2}.term{1}.loc = 1;


if GUI
    %%% Associated GUI save file
    app = PIETOOLS_PDE_GUI;
    load(fullfile(root,'PIETOOLS_PDE_Ex_Reaction_Diffusion_Ahmadi_GUI.mat'));
    logval = app.loadData(data);
%     if logvals                                                            % DB, MMP, 11/22/2026 (was)
    if logval                                                               % DB, MMP, 11/22/2026
        disp("Failed to load data object. Incorrect structure");
    end
end

% Example-specific settings                                                 % DB, MMP, 11/22/2026
% (stripped tier; 1-D storage without the upper-integral basis operator,
% s-length 0, and theta-lengths 1 (multiplier), 2 (lower integral))
st = lpisettings('stripped');                                               % DB, MMP, 11/22/2026
st.lpi.pos.len_com{1} = 0;                                                  % DB, MMP, 11/22/2026
st.lpi.pos.len_sos{1} = [1;2;2];                                            % DB, MMP, 11/22/2026

end
% @inproceedings{valmorbida2014semi,
%   title={Semi-definite programming and functional inequalities for distributed parameter systems},
%   author={Valmorbida, Giorgio and Ahmadi, Mohamadreza and Papachristodoulou, Antonis},
%   booktitle={53rd IEEE conference on decision and control},
%   pages={4304--4309},
%   year={2014},
%   organization={IEEE}
% }