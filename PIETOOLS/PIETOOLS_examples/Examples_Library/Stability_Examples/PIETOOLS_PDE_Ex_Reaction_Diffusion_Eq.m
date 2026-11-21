function [PDE_t,PDE_b,st] = PIETOOLS_PDE_Ex_Reaction_Diffusion_Eq(GUI,params) % MMP, 11/21/2026
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% PIETOOLS PDE Examples
% [PDE_t,PDE_b,ST] = PIETOOLS_PDE_EX_REACTION_DIFFUSION_EQ(GUI,PARAMS)
% returns the reaction-diffusion equation below with Dirichlet BCs, in the
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
% - st:         example-specific settings (lpisettings_sop structure)       % MMP, 11/21/2026
%               for the tests this file names                               % MMP, 11/21/2026
%
% %---------------------------------------------------------------------% %
% % Scalable Diffusion Equation on [0,1] adapted from Ahmadi 2015:
% % PDE        x_{t} = lam*x + x_{ss} 
% % with BCs   x(s=0) = 0
% %            x(s=1) = 0
% %
% % Parameters lam and ne (state size) can be set.
% % Stable for lam < pi^2 = 9.8696.
% %---------------------------------------------------------------------% %

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PIETOOLS - PIETOOLS_PDE_Ex_Reaction_Diffusion_Eq
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
% MMP, 11/21/2026: Changes since PIETOOLS 2025: new output st, the
%   recommended settings of the stability test; stability_form = 'P' set
%   in the base workspace; help and comments revised.
% function [PDE_t,PDE_b] = PIETOOLS_PDE_Ex_Reaction_Diffusion_Eq(GUI,params) % MMP, 11/21/2026 (was)

% Determine the location of this example file <-- DO NOT MOVE THE FILE
loc = mfilename('fullpath');
root = fileparts(loc);

% Initialize variables
pvar s

%%% Executive Function:
evalin('base','stability = 1;');
evalin('base','stability_form = ''P'';');  % form of the stability test     % MMP, 11/21/2026

% Specify the parameters
ne = 1;   lam = 9.86; %9.86
npars = length(params);
if npars~=0
    %%% Specify potential parameters
    for j=1:npars
        eval(params{j});
    end
end


% Construct the PDE.
%%% Batch input format
PDE_b.n0 = 0;   PDE_b.n1 = 0;   PDE_b.n2 = ne;
PDE_b.dom = [0,1];

PDE_b.A0 = lam*eye(ne);
PDE_b.A2 = 1*eye(ne);

PDE_b.B=[eye(ne) zeros(ne) zeros(ne) zeros(ne);
         zeros(ne) eye(ne) zeros(ne)  zeros(ne)];


%%% Term-based input format
PDE_t.x{1}.vars = s;   PDE_t.x{1}.dom = [0,1];

% PDE: x_{t} = lam*x
PDE_t.x{1}.term{1}.C = lam * eye(ne);

% PDE: x_{t} = ... + x_{ss}
PDE_t.x{1}.term{2}.D = 2;
PDE_t.x{1}.term{2}.C = eye(ne);

% BCs: 0 = x(s=0)
PDE_t.BC{1}.term{1}.loc = 0;

% BCs: 0 = x(s=1)
PDE_t.BC{2}.term{1}.loc = 1;


if GUI
    %%% Associated GUI save file
    app = PIETOOLS_PDE_GUI;
    load(fullfile(root,'PIETOOLS_PDE_Ex_Reaction_Diffusion_Eq_1_GUI.mat'));
    logval = app.loadData(data);
    if logval
        disp("Failed to load data object. Incorrect structure");
    end
end

% Example-specific settings                                                 % MMP, 11/21/2026
% (stripped tier; 1-D storage from the lower-integral basis operator only,
% theta-length 2: multiplier theta-length 0, upper-integral s-length 0)
st = lpisettings_sop('stripped');                                           % MMP, 11/21/2026
st.lpi.pos.len_com{1} = 0;                                                  % MMP, 11/21/2026
st.lpi.pos.len_sos{1} = [0;2;2];                                            % MMP, 11/21/2026

end
% @article{valmorbida2015stability,
%   title={Stability analysis for a class of partial differential equations via semidefinite programming},
%   author={Valmorbida, Giorgio and Ahmadi, Mohamadreza and Papachristodoulou, Antonis},
%   journal={IEEE Transactions on Automatic Control},
%   volume={61},
%   number={6},
%   pages={1649--1654},
%   year={2015},
%   publisher={IEEE}
% }
