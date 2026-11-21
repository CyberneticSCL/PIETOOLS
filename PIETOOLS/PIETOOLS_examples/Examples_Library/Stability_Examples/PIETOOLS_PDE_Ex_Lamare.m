function [PDE_t,PDE_b,st] = PIETOOLS_PDE_Ex_Lamare(index,GUI,params)        % MMP, 11/21/2026
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% PIETOOLS PDE Examples
% [PDE_t,PDE_b,ST] = PIETOOLS_PDE_EX_LAMARE(INDEX,GUI,PARAMS) returns the
% 2-state linear hyperbolic system below with the parameters of Lamare
% 2016, Example 5.1-5.4, selected by INDEX, in the term-based and in the
% batch format, and the recommended settings ST of its stability test. It
% sets stability = 1 and stability_form = 'P' (storage <T x, P T x>,
% P >= eppos I) in the base workspace. With GUI = 1 it also loads the
% associated GUI save file into a new PIETOOLS_PDE_GUI window.
%
% INPUT
% - index:      2 or 2.1, 2.2, 2.3, 2.4: Example 5.1, 5.2, 5.3, 5.4 of
%               Lamare 2016; any other value raises an error.
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
%               for the test this file names, by case                       % MMP, 11/21/2026
%
% %---------------------------------------------------------------------% %
% % Lamare 2016 Ex 5 (see reference below):
% % PDE          x_{t} = Fm*x - Lm*x_{s}
% % with BC [x_-(s=1)] = [Gm1, Gm2] [x_-(s=0)]
% %         [x_+(s=0)] = [Gm3, Gm4] [x_+(s=1)]
% %---------------------------------------------------------------------% %

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PIETOOLS - PIETOOLS_PDE_Ex_Lamare
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
%   recommended settings of the stability test by case; stability_form =
%   'P' set in the base workspace; help and comments revised.
% function [PDE_t,PDE_b] = PIETOOLS_PDE_Ex_Lamare(index,GUI,params)         % MMP, 11/21/2026 (was)

% Determine the location of this example file <-- DO NOT MOVE THE FILE
loc = mfilename('fullpath');
root = fileparts(loc);

% Initialize variables
pvar s

%%% Executive Function:
evalin('base','stability = 1;');
%evalin('base','stability_dual = 1;');
evalin('base','stability_form = ''P'';');  % form of the stability test     % MMP, 11/21/2026

%%%% There are several examples included here. Add a decimal to
%%%% your example input to specify a particular set of parameters.
if index==2 || index==2.1
    % % Example 5.1
    Gm1=[.2];
    Gm2=[-.3];
    Gm3=[.6];
    Gm4=[.1];
    Lm=[-3 0;0 1];
    Fm=[.2 -.3; .6 .1];
    dir = fullfile(root,'PIETOOLS_PDE_Ex_Lamare_Ex_5_1_GUI.mat');
elseif index==2.2
    % % Example 5.2
    Gm1=[.1];
    Gm2=[-.8];
    Gm3=[.6];
    Gm4=[-.4];
    Lm=[-1 0;0 1];
    Fm=[-.3 .1; .1 -.3];
    dir = fullfile(root,'PIETOOLS_PDE_Ex_Lamare_Ex_5_2_GUI.mat');
elseif index==2.3
    % % Example 5.3
    Gm1=[-.2202];
    Gm2=[1.3955];
    Gm3=[-.0596];
    Gm4=[.2090];
    Lm=[-1 0;0 2];
    Fm=[-.1 .1; .5 -.8];
    dir = fullfile(root,'PIETOOLS_PDE_Ex_Lamare_Ex_5_3_GUI.mat');
elseif index==2.4
    % % Example 5.4
    Gm1=[.5];
    Gm2=[-.4];
    Gm3=[.2];
    Gm4=[.8];
    Lm=[-2 0;0 1];
    Fm=[-.6565 -.3743; -.113 -.6485];
    dir = fullfile(root,'PIETOOLS_PDE_Ex_Lamare_Ex_5_4_GUI.mat');
else
    error(['Example ',num2str(index),'does not exist, please choose from 2.1, 2.2, 2.3 or 2.4']);
end

npars = length(params);
if npars~=0
    %%% Specify potential parameters
    for j=1:npars
        eval(params{j});
    end
end

% % % Construct the PDE.
%%% Batch input format
PDE_b.n0 = 0;   PDE_b.n1 = 2;   PDE_b.n2=0;
PDE_b.dom = [0,1];
PDE_b.A0 = Fm;
PDE_b.A1 = -Lm;
ny1 = size(Gm1,1);   ny2 = size(Gm3,1);
on0 = eye(ny1);   on1 = eye(ny2);   zer12 = zeros(ny1,ny2);
PDE_b.B = [-Gm1 zer12 on0 -Gm2;
           -Gm3 on1 zer12' -Gm4];


%%% Term-based input format
% A single state component.
PDE_t.x{1}.vars = s;
PDE_t.x{1}.dom = [0,1];

% PDE: x_{t} = Fm*x
PDE_t.x{1}.term{1}.x = 1;
PDE_t.x{1}.term{1}.C = Fm;

% PDE: x_{t} = ... -Lm*x_{s}
PDE_t.x{1}.term{2}.x = 1;
PDE_t.x{1}.term{2}.D = 1;
PDE_t.x{1}.term{2}.C = -Lm;

ny1 = size(Gm1,1);    ny2 = size(Gm3,1);
on0 = eye(ny1);       on1 = eye(ny2);     zer12 = zeros(ny1,ny2);

% BCs: 0 = [Gm1,0;-Gm3,1] * x(0)
PDE_t.BC{1}.term{1}.x = 1;
PDE_t.BC{1}.term{1}.loc = 0;
PDE_t.BC{1}.term{1}.C = [-Gm1, zer12; -Gm3, on1];

% BCs: 0 = ... + [1,-Gm2;0,-Gm4] * x(1)
PDE_t.BC{1}.term{2}.x = 1;
PDE_t.BC{1}.term{2}.loc = 1;
PDE_t.BC{1}.term{2}.C = [on0, -Gm2; zer12', -Gm4];

% Example-specific settings, by case                                        % MMP, 11/21/2026
% (stripped tier with the 1-D storage from the multiplier only, both
% integral s-lengths 0, and the slack lift degree at its minimum Dmin; for
% 2.2 and 2.4 the multiplier theta-length 1)
st = lpisettings_sop('stripped');                                           % MMP, 11/21/2026
st.lpi.pos.len_inc{1} = 0;                                                  % MMP, 11/21/2026
st.lpi.pos.len_com{1} = 0;                                                  % MMP, 11/21/2026
st.lpi.ineq.dD = 0;                                                         % MMP, 11/21/2026
if index==2.2 || index==2.4                                                 % MMP, 11/21/2026
    st.lpi.pos.len_sos{1} = [1;2;2];                                        % MMP, 11/21/2026
end                                                                         % MMP, 11/21/2026

if GUI
    %%% Associated GUI save file
    app = PIETOOLS_PDE_GUI;
    load(dir);
    logval = app.loadData(data);
    if logval
        disp("Failed to load data object. Incorrect structure");
    end
end

end
% @article{lamare2016optimisation,
%   title={An optimisation approach for stability analysis and controller synthesis of linear hyperbolic systems},
%   author={Lamare, Pierre-Olivier and Girard, Antoine and Prieur, Christophe},
%   journal={ESAIM: Control, Optimisation and Calculus of Variations},
%   volume={22},
%   number={4},
%   pages={1236--1263},
%   year={2016},
%   publisher={EDP Sciences}