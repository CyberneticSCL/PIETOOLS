function [PDE_t,PDE_b] = PIETOOLS_PDE_Example_Lamare(index,GUI,params)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% PIETOOLS PDE Examples
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
%
% %---------------------------------------------------------------------% %
% % Lamare 2016 Ex 5 (see reference below):
% % PDE          x_{t} = Fm*x - Lm*x_{s}
% % with BC [x_-(s=1)] = [Gm1, Gm2] [x_-(s=0)]
% %         [x_+(s=0)] = [Gm3, Gm4] [x_+(s=1)]
% %---------------------------------------------------------------------% %

% Determine the location of this example file <-- DO NOT MOVE THE FILE
loc = mfilename('fullpath');
root = fileparts(loc);

% Initialize variables
pvar s theta s1 s2 theta1 theta2

%%% Executive Function:
evalin('base','stability = 1;');
%evalin('base','stability_dual = 1;');

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
    dir = fullfile(root,'PIETOOLS_PDE_Example_Lamare_Ex_5_1_GUI.mat');
elseif index==2.2
    % % Example 5.2
    Gm1=[.1];
    Gm2=[-.8];
    Gm3=[.6];
    Gm4=[-.4];
    Lm=[-1 0;0 1];
    Fm=[-.3 .1; .1 -.3];
    dir = fullfile(root,'PIETOOLS_PDE_Example_Lamare_Ex_5_2_GUI.mat');
elseif index==2.3
    % % Example 5.3
    Gm1=[-.2202];
    Gm2=[1.3955];
    Gm3=[-.0596];
    Gm4=[.2090];
    Lm=[-1 0;0 2];
    Fm=[-.1 .1; .5 -.8];
    dir = fullfile(root,'PIETOOLS_PDE_Example_Lamare_Ex_5_3_GUI.mat');
elseif index==2.4
    % % Example 5.4
    Gm1=[.5];
    Gm2=[-.4];
    Gm3=[.2];
    Gm4=[.8];
    Lm=[-2 0;0 1];
    Fm=[-.6565 -.3743; -.113 -.6485];
    dir = fullfile(root,'PIETOOLS_PDE_Example_Lamare_Ex_5_4_GUI.mat');
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