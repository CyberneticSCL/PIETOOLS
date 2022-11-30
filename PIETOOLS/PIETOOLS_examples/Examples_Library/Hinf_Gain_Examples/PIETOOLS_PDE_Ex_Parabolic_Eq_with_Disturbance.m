function [PDE_t,PDE_b] = PIETOOLS_PDE_Ex_Parabolic_Eq_with_Disturbance(GUI,params)
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
% % gamma = 15.147 for lam = 4.6.
% %---------------------------------------------------------------------% %

% Determine the location of this example file <-- DO NOT MOVE THE FILE
loc = mfilename('fullpath');
root = fileparts(loc);

% Initialize variables
pvar s theta s1 s2 theta1 theta2

%%% Executive Function:
evalin('base','Hinf_gain = 1;')
%%% evalin('base','Hinf_gain_dual = 1;')

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
PDE_b.B21 = s;
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
disp('Warning: With sedumi solver, this example requires at least settings veryheavy.')
end
% @inproceedings{shivakumar2019computing,
%   title={Computing input-ouput properties of coupled linear pde systems},
%   author={Shivakumar, Sachin and Peet, Matthew M},
%   booktitle={2019 American Control Conference (ACC)},
%   pages={606--613},
%   year={2019},
%   organization={IEEE}
% }