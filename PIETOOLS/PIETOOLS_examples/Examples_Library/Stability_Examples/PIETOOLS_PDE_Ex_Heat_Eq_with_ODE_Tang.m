function [PDE_t,PDE_b] = PIETOOLS_PDE_Ex_Heat_Eq_with_ODE_Tang(GUI,params)
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
% % Coupled ODE-PDE system from Tang (see reference below):
% % ODE             xo_{t} = xo + x_{s}(s=0)
% % PDE             x_{t} = x_{ss}
% % With BCs        x(s=0) = -xo
% %                 x(s=1) = k * xo
% %
% % Parameter k can be set.
% %---------------------------------------------------------------------% %

% Determine the location of this example file <-- DO NOT MOVE THE FILE
loc = mfilename('fullpath');
root = fileparts(loc);

% Initialize variables
pvar s theta s1 s2 theta1 theta2

%%% Executive Function:
evalin('base','stability = 1;');

k = -2;
npars = length(params);
if npars~=0
    %%% Specify potential parameters
    for j=1:npars
        eval(params{j});
    end
end


% Construct the PDE.
%%% Batch input format
PDE_b.nx = 1;   PDE_b.nw = 0;   PDE_b.ny = 0;   PDE_b.nz = 0;   PDE_b.nu = 0;
PDE_b.n0 = 0;   PDE_b.n1 = 0;   PDE_b.n2 = 1;
PDE_b.dom = [0,1];

PDE_b.A = 1;   PDE_b.E0 = [0 0 1 0];
PDE_b.A2 = 1;
PDE_b.Bx = [-1; k];

PDE_b.B = [1 0 0 0;
    0 1 0 0];


%%% Term-based input format
% Initialize ODE state component.
PDE_t.x{1}.vars = [];
% Initialize 1D PDE state component.
PDE_t.x{2}.vars = s;   PDE_t.x{2}.dom = [0,1];

% ODE: xo_{t} = xo
PDE_t.x{1}.term{1}.x = 1;

% ODE: xo_{t} = ... + x_{s}(s=0)
PDE_t.x{1}.term{2}.x = 2;
PDE_t.x{1}.term{2}.D = 1;
PDE_t.x{1}.term{2}.loc = 0;

% PDE: x_{t} = x_{ss}
PDE_t.x{2}.term{1}.x = 2;
PDE_t.x{2}.term{1}.D = 2;

% BC 1: 0 = x(0)
PDE_t.BC{1}.term{1}.x = 2;
PDE_t.BC{1}.term{1}.loc = 0;

% BC 1: 0 = ... + xo
PDE_t.BC{1}.term{2}.x = 1;

% BC 2: 0 = x(1)
PDE_t.BC{2}.term{1}.x = 2;
PDE_t.BC{2}.term{1}.loc = 1;

% BC 2: 0 = ... -k * xo
PDE_t.BC{2}.term{2}.x = 1;
PDE_t.BC{2}.term{2}.C = -k;


if GUI
    %%% Associated GUI save file
    app = PIETOOLS_PDE_GUI;
    load(fullfile(root,'PIETOOLS_PDE_Ex_Heat_Eq_with_ODE_Tang_GUI.mat'));
    logval = app.loadData(data);
    if logval
        disp("Failed to load data object. Incorrect structure");
    end
end

end
% @article{tang2011state,
% title = {State and output feedback boundary control for a coupled PDE–ODE system},
% journal = {Systems & Control Letters},
% volume = {60},
% number = {8},
% pages = {540-545},
% year = {2011},
% issn = {0167-6911},
% doi = {https://doi.org/10.1016/j.sysconle.2011.04.011},
% url = {https://www.sciencedirect.com/science/article/pii/S0167691111000922},
% author = {Shuxia Tang and Chengkang Xie},
% keywords = {Coupled system, PDEs, Boundary control, Output feedback},
% abstract = {This note is devoted to stabilizing a coupled PDE–ODE system with interaction at the interface. First, a state feedback boundary controller is designed, and the system is transformed into an exponentially stable PDE–ODE cascade with an invertible integral transformation, where PDE backstepping is employed. Moreover, the solution to the resulting closed-loop system is derived explicitly. Second, an observer is proposed, which is proved to exhibit good performance in estimating the original coupled system, and then an output feedback boundary controller is obtained. For both the state and output feedback boundary controllers, exponential stability analyses in the sense of the corresponding norms for the resulting closed-loop systems are provided. The boundary controller and observer for a scalar coupled PDE–ODE system as well as the solutions to the closed-loop systems are given explicitly.}
% }