function [PDE_t,PDE_b] = PIETOOLS_PDE_Ex_Wave_Eq_Datko_Boundary_Damped_1(GUI,params)
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
% % Test 7.5c - Datko Boundary-Damped Wave equation (see reference below):
% % PDE        u_{tt} + 2*ad*u_{t} = -ad^2*u + u_{ss} 
% % with BCs   u(s=0) = 0
% %            u_{s}(s=1) = -k*u_{t}(s=1)
% %
% % We use states x1 = u_{t}, x2 = u.
% % Then x1(0) = 0, x2(0) = 0,  k*x1(1) + x2_{s}(1) = 0.
% % Parameters k and ad can be set.
% %---------------------------------------------------------------------% %

% Determine the location of this example file <-- DO NOT MOVE THE FILE
loc = mfilename('fullpath');
root = fileparts(loc);

% Initialize variables
pvar s theta s1 s2 theta1 theta2

%%% Executive Function:
evalin('base','stability = 1;');

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
PDE_t.PDE.A{4}.coeff = 1; PDE_t.PDE.A{4}.Lstate = 1; PDE_t.PDE.A{4}.Rstate = 2; PDE_t.PDE.A{4}.D = 2;
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