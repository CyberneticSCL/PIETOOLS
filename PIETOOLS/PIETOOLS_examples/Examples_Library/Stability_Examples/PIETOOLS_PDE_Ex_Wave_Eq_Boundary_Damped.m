function [PDE_t,PDE_b] = PIETOOLS_PDE_Ex_Wave_Eq_Boundary_Damped(GUI,params)
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
% % Boundary-Damped Wave equation (Hyperbolic) [8] (Example 8.2)
% % PDE        u_{tt} = u_{ss}
% % with BCs   u(s=0) = 0
% %            u_{s}(s=1) = -k*u_{t}(s=1)
% %
% % We use states x1 = u_{s}, x2 = u_{t}.
% % Then x2(0) = 0,  x1(1) + k*x2(1) = 0.
% % Parameter k can be set.
% %---------------------------------------------------------------------% %

% Determine the location of this example file <-- DO NOT MOVE THE FILE
loc = mfilename('fullpath');
root = fileparts(loc);

% Initialize variables
pvar s theta s1 s2 theta1 theta2

%%% Executive Function:
evalin('base','stability = 1;');

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