function [PDE_t,PDE_b] = PIETOOLS_PDE_Ex_Timoshenko_Beam_2(GUI,params)
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
% % Timoschenko Beam equation Example (hyperbolic/diffusive) from Peet:
% % PDE        r*aa * w_{tt} = k*aa*g * (-phi_{s} + w_{ss})
% %            r*II * phi_{tt} = E*II * phi_{ss} + k*aa*g * (w_{s} - phi)
% % with BCs   phi(s=0) = 0
% %            w(s=0) = 0
% %            phi_{s}(s=1) = 0
% %            w_{s}(s=1) - phi(s=1) = 0
% %
% % We now use states x1 = w_{t}, x2 = w_{s}, x3 = phi_{t}, x4 = phi.
% % We also assume all coefficients are 1. Then, our system becomes:
% % PDE        x1_{t} = x2_{s} - x4_{s}
% %            x2_{t} = x1_{s}
% %            x3_{t} = x2 - x4
% %            x4_{t} = x3
% % with BCs   x4(0) = 0, x4_{s}(1) = 0, x3(0) = 0, x1(0) = 0, x2(1)-x4(1) = 0
% %
% % Unstable.
% %---------------------------------------------------------------------% %

% Determine the location of this example file <-- DO NOT MOVE THE FILE
loc = mfilename('fullpath');
root = fileparts(loc);

% Initialize variables
pvar s theta s1 s2 theta1 theta2

%%% Executive Function:
evalin('base','stability = 1;');

% Specify the parameters
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
    0 1 0 -1;
    0 0 1 0];
PDE_b.A1 = [0 1 0 -1;
    1 0 0 0;
    0 0 0 0;
    0 0 0 0];
PDE_b.A2 = [0;0;1;0];

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

% PDE 3: x3_{t} = x2 - x4
PDE_t.x{3}.term{1}.x = [2;4];
PDE_t.x{3}.term{1}.C = [1,-1];

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