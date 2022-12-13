function [PDE_t,PDE_b] = PIETOOLS_PDE_Ex_Timoshenko_Beam_1(GUI,params)
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
% % Timoschenko beam equation (hyperbolic) from Peet (Example 8.1.0.2):
% % PDE        r*aa * w_{tt} = k*aa*g * (-phi_{s} + w_{ss})
% %            r*II * phi_{tt} = E*II * phi_{ss} + k*aa*g * (w_{s} - phi)
% % with BCs   phi(s=0) = 0
% %            w(s=0) = 0 
% %            phi_{s}(s=1) = 0 
% %            w_{s}(s=1) - phi(s=1) = 0 
%
% % We use states x1 = w_{t}, x2 = k*aa*g * (w_{s}-phi), x3 = phi_{t}, x4 = E*II * phi_{s}.
% % Then, our system becomes:
% % PDE        x1_{t} = (1/r/aa) * x2_{s}
% %            x2_{t} = k*aa*g * x1_{s} - k*aa*g * x3
% %            x3_{t} = (1/r/II) * x2 + (1/r/II) * x4_{s}
% %            x4_{t} = E*II * x3_{s}
% % with BCs   x1(0) = 0, x3(0) = 0, x4(1) = 0, x2(1) = 0
% % 
% % Parameters k, aa, II, g, E and r can be set.
% %---------------------------------------------------------------------% %

% Determine the location of this example file <-- DO NOT MOVE THE FILE
loc = mfilename('fullpath');
root = fileparts(loc);

% Initialize variables
pvar s theta s1 s2 theta1 theta2

%%% Executive Function:
evalin('base','stability = 1;');

% Specify the parameters
k = 1;    aa = 1;    II = 1;    g = 1;    E = 1;   r = 1;  %5-30kPa
npars = length(params);
if npars~=0
    %%% Specify potential parameters
    for j=1:npars
        eval(params{j});
    end
end


% Construct the PDE.
%%% Batch input format
PDE_b.n0 = 0;   PDE_b.n1 = 4;   PDE_b.n2 = 0;
PDE_b.dom = [0,1];

PDE_b.A0 = [0 0 0 0;
    0 0 -k*aa*g 0;
    0 1/r/II 0 0;
    0 0 0 0];
PDE_b.A1 = [0 1/r/aa 0 0;
    k*aa*g 0 0 0;
    0 0 0 1/r/II;
    0 0 E*II 0];

PDE_b.B = [ 1 0 0 0 0 0 0 0;
    0 0 1 0 0 0 0 0;
    0 0 0 0 0 0 0 1;
    0 0 0 0 0 1 0 0];


%%% Term-based input format
PDE_t.x{1}.vars = s;    PDE_t.x{1}.dom = [0,1];

% PDE: x2_{t} = -k*aa*g * x3,   x3_{t} = (1/r/II) * x2
PDE_t.x{1}.term{1}.C = [0, 0,      0,       0;
    0, 0,      -k*aa*g, 0;
    0, 1/r/II, 0,       0;
    0, 0,      0,       0];

% PDE: x1_{t} = (1/r/aa) * x2_{s},        x2_{t} = ... + k*aa*g * x1_{s}
%      x3_{t} = ... + (1/r/II) * x4_{s}   x4_{t} = E*II * x3_{s}
PDE_t.x{1}.term{2}.D = 1;
PDE_t.x{1}.term{2}.C = [0,      1/r/aa, 0,    0;
    k*aa*g, 0,      0,    0;
    0,      0,      0,    1/r/II;
    0,      0,      E*II, 0];

% BCs: 0 = x1(0),   0 = x3(0)
PDE_t.BC{1}.term{1}.loc = 0;
PDE_t.BC{1}.term{1}.C = [1,0,0,0;0,0,1,0;0,0,0,0;0,0,0,0];

% BCs: 0 = x4(1),   0 = x2(1)
PDE_t.BC{1}.term{2}.loc = 1;
PDE_t.BC{1}.term{2}.C = [0,0,0,0;0,0,0,0;0,0,0,1;0,1,0,0];


if GUI
    %%% Associated GUI save file
    app = PIETOOLS_PDE_GUI;
    load(fullfile(root,'PIETOOLS_PDE_Ex_Timoshenko_Beam_1_GUI.mat'));
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