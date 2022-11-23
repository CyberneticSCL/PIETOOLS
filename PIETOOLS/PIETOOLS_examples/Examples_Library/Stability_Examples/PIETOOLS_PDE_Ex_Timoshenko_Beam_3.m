function PDE_t = PIETOOLS_PDE_Ex_Timoshenko_Beam_3(GUI,params)
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
% % Assuming sufficiently smooth solutins, we may rewrite the system as a 
% % single equation
% % PDE        0 = alp*w_{ssss} + bet*w_{tt} - gam*w_{ttss} + w_{tttt}
% %
% % We use states x1 = w_{ttt}, x2 = w_{t}, x3 = w_{tt}, x4 = w,
% % and consider BCs
% % w(s=0) = w_{s}(s=0) = 0
% % w_{ss}(s=1) - w(s=1) = w_{sss}(s=1) - w_{s}(s=1) = 0
% % w_{t}(s=0) = 0      w_{ts}(s=0) = 0
% % w_{tt}(s=0) = 0     w_{tts}(s=0) = 0
% % 
% % Parameters alp, bet and gam can be set.
% % Timoschenko Beam equation 4th order implementation - stable
% %---------------------------------------------------------------------% %

% Determine the location of this example file <-- DO NOT MOVE THE FILE
loc = mfilename('fullpath');
root = fileparts(loc);

% Initialize variables
pvar s theta s1 s2 theta1 theta2

%%% Executive Function:
evalin('base','stability = 1;');

% Specify the parameters
alp = 1;     bet = 1;     gam = 1;
npars = length(params);
if npars~=0
    %%% Specify potential parameters
    for j=1:npars
        eval(params{j});
    end
end


% Construct the PDE.
%%% Term-based input format
PDE_t.vars = s;    PDE_t.dom = [0,1];

% PDE: x1 = -bet * x3 + gam * x3_{ss}
PDE_t.x{1}.term{1}.x = 3;
PDE_t.x{1}.term{1}.D = [0;2];
PDE_t.x{1}.term{1}.C = [-bet, gam];

% PDE: x1 = ... - alp * x4_{ssss}
PDE_t.x{1}.term{2}.x = 4;
PDE_t.x{1}.term{2}.D = 4;
PDE_t.x{1}.term{2}.C = -alp;

% PDE: x2 = x3
PDE_t.x{2}.term{1}.x = 3;

% PDE: x3 = x1
PDE_t.x{3}.term{1}.x = 1;

% PDE: x4 = x2
PDE_t.x{4}.term{1}.x = 2;

% BC 1: 0 = x4(0)
PDE_t.BC{1}.term{1}.x = 4;
PDE_t.BC{1}.term{1}.loc = 0;

% BC 2: 0 = x4(1)
PDE_t.BC{2}.term{1}.x = 4;
PDE_t.BC{2}.term{1}.loc = 1;

% BC 3: 0 = x4_{ss}(1) - x4(1)
PDE_t.BC{3}.term{1}.x = 4;
PDE_t.BC{3}.term{1}.D = [2;0];
PDE_t.BC{3}.term{1}.loc = 1;
PDE_t.BC{3}.term{1}.C = [1,-1];

% BC 4: 0 = x4_{sss}(1) - x4_{s}(1)
PDE_t.BC{4}.term{1}.x = 4;
PDE_t.BC{4}.term{1}.D = [3;1];
PDE_t.BC{4}.term{1}.loc = 1;
PDE_t.BC{4}.term{1}.C = [1,-1];

% BC 5: 0 = x2(0),             % BC 7: 0 = x3(0)
PDE_t.BC{5}.term{1}.x = 2;     PDE_t.BC{7}.term{1}.x = 3;
PDE_t.BC{5}.term{1}.loc = 0;   PDE_t.BC{7}.term{1}.loc = 0;

% BC 6: 0 = x2_{s}(0),         % BC 8: 0 = x3_{s}(0)
PDE_t.BC{6}.term{1}.x = 2;     PDE_t.BC{8}.term{1}.x = 3;
PDE_t.BC{6}.term{1}.D = 1;     PDE_t.BC{8}.term{1}.D = 1;
PDE_t.BC{6}.term{1}.loc = 0;   PDE_t.BC{8}.term{1}.loc = 0;


if GUI
    %%% Associated GUI save file
    app = PIETOOLS_PDE_GUI;
    load(fullfile(root,'PIETOOLS_PDE_Ex_Timoshenko_Beam_3_GUI.mat'));
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