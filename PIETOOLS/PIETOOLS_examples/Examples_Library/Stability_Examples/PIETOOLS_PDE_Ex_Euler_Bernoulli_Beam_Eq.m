function [PDE_t,PDE_b] = PIETOOLS_PDE_Ex_Euler_Bernoulli_Beam_Eq(GUI,params)
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
% % Euler-Bernoulli beam equation from Peet (Example 8.1.0.1)
% % PDE        u_{tt} = -c*u_{ssss}
% % with BCs   u(s=0) = 0
% %            u_{s}(s=0) = 0 
% %            u_{ss}(s=1) = 0 
% %            u_{sss}(s=1) = 0 
%
% % We use states x1 = u_{t}, x2 = u_{ss}, so that our system becomes
% % PDE        x1_{t} = -c * x2_{ss}
% %            x2_{t} = x1_{ss}
% % with BCs   x1(0) = 0, x2(1) = 0, x1_{s}(0) = 0, x2_{s}(1) = 0
% % 
% % Parameter c can be set.
% %---------------------------------------------------------------------% %

% Determine the location of this example file <-- DO NOT MOVE THE FILE
loc = mfilename('fullpath');
root = fileparts(loc);

% Initialize variables
pvar s theta s1 s2 theta1 theta2

%%% Executive Function:
evalin('base','stability = 1;');
evalin('base','stability_dual = 1;')

% Specify the parameters
c = 0.1;
npars = length(params);
if npars~=0
    %%% Specify potential parameters
    for j=1:npars
        eval(params{j});
    end
end


% Construct the PDE.
%%% Batch input format
PDE_b.n0 = 0;   PDE_b.n1 = 0;   PDE_b.n2 = 2;
PDE_b.dom = [0,1];

PDE_b.A2 = [0 -c; 1 0];

PDE_b.B = [ 1 0 0 0 0 0 0 0;
    0 0 0 1 0 0 0 0;
    0 0 0 0 1 0 0 0;
    0 0 0 0 0 0 0 1];


%%% Term-based input format
PDE_t.x{1}.vars = s;    PDE_t.x{1}.dom = [0,1];

% PDE: x_{t} = [0,-c;1,0] * x_{ss}
PDE_t.x{1}.term{1}.D = 2;
PDE_t.x{1}.term{1}.C = [0, -c; 1, 0];

% BCs: 0 = x1(0)
PDE_t.BC{1}.term{1}.loc = 0;
PDE_t.BC{1}.term{1}.C = [1,0];

% BCs: 0 = x2(1)
PDE_t.BC{2}.term{1}.loc = 1;
PDE_t.BC{2}.term{1}.C = [0,1];

% BCs: 0 = x1_{s}(0)
PDE_t.BC{3}.term{1}.D = 1;
PDE_t.BC{3}.term{1}.loc = 0;
PDE_t.BC{3}.term{1}.C = [1,0];

% BCs: 0 = x2_{s}(1)
PDE_t.BC{4}.term{1}.D = 1;
PDE_t.BC{4}.term{1}.loc = 1;
PDE_t.BC{4}.term{1}.C = [0,1];

%  % Note that BCs can also be coupled together:
%
%  PDE_t.BC{1}.term{1}.loc = 0;
%  PDE_t.BC{1}.term{1}.C = [1,0;0,0];
%
%  PDE_t.BC{1}.term{2}.loc = 1;
%  PDE_t.BC{1}.term{2}.C = [0,0;0,1];
%
%  PDE_t.BC{2}.term{1}.D = 1;
%  PDE_t.BC{2}.term{1}.loc = 0;
%  PDE_t.BC{2}.term{1}.C = [1,0;0,0];
%
%  PDE_t.BC{2}.term{2}.D = 1;
%  PDE_t.BC{2}.term{2}.loc = 1;
%  PDE_t.BC{2}.term{2}.C = [0,0;0,1];
% command-line formaty
% PDE_c=sys();
% x1 = state('pde');x2 = state('pde');
% eq_PDEs = [diff(x1,t) == -c*diff(x2,s,2); diff(x2,t)==diff(x1,s,2)];
% eq_BCs = [subs(x1,s,0) == 0; subs(x2,s,1) == 0; subs(diff(x1,s),s,0) ==0; subs(diff(x2,s),s,1)==0];
% PDE_c = addequation(PDE_c,[eq_PDEs;eq_BCs]);
if GUI
    %%% Associated GUI save file
    app = PIETOOLS_PDE_GUI;
    load(fullfile(root,'PIETOOLS_PDE_Ex_Euler_Bernoulli_Beam_Eq_GUI.mat'));
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