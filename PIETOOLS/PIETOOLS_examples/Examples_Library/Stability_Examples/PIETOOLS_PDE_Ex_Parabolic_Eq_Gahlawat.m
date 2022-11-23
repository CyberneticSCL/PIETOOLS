function [PDE_t,PDE_b] = PIETOOLS_PDE_Ex_Parabolic_Eq_Gahlawat(GUI,params)
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
% % Diffusion Equation Problem 2 from Gahlawat 2017:
% % PDE        x_{t} = c(s)*x + b(s)*x_{s} + a(s)*x_{ss} 
% % with BCs   x(s=0) = 0
% %            x_{s}(s=1) = 0
% %
% % and with   a(s) = s^3 - s^2 + 2
% %            b(s) = 3*s^2 - 2*s
% %            c(s) = -0.5*s^3 + 1.3*s^2 - 1.5*s + 0.7 + lam
% %
% % Parameters lam, ne, afun, bfun and cfun (state size) can be set.
% % Unstable for lam > 4.66.
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
ne = 1; lam = 4.66;
a = s^3 - s^2 + 2;
b = 3*s^2 - 2*s;
c = -0.5*s^3 + 1.3*s^2 - 1.5*s + 0.7 + lam;
 
npars = length(params);
if npars~=0
    %%% Specify potential parameters
    for j=1:npars
        eval(params{j});
    end
end


% Construct the PDE.
%%% Batch input format
PDE_b.n0 = 0;   PDE_b.n1 = 0;   PDE_b.n2 = ne;
PDE_b.dom = [0,1];

PDE_b.A2 = a*eye(ne);
PDE_b.A1 = b*eye(ne);
PDE_b.A0 = c*eye(ne);

PDE_b.B = [eye(ne), zeros(ne), zeros(ne), zeros(ne);
    zeros(ne), zeros(ne), zeros(ne), eye(ne)];


%%% Term-based input format
PDE_t.x{1}.vars = s;   PDE_t.x{1}.dom = [0,1];

% PDE: x_{t} = cfun * x
PDE_t.x{1}.term{1}.C = c*eye(ne);

% PDE: x_{t} = ... + bfun * x_{s}
PDE_t.x{1}.term{2}.D = 1;
PDE_t.x{1}.term{2}.C = b*eye(ne);

% PDE: x_{t} = ... + afun * x_{ss}
PDE_t.x{1}.term{3}.D = 2;
PDE_t.x{1}.term{3}.C = a;

% BCs: 0 = x(0)
PDE_t.BC{1}.term{1}.loc = 0;

% BCs: 0 = x_{s}(1)
PDE_t.BC{2}.term{1}.D = 1;
PDE_t.BC{2}.term{1}.loc = 1;


if GUI
    %%% Associated GUI save file
    app = PIETOOLS_PDE_GUI;
    load(fullfile(root,'PIETOOLS_PDE_Ex_Parabolic_Eq_Gahlawat_GUI.mat'));
    logval = app.loadData(data);
    if logval
        disp("Failed to load data object. Incorrect structure");
    end
end

end
% @article{gahlawat2016convex,
%   title={A convex sum-of-squares approach to analysis, state feedback and output feedback control of parabolic PDEs},
%   author={Gahlawat, Aditya and Peet, Matthew M},
%   journal={IEEE Transactions on Automatic Control},
%   volume={62},
%   number={4},
%   pages={1636--1651},
%   year={2016},
%   publisher={IEEE}
% }