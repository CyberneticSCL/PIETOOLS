function [PDE_t,PDE_b] = PIETOOLS_PDE_Ex_Heat_Eq_with_ODE_Das(GUI,params)
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
% % Heat Equation coupled with ODE at the boundary from Das (Example 2):
% % ODE        xo_{t} = A * xo + Bxr * x_{s}(s=a)
% % PDE        x_{t} = lam * x + x_{ss} + Bpv * xo
% % with BCs   x(s=a) = 0
% %            x(s=b) = 0
% %
% %
% % Parameters lam, a and b (domain) can be set.
% % Stable for lam<pi^2
% %---------------------------------------------------------------------% %

% Determine the location of this example file <-- DO NOT MOVE THE FILE
loc = mfilename('fullpath');
root = fileparts(loc);

% Initialize variables
pvar s theta s1 s2 theta1 theta2

%%% Executive Function:
evalin('base','stability = 1;');
evalin('base','stability_dual = 1;')

lam = pi^2-1;  a = 0;  b = 1;
npars = length(params);
if npars~=0
    %%% Specify potential parameters
    for j=1:npars
        eval(params{j});
    end
end


% Construct the PDE.
%%% Batch input format
PDE_b.n0 = 0;   PDE_b.n1 = 0;   PDE_b.n2 = 2;   PDE_b.nx = 4;
PDE_b.dom = [a,b];

% ODE
PDE_b.A = [-1.2142,  1.9649,  0.2232,  0.5616;
    -1.8042, -0.7260, -0.3479,  5.4355;
    -0.2898,  0.7381, -1.7606,  0.8294;
    -0.9417, -5.3399, -1.0704, -0.7590];
PDE_b.E0 = [-1.5368 0;0 0.8871;1.0656 0;1.1882 0] * [zeros(PDE_b.n2) zeros(PDE_b.n2) eye(PDE_b.n2) zeros(PDE_b.n2)];

% PDE
PDE_b.A0 = lam*eye(PDE_b.n2);   PDE_b.A2 = eye(PDE_b.n2);
PDE_b.E = [-2.5575 0 1.0368 0;-1.8067 0.4630 1.3621 0];

% BCs
PDE_b.B = [eye(PDE_b.n2), zeros(PDE_b.n2), zeros(PDE_b.n2),zeros(PDE_b.n2);
    zeros(PDE_b.n2), eye(PDE_b.n2), zeros(PDE_b.n2), zeros(PDE_b.n2)];


%%% Term-based input format
% Initialize ODE state component.
PDE_t.x{1}.vars = [];
% Initialize PDE state component.
PDE_t.x{2}.vars = s;   PDE_t.x{2}.dom = [a,b];

% ODE: xo_{t} = A * xo;
PDE_t.x{1}.term{1}.x = 1;
PDE_t.x{1}.term{1}.C = [-1.2142,  1.9649,  0.2232,  0.5616;
    -1.8042, -0.7260, -0.3479,  5.4355;
    -0.2898,  0.7381, -1.7606,  0.8294;
    -0.9417, -5.3399, -1.0704, -0.7590];

% ODE: xo_{t} = ... + Bxr * x_{s}(s=0)
PDE_t.x{1}.term{2}.x = 2;
PDE_t.x{1}.term{2}.loc = 0;
PDE_t.x{1}.term{2}.C = [-1.5368 0;0 0.8871;1.0656 0;1.1882 0];

% PDE: x_{t} = lam * x
PDE_t.x{2}.term{1}.x = 2;
PDE_t.x{2}.term{1}.C = lam*eye(2);

% PDE: x_{t} = ... + x_{ss}
PDE_t.x{2}.term{2}.x = 2;
PDE_t.x{2}.term{2}.D = 2;

% PDE: x_{t} = ... + Bpv * v
PDE_t.x{2}.term{3}.x = 1;
PDE_t.x{2}.term{3}.C = [-2.5575 0 1.0368 0;-1.8067 0.4630 1.3621 0];

% BCs: 0 = x(0)
PDE_t.BC{1}.term{1}.x = 2;
PDE_t.BC{1}.term{1}.loc = a;
PDE_t.BC{1}.term{1}.C = eye(2);

% BCs: 0 = x(1)
PDE_t.BC{2}.term{1}.x = 2;
PDE_t.BC{2}.term{1}.loc = b;
PDE_t.BC{2}.term{1}.C = eye(2);
% command-line format
% PDE_c=sys();
% xo=state('ode');x=state('pde');
% eq_ODEPDE=[diff(xo,t) == A*xo+Bxr*subs(diff(x,s),s,a)
%                         diff(x,t) == lam*x+diff(x,s,2)+Bpv*xo];
% eq_BCs=[subs(x,s,a)==0;subs(x,s,b)==0];
% PDE_c=addequation(PDE_c,[eq_ODEPDE;eq_BCs]);
if GUI
    %%% Associated GUI save file
    app = PIETOOLS_PDE_GUI;
    load(fullfile(root,'PIETOOLS_PDE_Ex_Heat_Eq_with_ODE_Das_GUI.mat'));
    logval = app.loadData(data);
    if logval
        disp("Failed to load data object. Incorrect structure");
    end
end

end
% @article{das2018representation,
% author = {Das, Amritam and Shivakumar, Sachin and Weiland, Siep and Peet, Matthew},
% year = {2018},
% pages = {},
% title = {Representation and Stability Analysis of PDE-ODE Coupled Systems}
% }