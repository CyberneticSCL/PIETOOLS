function [PDE_t,PDE_b] = PIETOOLS_PDE_Ex_Heat_Eq_with_Distributed_Disturbance(GUI,params)
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
% % Heat equation with distributed disturbance:
% % PDE         x_{t} = x_{ss} + s*w(t)
% % With BCs    x(s=0) = 0
% %             x_{s}(s=1) = 0
% % And Output  z(t) = int(x(t,s),s,0,1)
% %
% % gamma = 0.3333.
% %---------------------------------------------------------------------% %

% Determine the location of this example file <-- DO NOT MOVE THE FILE
loc = mfilename('fullpath');
root = fileparts(loc);

% Initialize variables
pvar s theta s1 s2 theta1 theta2

%%% Executive Function:
evalin('base','Hinf_gain = 1;')

ne = 1;
npars = length(params);
if npars~=0
    %%% Specify potential parameters
    for j=1:npars
        eval(params{j});
    end
end


% % % Construct the PDE.
on = eye(ne);     ze = zeros(ne);

%%% Batch input format
PDE_b.n0 = 0;   PDE_b.n1 = 0;   PDE_b.n2 = ne;   PDE_b.nw = ne;   PDE_b.nz = ne;
PDE_b.dom = [0,1];

PDE_b.A2 = on;    PDE_b.B21 = s*on;
PDE_b.Ca1 = on;

PDE_b.B = [on ze ze ze;
           ze ze ze on];


%%% Term-based input format
% Initialize 1D PDE state component.
PDE_t.x{1}.vars = s;   PDE_t.x{1}.dom = [0,1];
% Initialize finite-dimensional inputs and outputs.
PDE_t.w{1}.vars = [];  PDE_t.z{1}.vars = [];

% PDE: x_{t} = x_{ss}
PDE_t.x{1}.term{1}.x = 1;
PDE_t.x{1}.term{1}.D = 2;

% PDE: x_{t} = ... + s * w,
PDE_t.x{1}.term{2}.w = 1;
PDE_t.x{1}.term{2}.C = s*eye(ne);

% Output: z = int_{0}^{1} x(s) ds
PDE_t.z{1}.term{1}.x = 1;
PDE_t.z{1}.term{1}.I{1} = PDE_t.x{1}.dom;

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
    load(fullfile(root,'PIETOOLS_PDE_Ex_Heat_Eq_with_Distributed_Disturbance_GUI.mat'));
    logval = app.loadData(data);
    if logval
        disp("Failed to load data object. Incorrect structure");
    end
end

end