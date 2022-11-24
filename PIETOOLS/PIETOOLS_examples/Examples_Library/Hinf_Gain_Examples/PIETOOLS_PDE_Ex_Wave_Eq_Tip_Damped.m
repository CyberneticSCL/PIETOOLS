function [PDE_t,PDE_b] = PIETOOLS_PDE_Example_Wave_Eq_Tip_Damped(GUI,params)
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
% % Example tip damped wave equation:
% % PDE         phi_{tt} = phi_{ss} + w(t)
% % With BCs    phi(s=0) = 0
% %             phi_{s}(s=1) = -k*phi_{t}(s=1)
% % And output  z(t) = int(phi_{t}(t,s),s,0,1)
% %
% % We use states x1 = phi_{s}, x2 = phi_{t}
% % Then x2(0) = 0,  x1(1) + k*x2(1) = 0.
% %
% % Parameter k can be set.
% % gamma = 2 for k = 0.5.
% %---------------------------------------------------------------------% %

% Determine the location of this example file <-- DO NOT MOVE THE FILE
loc = mfilename('fullpath');
root = fileparts(loc);

% Initialize variables
pvar s theta s1 s2 theta1 theta2

%%% Executive Function:
evalin('base','Hinf_gain = 1;')
evalin('base','Hinf_gain_dual = 1;')

k = 0.5;
npars = length(params);
if npars~=0
    %%% Specify potential parameters
    for j=1:npars
        eval(params{j});
    end
end


% % % Construct the PDE.
%%% Batch input format
PDE_b.n0 = 0;   PDE_b.n1 = 2;   PDE_b.n2 = 0;   PDE_b.nw = 1;   PDE_b.nz = 1;
PDE_b.dom = [0,1];
np = PDE_b.n0 + PDE_b.n1 + PDE_b.n2;

PDE_b.A1= [0 1; 1 0];
PDE_b.B21 = [0; 1];
PDE_b.Ca1 = [0 1];

PDE_b.B = [0 1 0 0; 0 0 1 k];


%%% Term-based input format
% Initialize 1D PDE state components.
PDE_t.x{1}.vars = s;   PDE_t.x{1}.dom = [0,1];
PDE_t.x{2}.vars = s;   PDE_t.x{2}.dom = [0,1];
% Initialize finite-dimensional inputs and outputs.
PDE_t.w{1}.vars = [];  PDE_t.z{1}.vars = [];

% PDE: x1_{t} = x2_{s}
PDE_t.x{1}.term{1}.x = 2;
PDE_t.x{1}.term{1}.D = 1;

% PDE: x2_{t} = x1_{s}
PDE_t.x{2}.term{1}.x = 1;
PDE_t.x{2}.term{1}.D = 1;

% PDE: x2_{t} = ... + w
PDE_t.x{2}.term{2}.w = 1;

% Output: z = int_{0}^{1} x2(s) ds
PDE_t.z{1}.term{1}.x = 2;
PDE_t.z{1}.term{1}.I{1} = PDE_t.x{2}.dom;

% BC 1: 0 = x2(0)
PDE_t.BC{1}.term{1}.x = 2;
PDE_t.BC{1}.term{1}.loc = 0;

% BC 2: 0 = x1(1)...
PDE_t.BC{2}.term{1}.x = 1;
PDE_t.BC{2}.term{1}.loc = 1;
%...+ k*x2(1)
PDE_t.BC{2}.term{2}.x=2;
PDE_t.BC{2}.term{2}.loc=1;
PDE_t.BC{2}.term{2}.C=k;

if GUI
    %%% Associated GUI save file
    app = PIETOOLS_PDE_GUI;
    load(fullfile(root,'PIETOOLS_PDE_Ex_Wave_Eq_Tip_Damped_GUI.mat'));
    logval = app.loadData(data);
    if logval
        disp("Failed to load data object. Incorrect structure");
    end
end

end