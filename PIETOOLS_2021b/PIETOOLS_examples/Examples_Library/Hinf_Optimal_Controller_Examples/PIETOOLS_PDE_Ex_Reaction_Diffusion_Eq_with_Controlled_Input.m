function [PDE_t,PDE_b] = PIETOOLS_PDE_Ex_Reaction_Diffusion_Eq_with_Controlled_Input(GUI,params)
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
% % Example pure transport equation 1D: 
% % % Stabilizing controller for Heat Equation with z=0 and w=0:
% % PDE             x_{t} = lam*x + x_{ss} + u(t)
% % With BCs        x(s=0) = 0
% %                 x(s=1) = 0
% %
% % Parameter lam can be set.
% %---------------------------------------------------------------------% %

% Determine the location of this example file <-- DO NOT MOVE THE FILE
loc = mfilename('fullpath');
root = fileparts(loc);

% Initialize variables
pvar s theta s1 s2 theta1 theta2

%%% Executive Function:
evalin('base','Hinf_control = 1;')

lam = 10;
npars = length(params);
if npars~=0
    %%% Specify potential parameters
    for j=1:npars
        eval(params{j});
    end
end


% % % Construct the PDE.
%%% Batch input format
PDE_b.nw = 0;   PDE_b.ny = 0;   PDE_b.nz = 0;   PDE_b.nx = 0;   PDE_b.nu = 1;
PDE_b.n0 = 0;   PDE_b.n1 = 0;   PDE_b.n2 = 1;
PDE_b.dom = [0,1];

PDE_b.A0 = lam;   PDE_b.A2 = 1;
PDE_b.B22 = 1;    %PDE_b.D12 = 1;

PDE_b.B = [1 0 0 0;
           0 1 0 0];


%%% Term-based input format
% Initialize 1D PDE state component.
PDE_t.x{1}.vars = s;   PDE_t.x{1}.dom = [0,1];
% Initialize a finite-dimensional input.
PDE_t.u{1}.vars = [];

% PDE: x_{t} = lam * x + x_{ss}
PDE_t.x{1}.term{1}.x = 1;
PDE_t.x{1}.term{1}.D = [0; 2];
PDE_t.x{1}.term{1}.C = [lam, 1];

% PDE: x_{t} = ... + w
PDE_t.x{1}.term{2}.u = 1;

% BC 1: 0 = x(0)
PDE_t.BC{1}.term{1}.x = 1;
PDE_t.BC{1}.term{1}.loc = 0;

% BC 2: 0 = x(1)
PDE_t.BC{2}.term{1}.x = 1;
PDE_t.BC{2}.term{1}.loc = 1;


if GUI
    %%% Associated GUI save file
    app = PIETOOLS_PDE_GUI;
    load(fullfile(root,'Examples_Libary/Hinf_Optimal_Controller_Examples/Reaction_Diffusion_Eq_with_Controlled_Input.mat'));
    logval = app.loadData(data);
    if logval
        disp("Failed to load data object. Incorrect structure");
    end
end

end