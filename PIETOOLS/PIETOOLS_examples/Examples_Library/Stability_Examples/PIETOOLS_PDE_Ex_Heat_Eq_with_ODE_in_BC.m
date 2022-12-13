function [PDE_t,PDE_b] = PIETOOLS_PDE_Ex_Heat_Eq_with_ODE_in_BC(GUI,params)
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
% % Heat Equation coupled with ODE at the boundary:
% % ODE        xo_{t} = k * xo
% % PDE        x_{t} = x_{ss} 
% % with BCs   x(s=0) = 0
% %            x(s=1) = xo
% %
% % Parameter k can be set.
% % Stable for k<0
% %---------------------------------------------------------------------% %

% Determine the location of this example file <-- DO NOT MOVE THE FILE
loc = mfilename('fullpath');
root = fileparts(loc);

% Initialize variables
pvar s theta s1 s2 theta1 theta2

%%% Executive Function:
evalin('base','stability = 1;');
evalin('base','stability_dual = 1;')

k = -1;
npars = length(params);
if npars~=0
    %%% Specify potential parameters
    for j=1:npars
        eval(params{j});
    end
end


% Construct the PDE.
%%% Batch input format
PDE_b.n0 = 0;   PDE_b.n1 = 0;   PDE_b.n2 = 1;   PDE_b.nx = 1;
PDE_b.dom = [0,1];

PDE_b.A2 = 1;  PDE_b.A = k;

PDE_b.B = [1 0 0 0; 0 1 0 0]; PDE_b.Bx = [0;1];


%%% Term-based input format
% Initialize ODE state component.
PDE_t.x{1}.vars = [];
% Initialize 1D PDE state component.
PDE_t.x{2}.vars = s;   PDE_t.x{2}.dom = [0,1];

% ODE: xo_{t} = k * xo;
PDE_t.x{1}.term{1}.x = 1;
PDE_t.x{1}.term{1}.C = k;

% PDE: x_{t} = x_{ss}
PDE_t.x{2}.term{1}.x = 2;
PDE_t.x{2}.term{1}.D = 2;

% BCs: 0 = x(0)
PDE_t.BC{1}.term{1}.x = 2;
PDE_t.BC{1}.term{1}.loc = 0;

% BCs: 0 = x(1) - xo
PDE_t.BC{2}.term{1}.x = 2;
PDE_t.BC{2}.term{1}.loc = 1;

PDE_t.BC{2}.term{2}.x = 1;
PDE_t.BC{2}.term{2}.C = -1;


if GUI
    %%% Associated GUI save file
    app = PIETOOLS_PDE_GUI;
    load(fullfile(root,'PIETOOLS_PDE_Ex_Heat_Eq_with_ODE_in_BC_GUI.mat'));
    logval = app.loadData(data);
    if logval
        disp("Failed to load data object. Incorrect structure");
    end
end

end