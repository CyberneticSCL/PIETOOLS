function [PDE_t,PDE_b] = PIETOOLS_PDE_Ex_Transport_Eq(GUI,params)
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
% % Transport equation: 
% % PDE:        x_{t}  = x_{s}
% % with BC     x(s=0) = 0
% %---------------------------------------------------------------------% %

% Determine the location of this example file <-- DO NOT MOVE THE FILE
loc = mfilename('fullpath');
root = fileparts(loc);

% Initialize variables
pvar s theta s1 s2 theta1 theta2

%%% Executive Function:
evalin('base','stability = 1;');
evalin('base','stability_dual = 1;');

npars = length(params);
if npars~=0
    %%% Specify potential parameters
    for j=1:npars
        eval(params{j});
    end
end


% % % Construct the PDE.
%%% Batch input format
PDE_b.n0 = 0;   PDE_b.n1 = 1;   PDE_b.n2 = 0;
PDE_b.dom = [0,1];

PDE_b.A1= 1;

on = eye(PDE_b.n1);   ze = zeros(PDE_b.n1);
PDE_b.B = [ze on];


%%% Term-based input format
% Single 1D state component
PDE_t.x{1}.vars = s;
PDE_t.x{1}.dom = [0,1];

% PDE: x_{t} = x_{s}
PDE_t.x{1}.term{1}.x = 1;
PDE_t.x{1}.term{1}.D = 1;

% BCs: 0 = x(0)
PDE_t.BC{1}.term{1}.loc = 0;


if GUI
    %%% Associated GUI save file
    app = PIETOOLS_PDE_GUI;
    load(fullfile(root,'PIETOOLS_PDE_Ex_Transport_Eq_GUI.mat'));
    logval = app.loadData(data);
    if logval
        disp("Failed to load data object. Incorrect structure");
    end
end

end