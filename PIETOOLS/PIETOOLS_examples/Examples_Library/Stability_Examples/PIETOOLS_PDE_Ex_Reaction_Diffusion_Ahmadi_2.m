function [PDE_t,PDE_b] = PIETOOLS_PDE_Ex_Reaction_Diffusion_Ahmadi_2(GUI,params)
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
% % Adapted Example B (Diffusion-Reaction equation) from Ahmadi:
% % PDE        x_{t} = Cm*x + (1/R)*x_{ss} 
% % with BCs   x(s=0) = 0
% %            x_{s}(s=1) = 0
% %
% % Parameters Cm and R can be set.
% % Stable for R<21.
% %---------------------------------------------------------------------% %

% Determine the location of this example file <-- DO NOT MOVE THE FILE
loc = mfilename('fullpath');
root = fileparts(loc);

% Initialize variables
pvar s theta s1 s2 theta1 theta2

%%% Executive Function:
evalin('base','stability = 1;');
evalin('base','stability_dual = 0;')

% Specify the parameters
Cm = [0 0 0; s 0 0; s -s^2 0];      R = (21+9);
npars = length(params);
if npars~=0
    %%% Specify potential parameters
    for j=1:npars
        eval(params{j});
    end
end


% Construct the PDE.
ne = size(Cm,2);
on = eye(ne);      ze = zeros(ne);


%%% Batch input format
PDE_b.n0 = 0;   PDE_b.n1 = 0;   PDE_b.n2 = ne;
PDE_b.dom = [0,1];

PDE_b.A0 = Cm;
PDE_b.A2 = (1/R)*on;

PDE_b.B = [on ze ze ze;
    ze on ze ze];


%%% Term-based input format
PDE_t.x{1}.vars = s;   PDE_t.x{1}.dom = [0,1];

% PDE: x_{t} = Cm * x
PDE_t.x{1}.term{1}.C = Cm;

% PDE: x_{t} = ... + (1/R) * x_{ss}
PDE_t.x{1}.term{2}.D = 2;
PDE_t.x{1}.term{2}.C = (1/R)*eye(ne);

% BCs: 0 = x(0)
PDE_t.BC{1}.term{1}.loc = 0;

% BCs: 0 = x_{s}(1)
PDE_t.BC{2}.term{1}.D = 1;
PDE_t.BC{2}.term{1}.loc = 1;


if GUI
    %%% Associated GUI save file
    app = PIETOOLS_PDE_GUI;
    load(fullfile(root,'PIETOOLS_PDE_Ex_Reaction_Diffusion_Ahmadi_2_GUI.mat'));
    logval = app.loadData(data);
    if logval
        disp("Failed to load data object. Incorrect structure");
    end
end

end
% @article{valmorbida2015stability,
%   title={Stability analysis for a class of partial differential equations via semidefinite programming},
%   author={Valmorbida, Giorgio and Ahmadi, Mohamadreza and Papachristodoulou, Antonis},
%   journal={IEEE Transactions on Automatic Control},
%   volume={61},
%   number={6},
%   pages={1649--1654},
%   year={2015},
%   publisher={IEEE}
% }