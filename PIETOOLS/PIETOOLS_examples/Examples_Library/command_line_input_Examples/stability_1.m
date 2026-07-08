% 	PDE: x_{t} = v*x_{s}                            | v=-1;   (Exponentially PDE stable ) 
%   BCs: x(s=0) = 0                                 |                (Finite Energy PDE stable ) 
%clc; clear; 
clear stateNameGenerator
pvar s t; % define independent variables
%% Define dependent variables and system variable
% 	PDE: x_{t} = x_{s}                                 
%   BCs: x(s=0) = 0 
x = pde_var(s,[0,1]);
%% Specify the parameters
v=-1;
if exist('params','var')
    npars = length(params);
    %%% Specify potential parameters
    for j=1:npars
        eval(params{j});
    end
end
%% Define equations
eq_PDE = diff(x,t)==v*diff(x,s); % 	PDE: x_{t} = -x_{s}
eq_BC = subs(x,s,0)==0;        %   BCs: x(s=0) = 0 
%% initialize pde system;
PDE = initialize([eq_PDE;eq_BC]);
%% display pde to verify and convert to pie
display(PDE);
PIE = convert(PDE,'pie');
display(PIE);
%%
if exist('GUI','var')
    if GUI
            loc = mfilename('fullpath');
        if isempty(loc)
            loc = which('stability_1');
        end
        script_dir = fileparts(loc);
        examples_dir = fileparts(script_dir);
        gui_dir = fullfile(examples_dir,'Stability_Examples');
        %%% Associated GUI save file
        app = PIETOOLS_PDE_GUI;
        load(fullfile(gui_dir,'PIETOOLS_PDE_Example_Transport_Eq_GUI.mat'));
        logval = app.loadData(data);
        if logval
            disp("Failed to load data object. Incorrect structure");
        end
    end
end


% %% Define dependent variables and system variable
% % 	PDE: x_{t} = -x_{s}                                 
% %   BCs: x(s=0) = 0 
% x = state('pde');
% pde = sys();
% %% Define equations
% eq_PDE = diff(x,t)==-diff(x,s); % 	PDE: x_{t} = x_{s}
% eq_BC = subs(x,s,0)==0;        %   BCs: x(s=0) = 0 
% %% addequations to pde system; set control inputs/observed inputs, if any
% pde = addequation(pde,[eq_PDE;eq_BC]);
% %% display pde to verify and convert to pie
% display(pde);
% pie = convert(pde,'pie');
% display(pie);
