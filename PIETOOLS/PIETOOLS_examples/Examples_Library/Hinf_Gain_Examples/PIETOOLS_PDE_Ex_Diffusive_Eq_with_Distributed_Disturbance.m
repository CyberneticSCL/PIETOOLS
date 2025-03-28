function [PDE_t,PDE_b] = PIETOOLS_PDE_Ex_Diffusive_Eq_with_Distributed_Disturbance(GUI,params)
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
% % Diffusion-reaction PDE from Valmorbida (2014):
% % PDE         x_{t} = Cm*x + (1/R)*x_{ss} + [s;s]*w(t)
% % with BCs    x(s=0) = 0
% %             x(s=1) = 0
% % And Output  z(t) = int(x1(t,s),s,0,1)
% %
% % Parameters Cm, R can be set.
% % gamma = 0.8102 for Cm = [1, 1.5; 5, 0.2], R = 2.6.
% %---------------------------------------------------------------------% %

% Determine the location of this example file <-- DO NOT MOVE THE FILE
loc = mfilename('fullpath');
root = fileparts(loc);

% Initialize variables
pvar s

%%% Executive Function:
evalin('base','Hinf_gain = 1;')
%%% evalin('base','Hinf_gain_dual = 1;')

Cm = [1, 1.5; 5, 0.2];      R = 2.6;
npars = length(params);
if npars~=0
    %%% Specify potential parameters
    for j=1:npars
        eval(params{j});
    end
end


% % % Construct the PDE.
ne = size(Cm,2);
on = eye(ne);     ze = zeros(ne);

%%% Batch input format
PDE_b.n0 = 0;   PDE_b.n1 = 0;   PDE_b.n2 = ne;   PDE_b.nw = 1;   PDE_b.nz = 1;
PDE_b.dom = [0,1];

PDE_b.A0 = Cm;
PDE_b.A2 = (1/R)*eye(ne);
PDE_b.B21 = s*ones(ne,1);
PDE_b.Ca1 = [1 zeros(1,ne-1)];

PDE_b.B = [on ze ze ze;
           ze on ze ze];


%%% Term-based input format
% Initialize 1D PDE state component.
PDE_t.x{1}.vars = s;   PDE_t.x{1}.dom = [0,1];
% Initialize finite-dimensional inputs and outputs.
PDE_t.w{1}.vars = [];  PDE_t.z{1}.vars = [];

% PDE: x_{t} = Cm * x + (1/R) * x_{ss}
PDE_t.x{1}.term{1}.x = 1;
PDE_t.x{1}.term{1}.D = [0; 2];
PDE_t.x{1}.term{1}.C = [Cm, (1/R)*eye(ne)];

% PDE: x_{t} = ... + [s;s] * w
PDE_t.x{1}.term{2}.w = 1;
PDE_t.x{1}.term{2}.C = s*ones(ne,1);

% Output: z = int_{0}^{1} x(s) ds
PDE_t.z{1}.term{1}.x = 1;
PDE_t.z{1}.term{1}.C = [1,zeros(1,ne-1)];
PDE_t.z{1}.term{1}.I{1} = PDE_t.x{1}.dom;

% BC 1: 0 = x(0)
PDE_t.BC{1}.term{1}.x = 1;
PDE_t.BC{1}.term{1}.loc = 0;

% BC 2: 0 = x(1)
PDE_t.BC{2}.term{1}.x = 1;
PDE_t.BC{2}.term{1}.loc = 1;


if GUI
    %%% Associated GUI save file
    app = PIETOOLS_PDE_GUI;
    load(fullfile(root,'PIETOOLS_PDE_Ex_Diffusive_Eq_with_Distributed_Disturbance_GUI.mat'));
    logval = app.loadData(data);
    if logval
        disp("Failed to load data object. Incorrect structure");
    end
end

end
% % [6] - 
% @inproceedings{valmorbida2014semi,
%   title={Semi-definite programming and functional inequalities for distributed parameter systems},
%   author={Valmorbida, Giorgio and Ahmadi, Mohamadreza and Papachristodoulou, Antonis},
%   booktitle={53rd IEEE conference on decision and control},
%   pages={4304--4309},
%   year={2014},
%   organization={IEEE}
% }