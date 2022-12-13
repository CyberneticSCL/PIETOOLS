function [PDE_t,PDE_b] = PIETOOLS_PDE_Ex_Diffusive_Eq_with_Distributed_Disturbance_2(GUI,params)
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
% % Diffusion-reaction PDE from Shivakumar (2019):
% % PDE         x_{t} = Cm(s)*x + (1/R)*x_{ss} + [s;s;s]*w(t)
% % with BCs    x(s=0) = 0
% %             x(s=1) = 0
% % And Output  z(t) = int(x(t,s),s,0,1)
% %
% % Parameters Cm, R can be set.
% % gamma =  4.23, for R = 21-1e-3 and Cm = [0,0,0; s,0,0; s^2,-s^3,0];
% %---------------------------------------------------------------------% %

% Determine the location of this example file <-- DO NOT MOVE THE FILE
loc = mfilename('fullpath');
root = fileparts(loc);

% Initialize variables
pvar s theta s1 s2 theta1 theta2

%%% Executive Function:
evalin('base','Hinf_gain = 1;')
%%% evalin('base','Hinf_gain_dual = 1;')

Cm = [0,0,0; s,0,0; s^2,-s^3,0];    R = (21-1e-3); 
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
PDE_b.n0 = 0;   PDE_b.n1 = 0;   PDE_b.n2 = ne;   PDE_b.nw = 1;   PDE_b.nz = ne;
PDE_b.dom = [0,1];

PDE_b.A0 = Cm;
PDE_b.A2= (1/R)*on;
PDE_b.B21 = s*ones(ne,1);
PDE_b.Ca1 = eye(PDE_b.n2);

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

% PDE: x_{t} = ... + [s;s;s] * w
PDE_t.x{1}.term{2}.w = 1;
PDE_t.x{1}.term{2}.C = s*ones(ne,1);

% Output: z = int_{0}^{1} x(s) ds
PDE_t.z{1}.term{1}.x = 1;
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
    load(fullfile(root,'PIETOOLS_PDE_Ex_Diffusive_Eq_with_Distributed_Disturbance_2_GUI.mat'));
    logval = app.loadData(data);
    if logval
        disp("Failed to load data object. Incorrect structure");
    end
end

end
% @inproceedings{shivakumar2019computing,
%   title={Computing input-ouput properties of coupled linear pde systems},
%   author={Shivakumar, Sachin and Peet, Matthew M},
%   booktitle={2019 American Control Conference (ACC)},
%   pages={606--613},
%   year={2019},
%   organization={IEEE}
% }