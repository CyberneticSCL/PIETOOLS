function [PDE_t] = PIETOOLS_PDE_Ex_2D_Heat_Eq_with_ODE(GUI,params)
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
% % 2D heat equation coupled to ODE
% % ODE         x1_{t}  = (A+BK)*x1 + B*x2(s1=0,s2=0)
% % PDE         x2_{t}  = c1*x_{(2,0)} + c2*x_{(0,2)}
% % With BCs    x2_{(1,0)}(s1=0) = 0;     x2(s1=1) = 0;
% %             x2_{(0,1)}(s2=0) = 0;     x2(s2=1) = 0;
% % 
% % Parameters A, B and K can be set.
% % Stable if A+BK is Hurwitz.
% %---------------------------------------------------------------------% %

% Determine the location of this example file <-- DO NOT MOVE THE FILE
loc = mfilename('fullpath');
root = fileparts(loc);

% Initialize variables
pvar s theta s1 s2 theta1 theta2

%%% Executive Function:
evalin('base','stability = 1;');
%evalin('base','stability_dual = 1;');

c1 = 1;  c2 = 1;   k = -2;
nx = 1;    ne = 1;
A = eye(nx);    B = eye(nx,ne);  K = k*eye(ne,nx);
npars = length(params);
if npars~=0
    %%% Specify potential parameters
    for j=1:npars
        eval(params{j});
    end
end


% % % Construct the PDE.
ne = size(K,1);

%%% Term-based input format
% Initialize an ODE state component.
PDE_t.x{1}.vars = [];
% Initialize a 2D PDE state component.
PDE_t.x{2}.vars = [s1;s2];   PDE_t.x{2}.dom = [0,1;0,1];

% ODE: x1_{t} = (A+B*K)*x1
PDE_t.x{1}.term{1}.x = 1;
PDE_t.x{1}.term{1}.C = A+B*K;

% ODE: x1_{t} = ... + B*x2(s1=0,s2=0)
PDE_t.x{1}.term{2}.x =2;
PDE_t.x{1}.term{2}.loc = [0,0];
PDE_t.x{1}.term{2}.C = B;

% PDE: x2_{t} = [c1, c2] * [x2_{s1s1}; x2_{s2s2}]
PDE_t.x{2}.term{1}.x = 2;
PDE_t.x{2}.term{1}.D = [2,0; 0,2];
PDE_t.x{2}.term{1}.C = [c1*eye(ne), c2*eye(ne)];

% BC1: 0 = x2_{s2}(s1,0)
PDE_t.BC{1}.term{1}.x = 2;
PDE_t.BC{1}.term{1}.D = [0,1];
PDE_t.BC{1}.term{1}.loc = [s1,0];
% BC2: 0 = x2_{s1}(0,s2)
PDE_t.BC{2}.term{1}.x = 2;
PDE_t.BC{2}.term{1}.D = [1,0];
PDE_t.BC{2}.term{1}.loc = [0,s2];
% BC3: 0 = x2(s1,1)
PDE_t.BC{3}.term{1}.x = 2;
PDE_t.BC{3}.term{1}.loc = [s1,1];
% BC4: 0 = x2(1,s2)
PDE_t.BC{4}.term{1}.x = 2;
PDE_t.BC{4}.term{1}.loc = [1,s2];


if GUI~=0
    disp('No GUI representation available for this system.')
end

end