function [PDE_t] = PIETOOLS_PDE_Ex_2D_KISS_Model(GUI,params)
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
% % KISS Model from Holmes (1994):
% % PDE         x_{t}   = lam*x + c1*x_{s1s1} + c2*x_{s2s2}
% % With BCs    x(s1=0) = 0;   x(s2=0) = 0;
% %             x(s1=1) = 0;   x(s2=1) = 0;
% %
% % Parameters c1, c2 and lam can be set.
% % Stable for lam <= 2*pi^2 when c1=c2=1.
% %---------------------------------------------------------------------% %

% Determine the location of this example file <-- DO NOT MOVE THE FILE
loc = mfilename('fullpath');
root = fileparts(loc);

% Initialize variables
pvar s theta s1 s2 theta1 theta2

%%% Executive Function:
evalin('base','stability = 1;');
%evalin('base','stability_dual = 1;');

c1 = 1; c2 = 1; lam = 19;   ne = 1;
npars = length(params);
if npars~=0
    %%% Specify potential parameters
    for j=1:npars
        eval(params{j});
    end
end


% % % Construct the PDE.
%%% Term-based input format
PDE_t.x{1}.vars = [s1;s2];   PDE_t.x{1}.dom = [0,1;0,1];

% PDE: x_{t} = [lam, c1, c2] * [x; x_{s1s1}; x_{s2s2}]
PDE_t.x{1}.term{1}.D = [0,0; 2,0; 0,2];
PDE_t.x{1}.term{1}.C = [lam*eye(ne), c1*eye(ne), c2*eye(ne)];

% BC1: 0 = x(s1,0)                     % BC3: 0 = x(s1,1)
PDE_t.BC{1}.term{1}.loc = [s1,0];      PDE_t.BC{3}.term{1}.loc = [s1,1];
% BC2: 0 = x(0,s2)                     % BC4: 0 = x(1,s2)
PDE_t.BC{2}.term{1}.loc = [0,s2];      PDE_t.BC{4}.term{1}.loc = [1,s2];


if GUI~=0
    disp('No GUI representation available for this system.')
end

end
% E. E. Holmes, M. Lewis, J. Banks, and R. R. Veit, “Partial differential
% equations in ecology: Spatial interactions and population dynamics,”
% Ecology, vol. 75, pp. 17–29, 1994.