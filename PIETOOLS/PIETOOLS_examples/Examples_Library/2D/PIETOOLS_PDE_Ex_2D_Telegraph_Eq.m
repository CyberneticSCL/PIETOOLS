function [PDE_t] = PIETOOLS_PDE_Ex_2D_Telegraph_Eq(GUI,params)
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
% % % Telegraph equation (Holmes, 1994):
% % % PDE         x_{t}   = -1/(2*lam) * x_{tt} + (C^2/2*lam)*(x_{s1s1} + x_{s2s2})
% % % With BCs    x(s1=0) = 0;   x(s2=0) = 0;
% % %             x(s1=1) = 0;   x(s2=1) = 0;
% % %
% % % Use states x1 = x, x2 = x_{t}, then:
% % %
% % % PDE         x1_{t}   = x2
% % %             x2_{t}   = -2*lam * x2 + C^2*(x1_{s1s1} + x1_{s2s2})
% % % With BCs    x1(s1=0) = 0;   x1(s2=0) = 0;
% % %             x1(s1=1) = 0;   x1(s2=1) = 0;
% % %             x2(s1=0) = 0;   x2(s2=0) = 0;
% % %             x2(s1=1) = 0;   x2(s2=1) = 0;
% %
% % Parameters C, lam can be set.
% %---------------------------------------------------------------------% %

% Determine the location of this example file <-- DO NOT MOVE THE FILE
loc = mfilename('fullpath');
root = fileparts(loc);

% Initialize variables
pvar s theta s1 s2 theta1 theta2

%%% Executive Function:
evalin('base','stability = 1;');
%evalin('base','stability_dual = 1;');

C = 1;  lam = 1;    ne = 1;
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
PDE_t.x{2}.vars = [s1;s2];   PDE_t.x{2}.dom = [0,1;0,1];
PDE_t.x{2}.diff = [2,2];    % x2 must be 2nd order differentiable wrt both s1 and s2

% PDE: x1_{t} = x2;
PDE_t.x{1}.term{1}.x = 2;

% PDE: x2_{t} = -2*lam*x2
PDE_t.x{2}.term{1}.x = 2;
PDE_t.x{2}.term{1}.C = -2*lam;

% PDE: x2_{t} = ... + [C^2, C^2] * [x1_{s1s1}; x1_{s2s2}];
PDE_t.x{2}.term{2}.x = 1;
PDE_t.x{2}.term{2}.D = [2,0; 0,2];
PDE_t.x{2}.term{2}.C = [C^2 * eye(ne), C^2 * eye(ne)];

% BC1: 0 = x1(s1,0)                    % BC3: 0 = x1(s1,1)
PDE_t.BC{1}.term{1}.x = 1;             PDE_t.BC{3}.term{1}.x = 1;
PDE_t.BC{1}.term{1}.loc = [s1,0];      PDE_t.BC{3}.term{1}.loc = [s1,1];
% BC2: 0 = x1(0,s2)                    % BC4: 0 = x1(1,s2)
PDE_t.BC{2}.term{1}.x = 1;             PDE_t.BC{4}.term{1}.x = 1;
PDE_t.BC{2}.term{1}.loc = [0,s2];      PDE_t.BC{4}.term{1}.loc = [1,s2];

% BC5: 0 = x1(s1,0)                    % BC7: 0 = x1(s1,1)
PDE_t.BC{5}.term{1}.x = 2;             PDE_t.BC{7}.term{1}.x = 2;
PDE_t.BC{5}.term{1}.loc = [s1,0];      PDE_t.BC{7}.term{1}.loc = [s1,1];
% BC6: 0 = x1(0,s2)                    % BC8: 0 = x1(1,s2)
PDE_t.BC{6}.term{1}.x = 2;             PDE_t.BC{8}.term{1}.x = 2;
PDE_t.BC{6}.term{1}.loc = [0,s2];      PDE_t.BC{8}.term{1}.loc = [1,s2];


if GUI~=0
    disp('No GUI representation available for this system.')
end

end
% E. E. Holmes, M. Lewis, J. Banks, and R. R. Veit, “Partial differential
% equations in ecology: Spatial interactions and population dynamics,”
% Ecology, vol. 75, pp. 17–29, 1994.