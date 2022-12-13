function [PDE_t] = PIETOOLS_PDE_Ex_Heat_Eq_w_Delayed_Boundary_Input(GUI,params)
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
% % Diffusive Equation with Delayed Boundary Input (Kristic 2009):
% % PDE         x_{t}     = x_{s1s1} + lam * x;     s1 in [0,1]
% % With BCs    x(t,s1=0) = 0;   
% %             x(t,s1=1) = u(t-tau);
% %
% % Introduce transport equation:
% %     u(t) = x2(t,1+tau),   x2_{t} = x2_{s2}        s2 in [1,1+tau]
% % Then
% %
% % PDE:    x1_{t}         = x1_{s1s1} + lam * x1;
% %         x2_{t}         = x2_{s2};
% % BCs:    x1(t,s1=0)     = 0;
% %         x1(t,s1=1)     = x2(t,s2=1);
% %         x2(t,s2=1+tau) = u(t);
% %
% % Parameters tau and lam can be set.
% %---------------------------------------------------------------------% %

% Determine the location of this example file <-- DO NOT MOVE THE FILE
loc = mfilename('fullpath');
root = fileparts(loc);

% Initialize variables
pvar s theta s1 s2 theta1 theta2

%%% Executive Function:
evalin('base','stability = 1;');
%evalin('base','stability_dual = 1;');

ne = 1; tau = 1;  lam = 1;
npars = length(params);
if npars~=0
    %%% Specify potential parameters
    for j=1:npars
        eval(params{j});
    end
end


% % % Construct the PDE.
expand_delays = false;
if ~expand_delays
    %%% Term-baed input format with delay.
    PDE_t.x{1}.vars = s1;       PDE_t.x{1}.dom = [0,1];
    PDE_t.u{1}.vars = [];   % Input does not vary in space.
    
    % PDE: x1_{t} = x1_{s1s1} + lam * x1;
    PDE_t.x{1}.term{1}.x = [1;1];
    PDE_t.x{1}.term{1}.D = [2; 0];
    PDE_t.x{1}.term{1}.C = [eye(ne), lam*eye(ne)];
    
    % BC1: 0 = x1(t,s1=0)
    PDE_t.BC{1}.term{1}.x = 1;
    PDE_t.BC{1}.term{1}.loc = 0;
    % BC2: 0 = x1(t,s1=1) - u(t-tau)
    PDE_t.BC{2}.term{1}.x = 1;          PDE_t.BC{2}.term{2}.u = 1;
    PDE_t.BC{2}.term{1}.loc = 1;        PDE_t.BC{2}.term{2}.delay = tau;
                                        PDE_t.BC{2}.term{2}.C = -1;
    
else
    %%% Term-based input format
    % Initialize the state variables and input.
    PDE_t.x{1}.vars = s1;       PDE_t.x{1}.dom = [0,1];
    PDE_t.x{2}.vars = s2;       PDE_t.x{2}.dom = [1,1+tau];
    PDE_t.u{1}.vars = [];   % Input does not vary in space.

    % PDE: x1_{t} = x1_{s1s1} + lam * x1;
    PDE_t.x{1}.term{1}.x = [1;1];
    PDE_t.x{1}.term{1}.D = [2; 0];
    PDE_t.x{1}.term{1}.C = [eye(ne), lam*eye(ne)];

    % PDE: x2_{t} = x2_{s2}
    PDE_t.x{2}.term{1}.x = 2;
    PDE_t.x{2}.term{1}.D = 1;

    % BC1: 0 = x1(s1=0)
    PDE_t.BC{1}.term{1}.x = 1;
    PDE_t.BC{1}.term{1}.loc = 0;
    % BC2: 0 = x1(s1=1) - x2(s2=1)
    PDE_t.BC{2}.term{1}.x = 1;          PDE_t.BC{2}.term{2}.x = 2;
    PDE_t.BC{2}.term{1}.loc = 1;        PDE_t.BC{2}.term{2}.loc = 1;
                                        PDE_t.BC{2}.term{2}.C = -1;
    % BC3: 0 = x2(s2=1+D) - u;
    PDE_t.BC{3}.term{1}.x = 2;          PDE_t.BC{3}.term{2}.u = 1;
    PDE_t.BC{3}.term{1}.loc = 1+tau;      PDE_t.BC{3}.term{2}.C = -1;
end

if GUI~=0
    disp('No GUI representation available for this system.')
end

end
% @article{krstic2009control,
%   title={Control of an unstable reaction--diffusion PDE with long input delay},
%   author={Krstic, Miroslav},
%   journal={Systems \& Control Letters},
%   volume={58},
%   number={10-11},
%   pages={773--782},
%   year={2009},
%   publisher={Elsevier}
% }