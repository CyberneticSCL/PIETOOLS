function [PDE_t] = PIETOOLS_PDE_Ex_Coupled_ODE_Heat_Eq_w_Delay(GUI,params)
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
% % ODE coupled to heat equation with delay (Kang, 2017):
% % ODE         X_{t} = A*X(t) + A1*X(t-tau) + B*x(t,s1=0)
% % PDE         x_{t} = x_{s1s1} + a*x + a2*x(t-tau);     s1 in [0,1]
% % With BCs    x_{s1}(t,s1=0) = 0;   
% %             x(t,s1=1) = u(t);   or x_{s1}(t,s1=1) = u(t);
% %
% % Let x1(t)=X(t) and x2(t,s1)=x(t,s1). Let further x3(t,s2) and x4(t,s1,s2)
% % be such that:
% %     x3(t,s2=1+tau)    = x1(t)     with x3_{t} = x3_{s2}   s2 in [1,1+tau]
% %     x4(t,s1,s2=1+tau) = x2(t,s1)  with x4_{t} = x4_{s2}
% % Then we can write the system as
% %
% % ODE:    x1_{t} = A*x1(t) + A1*x3(t,1) + B*x2(t,s1=0);
% % PDE:    x2_{t} = x2_{s1s1} + a*x2 + a2*x4(t,s1,s2=1);
% %         x3_{t} = x3_{s2};
% %         x4_{t} = x4_{s2};
% % BCs:    x2_{s1}(t,s1=0) = 0;
% %         x2(t,s1=1) = u(t);              or x2_{s1}(t,s1=1) = u(t);
% %         x3(t,s2=1+tau) = x1(t);
% %         x4(t,s1,s2=1+tau) = x2(t,s1);
% %         x4_{s1}(t,s1=0,s2) = 0;
% %         x4(t,s1=1,s2) = u(t);
% %
% % Parameters tau, A, A1, B, a and a2 can be set.
% %---------------------------------------------------------------------% %

% Determine the location of this example file <-- DO NOT MOVE THE FILE
loc = mfilename('fullpath');
root = fileparts(loc);

% Initialize variables
pvar s theta s1 s2 theta1 theta2

%%% Executive Function:
evalin('base','stability = 1;');
%evalin('base','stability_dual = 1;');

ne = 1; tau = 1; A = 1; A1 = 0.5; B = -1; a = 1; a2 = 0.5;
npars = length(params);
if npars~=0
    %%% Specify potential parameters
    for j=1:npars
        eval(params{j});
    end
end


% % % Construct the PDE.
% Initialize the state variables and input.
PDE_t.x{1}.vars = [];
PDE_t.x{2}.vars = s1;       PDE_t.x{2}.dom = [0,1];
PDE_t.x{3}.vars = s2;       PDE_t.x{3}.dom = [1,1+tau];
PDE_t.x{4}.vars = [s1;s2];  PDE_t.x{4}.dom = [0,1; 1,1+tau];
PDE_t.u{1}.vars = [];   % Input does not vary in space.

% ODE: x1_{t} = A*x1(t)
PDE_t.x{1}.term{1}.x = 1;
PDE_t.x{1}.term{1}.C = A;
% ODE: x1_{t} = ... + A1*x3(t,1)
PDE_t.x{1}.term{2}.x = 3;
PDE_t.x{1}.term{2}.loc = 1;
PDE_t.x{1}.term{2}.C = A1;
% ODE: x1_{t} = ... + B*x2(t,s1=0);
PDE_t.x{1}.term{3}.x = 2;
PDE_t.x{1}.term{3}.loc = 0;
PDE_t.x{1}.term{3}.C = B;

% PDE: x2_{t} = x2_{s1s1} + a*x2
PDE_t.x{2}.term{1}.x = [2; 2];
PDE_t.x{2}.term{1}.D = [2; 0];
PDE_t.x{2}.term{1}.C = [1, a];
% PDE: x2_{t} = ... + a2*x4(t,s1,s2=1);
PDE_t.x{2}.term{2}.x = 4;
PDE_t.x{2}.term{2}.loc = [s1,1];
PDE_t.x{2}.term{2}.C = a2;

% PDE: x3_{t} = x3_{s2};
PDE_t.x{3}.term{1}.x = 3;
PDE_t.x{3}.term{1}.D = 1;

% PDE: x4_{t} = x4_{s2};
PDE_t.x{4}.term{1}.x = 4;
PDE_t.x{4}.term{1}.D = [0,1];


% BC1: 0 = x2_{s1}(t,s1=0);
PDE_t.BC{1}.term{1}.x = 2;
PDE_t.BC{1}.term{1}.D = 1;
PDE_t.BC{1}.term{1}.loc = 0;

% BC2: 0 = x2(t,s1=1) - u(t)
PDE_t.BC{2}.term{1}.x = 2;          PDE_t.BC{2}.term{2}.u = 1;
PDE_t.BC{2}.term{1}.loc = 1;        PDE_t.BC{2}.term{2}.C = -1;
PDE_t.BC{2}.term{1}.D = 0;

% BC3: 0 = x3(t,s2=1+tau) - x1(t);
PDE_t.BC{3}.term{1}.x = 3;          PDE_t.BC{3}.term{2}.x = 1;
PDE_t.BC{3}.term{1}.loc = 1+tau;    PDE_t.BC{3}.term{2}.C = -1;

% BC4: 0 = x4(t,s1,s2=1+tau) - x2(t,s1);
PDE_t.BC{4}.term{1}.x = 4;              PDE_t.BC{4}.term{2}.x = 2;
PDE_t.BC{4}.term{1}.loc = [s1,1+tau];   PDE_t.BC{4}.term{2}.loc = s1;          
                                        PDE_t.BC{4}.term{2}.C = -1; 
                                        
% BC5: 0 = x4_{s1}(t,s1=0,s2);
PDE_t.BC{5}.term{1}.x = 4;
PDE_t.BC{5}.term{1}.D = [1,0];
PDE_t.BC{5}.term{1}.loc = [0,s2];

% BC6: 0 = x4(t,s1=1,s2) - u(t)
PDE_t.BC{6}.term{1}.x = 4;          PDE_t.BC{6}.term{2}.u = 1;
PDE_t.BC{6}.term{1}.loc = [1,s2];   PDE_t.BC{6}.term{2}.C = -1;
PDE_t.BC{6}.term{1}.D = [0,0];              
                                        
                                        
if GUI~=0
    disp('No GUI representation available for this system.')
end

end
% @article{kang2017boundary,
%   title={Boundary control of delayed ODE--heat cascade under actuator saturation},
%   author={Kang, Wen and Fridman, Emilia},
%   journal={Automatica},
%   volume={83},
%   pages={252--261},
%   year={2017},
%   publisher={Elsevier}
% }