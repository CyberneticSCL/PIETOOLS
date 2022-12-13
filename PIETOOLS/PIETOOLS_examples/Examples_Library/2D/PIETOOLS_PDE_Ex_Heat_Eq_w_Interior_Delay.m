function [PDE_t] = PIETOOLS_PDE_Ex_Heat_Eq_w_Interior_Delay(GUI,params)
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
% % Diffusive Equation with Delay in Dynamics (Caliskan 2009):
% % PDE         x_{t}  = c*x_{s1s1}(t,s1) + a0*x(t,s1) - a1*x(t-tau,s1);    s1 in [0,pi]   
% % With BCs    x(t,s1=0) = 0;   
% %             x(t,s1=pi) = 0;
% %
% % Let x1=x, and x2(t,s1,s2)=x1(t-s2,s1), introducing transport equation:
% %         x1(t) = x2(t,0),   x2_{t} = -x2_{s2}        s2 in [0,tau]
% % Then
% %
% % PDE:    x1_{t}         = c*x1_{s1s1}(t,s1) + a0*x1(t,s1) - a1*x2(t,s1,s2=tau);
% %         x2_{t}         = -x2_{s2};
% % BCs:    x1(t,s1=0)     = 0;
% %         x1(t,s1=pi)    = 0;
% %         x2(t,s1,s2=0)  = x1(t,s1);
% %         x2(t,s1=0,s2)  = 0;
% %         x2(t,s1=pi,s2) = 0;
% %
% % Parameters tau, c, a0 and a1 can be set.
% % For c = 1, a0 = 1.9, and a1 = 1, stable for tau<1.0347.
% %---------------------------------------------------------------------% %

% Determine the location of this example file <-- DO NOT MOVE THE FILE
loc = mfilename('fullpath');
root = fileparts(loc);

% Initialize variables
pvar s theta s1 s2 theta1 theta2

%%% Executive Function:
evalin('base','stability = 1;');
%evalin('base','stability_dual = 1;');

c = 1;      a0 = 1.9;       a1 = 1;
ne = 1;     tau = 1; 
a = 0;      b = pi;
npars = length(params);
if npars~=0
    %%% Specify potential parameters
    for j=1:npars
        eval(params{j});
    end
end


% % % Construct the PDE.
%%% Term-based input format
% Initialize the state variables and input.
pde_struct PDE_t;
PDE_t.x{1}.vars = s1;           PDE_t.x{1}.dom = [a,b];
PDE_t.x{2}.vars = [s1;s2];      PDE_t.x{2}.dom = [a,b;0,tau];
PDE_t.x{2}.diff = [2,1];    % x2 must be 2nd order differentiable wrt s1, and 1st order wrt s2

% PDE: x1_{t}(t,s1) = c*x1_{s1s1}(t,s1)
PDE_t.x{1}.term{1}.x = 1;
PDE_t.x{1}.term{1}.D = 2;
PDE_t.x{1}.term{1}.C = c;

% PDE: x1_{t}(t,s1) = ... + a0*x1(t,s1)
PDE_t.x{1}.term{2}.x = 1;
PDE_t.x{1}.term{2}.C = a0;

% PDE: x1_{t}(t,s1) = ... + a1*x2(t,s1,tau)
PDE_t.x{1}.term{3}.x = 2;
PDE_t.x{1}.term{3}.loc = [s1,tau];
PDE_t.x{1}.term{3}.C = -a1;

% PDE: x2_{t}(t,s1,s2) = -x2_{s2}(t,s1,s2)
PDE_t.x{2}.term{1}.x = 2;
PDE_t.x{2}.term{1}.D = [0,1];
PDE_t.x{2}.term{1}.C = -1;


% BC1: x1(s1=a) = 0;                BC2: x1(s1=b)=0;
PDE_t.BC{1}.term{1}.x = 1;          PDE_t.BC{2}.term{1}.x = 1;
PDE_t.BC{1}.term{1}.loc = a;        PDE_t.BC{2}.term{1}.loc = b;

% BC3: (Sewing condition) x2(s1,s2=0) = x1(s1)
PDE_t.BC{3}.term{1}.x = 2;          PDE_t.BC{3}.term{2}.x = 1;
PDE_t.BC{3}.term{1}.loc = [s1,0];   PDE_t.BC{3}.term{2}.C = -1;

% BC4: x2(s1=a,s2) = 0;             BC5: x2(s1=b,s2) = 0;
PDE_t.BC{4}.term{1}.x = 2;          PDE_t.BC{5}.term{1}.x = 2;
PDE_t.BC{4}.term{1}.loc = [a,s2];   PDE_t.BC{5}.term{1}.loc = [b,s2];


if GUI~=0
    disp('No GUI representation available for this system.')
end

end
% @article{CALISKAN2009,
% title = {Stability Analysis of the Heat Equation with Time-Delayed Feedback},
% journal = {IFAC Proceedings Volumes},
% volume = {42},
% number = {6},
% pages = {220-224},
% year = {2009},
% note = {6th IFAC Symposium on Robust Control Design},
% issn = {1474-6670},
% doi = {https://doi.org/10.3182/20090616-3-IL-2002.00038},
% url = {https://www.sciencedirect.com/science/article/pii/S1474667015404057},
% author = {Sina Yamaç çalişkan and Hitay özbay},
% }