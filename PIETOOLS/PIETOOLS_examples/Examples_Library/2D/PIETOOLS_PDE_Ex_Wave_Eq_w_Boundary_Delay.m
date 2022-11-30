function [PDE_t] = PIETOOLS_PDE_Ex_Wave_Eq_w_Boundary_Delay(GUI,params)
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
% % Wave Equation with Delay in BCs (Xu 2006):
% % PDE         x_{tt}  = x_{s1s1}(t,s1);       s1 in [0,1]   
% % With BCs    x(t,0) = 0;   
% %             x_{s1}(t,1) = -k*(1-mu)*x_{t}(t,1) - k*mu*x_{t}(t-tau,1);
% %
% % Introduce:
% %     x1(t) = - k * (1-mu) * \dot{u}(1,t) - k*mu*\dot{u}(1,t-tau);
% %     x2(t,s1) = x(t,s1);
% %     x3(t,s1) = x_{t}(t,s1);
% %     x4(t,s2) = x(t-tau*s2,1);
% %     x5(t,s2) = x_{t}(t-tau*s2,1);
% % Then:
% % ODE:    x1_{t}(t)    = x2_{s1}(t,1);
% % PDE:    x2_{t}(t,s1) = x3(t,s1);                    s1 in (0,1)
% %         x3_{t}(t,s1) = x2_{s1s1}(t,s1);                
% %         x4_{t}(t,s2) = -(1/tau) x4_{s2}(t,s2);      s2 in (0,1)
% %         x5_{t}(t,s2) = -(1/tau) x5_{s2}(t,s2);
% % BCs:    x2(t,0) = 0;
% %         x3(t,0) = 0;
% %         x4(t,0) = x2(t,1);
% %         x5(t,0) = x3(t,1);
% %         x1(t) = - k * (1-mu) * x2(t,1) - k* mu * x4(t,1);
% %         x2_{s1}(t,1) = -k*(1-mu)*x3_{t}(t,1) - k*mu*x5_{t}(t,1);
% %
% % Parameters tau, k, and mu can be set.
% % Stable whenever mu < 1/2
% %---------------------------------------------------------------------% %

% Determine the location of this example file <-- DO NOT MOVE THE FILE
loc = mfilename('fullpath');
root = fileparts(loc);

% Initialize variables
pvar s theta s1 s2 theta1 theta2

%%% Executive Function:
evalin('base','stability = 1;');
%evalin('base','stability_dual = 1;');

a = 0;      b = 1;  
k = 1;
tau = 1;    mu = 0.5-1e-1;
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
PDE_t.x{1}.vars = [];
PDE_t.x{2}.vars = [s1,theta1];    PDE_t.x{2}.dom = [a,b];
PDE_t.x{3}.vars = [s1,theta1];    PDE_t.x{3}.dom = [a,b];
PDE_t.x{3}.diff = 2;    % x3 must be 2nd order differentiable wrt s1
PDE_t.x{4}.vars = [s2,theta2];    PDE_t.x{4}.dom = [0,1];
PDE_t.x{5}.vars = [s2,theta2];    PDE_t.x{5}.dom = [0,1];

% PDE: x1_{t}(t) = x2_{s1}(t,1)
PDE_t.x{1}.term{1}.x = 2;
PDE_t.x{1}.term{1}.D = 1;
PDE_t.x{1}.term{1}.loc = b;

% PDE: x2_{t}(t,s1) = x3(t,s1)
PDE_t.x{2}.term{1}.x = 3;

% PDE: x3_{t}(t,s1) = x2_{s1s1}(t,s1);
PDE_t.x{3}.term{1}.x = 2;
PDE_t.x{3}.term{1}.D = 2;

% PDE: x4_{t}(t,s2) = -(1/tau)*x4_{s2}(t,s2);
PDE_t.x{4}.term{1}.x = 4;
PDE_t.x{4}.term{1}.D = 1;
PDE_t.x{4}.term{1}.C = -1/tau;

% PDE: x5_{t}(t,s2) = -(1/tau)*x5_{s2}(t,s2);
PDE_t.x{5}.term{1}.x = 5;
PDE_t.x{5}.term{1}.D = 1;
PDE_t.x{5}.term{1}.C = -1/tau;


% BC1: x2(t,s1=a) = 0;              BC4: x3(t,s1=a) = 0;
PDE_t.BC{1}.term{1}.x = 2;          PDE_t.BC{4}.term{1}.x = 3;
PDE_t.BC{1}.term{1}.loc = a;        PDE_t.BC{4}.term{1}.loc = a;

% BC2: x1(t) = -k*(1-mu)*x2(t,s1=b) - k*mu*x4(t,s2=1)
PDE_t.BC{2}.term{1}.x = 1;  PDE_t.BC{2}.term{2}.x = 2;        PDE_t.BC{2}.term{3}.x = 4;
                            PDE_t.BC{2}.term{2}.C = k*(1-mu); PDE_t.BC{2}.term{3}.C = k*mu;
                            PDE_t.BC{2}.term{2}.loc = b;      PDE_t.BC{2}.term{3}.loc = 1;

% BC3: x2_{s1}(t,s1=b) = -k*(1-mu)*x3(t,s1=b) - k*mu*x5(t,s2=1)
PDE_t.BC{3}.term{1}.x = 2;    PDE_t.BC{3}.term{2}.x = 3;        PDE_t.BC{3}.term{3}.x = 5;
PDE_t.BC{3}.term{1}.D = 1;    PDE_t.BC{3}.term{2}.C = k*(1-mu); PDE_t.BC{3}.term{3}.C = k*mu;
PDE_t.BC{3}.term{1}.loc = b;  PDE_t.BC{3}.term{2}.loc = b;      PDE_t.BC{3}.term{3}.loc = 1;

% BC5: (Sewing condition) x4(t,s2=0) = x2(t,s1=b);
PDE_t.BC{5}.term{1}.x = 4;      PDE_t.BC{5}.term{2}.x = 2;
PDE_t.BC{5}.term{1}.loc = 0;    PDE_t.BC{5}.term{2}.loc = b;
                                PDE_t.BC{5}.term{2}.C = -1;

% BC6: (Sewing condition) x5(t,s2=0) = x3(t,s1=b);
PDE_t.BC{6}.term{1}.x = 5;      PDE_t.BC{6}.term{2}.x = 3;
PDE_t.BC{6}.term{1}.loc = 0;    PDE_t.BC{6}.term{2}.loc = b;
                                PDE_t.BC{6}.term{2}.C = -1;

                                
if GUI~=0
    disp('No GUI representation available for this system.')
end

end
% @article{XU_2006,
%      author = {Xu, Gen Qi and Yung, Siu Pang and Li, Leong Kwan},
%      title = {Stabilization of wave systems with input delay in the boundary control},
%      journal = {ESAIM: Control, Optimisation and Calculus of Variations},
%      pages = {770--785},
%      publisher = {EDP-Sciences},
%      volume = {12},
%      number = {4},
%      year = {2006},
%      doi = {10.1051/cocv:2006021},
%      zbl = {1105.35016},
%      mrnumber = {2266817},
%      language = {en},
%      url = {http://www.numdam.org/articles/10.1051/cocv:2006021/}
% }