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
% % PDE         x_{tt}  = x_{ss}(t,s);       s1 in [0,1]   
% % With BCs    x(t,0) = 0;   
% %             x_{s}(t,1) = -k*(1-mu)*x_{t}(t,1) - k*mu*x_{t}(t-tau,1);
% %
% % Introduce:
% %     x1(t)   = -k*(1-mu)*x(t,1) -k*mu*x(t-tau,1);
% %     x2(t,s) = x(t,s);
% %     x3(t,s) = x_{t}(t,s);
% %     x4(t,s) = x(t-tau*s,1);
% %     x5(t,s) = x_{t}(t-tau*s,1);
% % Then:
% % ODE:    x1_{t}(t)   = x2_{s}(t,1);
% % PDE:    x2_{t}(t,s) = x3(t,s);                    s in (0,1)
% %         x3_{t}(t,s) = x2_{ss}(t,s);                
% %         x4_{t}(t,s) = -(1/tau) x4_{s}(t,s);     
% %         x5_{t}(t,s) = -(1/tau) x5_{s}(t,s);
% % BCs:    x2(t,0) = 0;
% %         x3(t,0) = 0;
% %         x4(t,0) = x2(t,1);
% %         x5(t,0) = x3(t,1);
% %         x1(t)   = - k*(1-mu)*x2(t,1) - k*mu*x4(t,1);
% %         x2_{s}(t,1) = - k*(1-mu)*x3(t,1) - k*mu*x5(t,1);
% %
% % Parameters tau, k, and mu can be set.
% % Stable whenever mu < 1/2
% %---------------------------------------------------------------------% %

% Determine the location of this example file <-- DO NOT MOVE THE FILE
loc = mfilename('fullpath');
root = fileparts(loc);

% Initialize variables
pvar s

%%% Executive Function:
evalin('base','stability = 1;');
%evalin('base','stability_dual = 1;');

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
%%% pde_var input format
clear stateNameGenerator
x1 = pde_var();
x2 = pde_var(s,[0,1]);
x3 = pde_var(s,[0,1],2);       % x3 must be 2nd order differentiable wrt s1
x4 = pde_var(s,[0,1]);
x5 = pde_var(s,[0,1]);
PDE_t = [diff(x1,'t')==subs(diff(x2,s),s,1);
         diff(x2,'t')==x3;
         diff(x3,'t')==diff(x2,s,2);
         diff(x4,'t')==-(1/tau)*diff(x4,s);
         diff(x5,'t')==-(1/tau)*diff(x5,s);
         subs([x2;x3],s,0)==0;
         x1==-k*(1-mu)*subs(x2,s,1)-k*mu*subs(x4,s,1);
         subs(diff(x2,s),s,1)==-k*(1-mu)*subs(x3,s,1)-k*mu*subs(x5,s,1);
         subs(x4,s,0)==subs(x2,s,1);
         subs(x5,s,0)==subs(x3,s,1)];


% %%% Term-based input format
% % Initialize the state variables and input.
% pde_struct PDE_t;
% PDE_t.x{1}.vars = [];
% PDE_t.x{2}.vars = s;    PDE_t.x{2}.dom = [0,1];
% PDE_t.x{3}.vars = s;    PDE_t.x{3}.dom = [0,1];
% PDE_t.x{3}.diff = 2;    % x3 must be 2nd order differentiable wrt s1
% PDE_t.x{4}.vars = s;    PDE_t.x{4}.dom = [0,1];
% PDE_t.x{5}.vars = s;    PDE_t.x{5}.dom = [0,1];
% 
% % PDE: x1_{t}(t) = x2_{s}(t,1)
% PDE_t.x{1}.term{1}.x = 2;
% PDE_t.x{1}.term{1}.D = 1;
% PDE_t.x{1}.term{1}.loc = 1;
% 
% % PDE: x2_{t}(t,s) = x3(t,s)
% PDE_t.x{2}.term{1}.x = 3;
% 
% % PDE: x3_{t}(t,s) = x2_{ss}(t,s);
% PDE_t.x{3}.term{1}.x = 2;
% PDE_t.x{3}.term{1}.D = 2;
% 
% % PDE: x4_{t}(t,s) = -(1/tau)*x4_{s}(t,s);
% PDE_t.x{4}.term{1}.x = 4;
% PDE_t.x{4}.term{1}.D = 1;
% PDE_t.x{4}.term{1}.C = -1/tau;
% 
% % PDE: x5_{t}(t,s) = -(1/tau)*x5_{s}(t,s);
% PDE_t.x{5}.term{1}.x = 5;
% PDE_t.x{5}.term{1}.D = 1;
% PDE_t.x{5}.term{1}.C = -1/tau;
% 
% 
% % BC1: x2(t,s=0) = 0;              BC4: x3(t,s=0) = 0;
% PDE_t.BC{1}.term{1}.x = 2;          PDE_t.BC{4}.term{1}.x = 3;
% PDE_t.BC{1}.term{1}.loc = 0;        PDE_t.BC{4}.term{1}.loc = 0;
% 
% % BC2: x1(t) = -k*(1-mu)*x2(t,s=1) - k*mu*x4(t,s=1)
% PDE_t.BC{2}.term{1}.x = 1;  PDE_t.BC{2}.term{2}.x = 2;        PDE_t.BC{2}.term{3}.x = 4;
%                             PDE_t.BC{2}.term{2}.C = k*(1-mu); PDE_t.BC{2}.term{3}.C = k*mu;
%                             PDE_t.BC{2}.term{2}.loc = 1;      PDE_t.BC{2}.term{3}.loc = 1;
% 
% % BC3: x2_{s1}(t,s=1) = -k*(1-mu)*x3(t,s=1) - k*mu*x5(t,s=1)
% PDE_t.BC{3}.term{1}.x = 2;    PDE_t.BC{3}.term{2}.x = 3;        PDE_t.BC{3}.term{3}.x = 5;
% PDE_t.BC{3}.term{1}.D = 1;    PDE_t.BC{3}.term{2}.C = k*(1-mu); PDE_t.BC{3}.term{3}.C = k*mu;
% PDE_t.BC{3}.term{1}.loc = 1;  PDE_t.BC{3}.term{2}.loc = 1;      PDE_t.BC{3}.term{3}.loc = 1;
% 
% % BC5: (Sewing condition) x4(t,s=0) = x2(t,s=1);
% PDE_t.BC{5}.term{1}.x = 4;      PDE_t.BC{5}.term{2}.x = 2;
% PDE_t.BC{5}.term{1}.loc = 0;    PDE_t.BC{5}.term{2}.loc = 1;
%                                 PDE_t.BC{5}.term{2}.C = -1;
% 
% % BC6: (Sewing condition) x5(t,s=0) = x3(t,s=1);
% PDE_t.BC{6}.term{1}.x = 5;      PDE_t.BC{6}.term{2}.x = 3;
% PDE_t.BC{6}.term{1}.loc = 0;    PDE_t.BC{6}.term{2}.loc = 1;
%                                 PDE_t.BC{6}.term{2}.C = -1;

                                
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