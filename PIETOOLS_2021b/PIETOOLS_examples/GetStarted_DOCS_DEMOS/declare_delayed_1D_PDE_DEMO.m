% % This script illustrates the process of declaring a 1D ODE-PDE with
% % delayed state variables in PIETOOLS. Delayed PDEs cannot be directly
% % declared in PIETOOLS, but we can represent such systems as 2D PDEs,
% % which can be declared.
% %
% % As an example, we implement an ODE coupled to a heat equation with
% % delay, from (Kang, 2017) (see reference below):
% %
% % ODE         X_{t} = A*X(t) + A1*X(t-tau) + B*x(t,s1=0)
% % PDE         x_{t} = x_{s1s1} + a*x + a2*x(t-tau);     s1 in [0,1]
% % With BCs    x_{s1}(t,s1=0) = 0;   
% %             x(t,s1=1) = u(t);   or x_{s1}(t,s1=1) = u(t);
% %
% % To declare this system, we let x1(t)=X(t) and x2(t,s1)=x(t). We 
% % implement the delay by declaring x3(t,s2) and x4(t,s1,s2) such that:
% % 
% %     x3(t,s2=1+tau)    = x1(t)     with x3_{t} = x3_{s2}   s2 in [1,1+tau]
% %     x4(t,s1,s2=1+tau) = x2(t,s1)  with x4_{t} = x4_{s2}
% % 
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
% %
% % =================================================================== % %
% @article{kang2017boundary,
%   title={Boundary control of delayed ODE--heat cascade under actuator saturation},
%   author={Kang, Wen and Fridman, Emilia},
%   journal={Automatica},
%   volume={83},
%   pages={252--261},
%   year={2017},
%   publisher={Elsevier}
% }
% %---------------------------------------------------------------------% %
%%
clc; clear;
echo on
%%%%%%%%%%%%%%%%%% Start Code Snippet %%%%%%%%%%%%%%%%%%

%%%%%   Step 0: Initialize PDE structure and spatial variables        %%%%%
tau = 1;  A = 1;  A1 = 0.5;  B = -1;  a = 1;  a2 = 0.5;
pvar s1 s2
PDE = pde_struct();


%%%%%   Step 1: Declare the state variables and input                 %%%%%
PDE.x{1}.vars = [];                                         % ODE state variable x1(t)
PDE.x{2}.vars = s1;       PDE.x{2}.dom = [0,1];             % PDE state variable x2(t,s1) with s1 in [0,1]
PDE.x{3}.vars = s2;       PDE.x{3}.dom = [1,1+tau];         % PDE state variable x3(t,s2) with s2 in [1,1+tau]
PDE.x{4}.vars = [s1;s2];  PDE.x{4}.dom = [0,1; 1,1+tau];    % PDE state variable x4(t,s1,s2)
PDE.u{1}.vars = [];                                         % Input u(t)


%%%%%   Step 2: Declare the ODE and PDE equations                     %%%%%
% ODE: x1_{t} = A*x1(t)
PDE.x{1}.term{1}.C = A;     PDE.x{1}.term{1}.x = 1;

% ODE: x1_{t} = ... + A1*x3(t,s2=1)
PDE.x{1}.term{2}.C = A1;    PDE.x{1}.term{2}.x = 3;     
PDE.x{1}.term{2}.loc = 1;       % Evaluate state at s2=1

% ODE: x1_{t} = ... + B*x2(t,s1=0);
PDE.x{1}.term{3}.C = B;     PDE.x{1}.term{3}.x = 2;     
PDE.x{1}.term{3}.loc = 0;       % Evaluate state at s1=0


% PDE: x2_{t} = x2_{s1s1} + a*x2
PDE.x{2}.term{1}.C = [1, a];    PDE.x{2}.term{1}.x = [2; 2];    
PDE.x{2}.term{1}.D = [2; 0];    % Take second order derivative wrt s1

% PDE: x2_{t} = ... + a2*x4(t,s1,s2=1);
PDE.x{2}.term{2}.C = a2;        PDE.x{2}.term{2}.x = 4;
PDE.x{2}.term{2}.loc = [s1,1];  % Evaluate state at s1=s1, s2=1


% PDE: x3_{t} = x3_{s2};
PDE.x{3}.term{1}.x = 3;
PDE.x{3}.term{1}.D = 1;         % Take first order derivative wrt s2


% PDE: x4_{t} = x4_{s2};
PDE.x{4}.term{1}.x = 4;
PDE.x{4}.term{1}.D = [0,1];     % Take first order derivative wrt s2


%%%%%   Step 3: Declare boundary conditions                           %%%%%
% BC1: 0 = x2_{s1}(t,s1=0);
PDE.BC{1}.term{1}.x = 2;
PDE.BC{1}.term{1}.D = 1;
PDE.BC{1}.term{1}.loc = 0;

% BC2: 0 = x2(t,s1=1) - u(t)
PDE.BC{2}.term{1}.x = 2;          PDE.BC{2}.term{2}.u = 1;
PDE.BC{2}.term{1}.loc = 1;        PDE.BC{2}.term{2}.C = -1;
PDE.BC{2}.term{1}.D = 0;

% BC3: 0 = x3(t,s2=1+tau) - x1(t);
PDE.BC{3}.term{1}.x = 3;          PDE.BC{3}.term{2}.x = 1;
PDE.BC{3}.term{1}.loc = 1+tau;    PDE.BC{3}.term{2}.C = -1;

% BC4: 0 = x4(t,s1,s2=1+tau) - x2(t,s1);
PDE.BC{4}.term{1}.x = 4;                PDE.BC{4}.term{2}.x = 2;
PDE.BC{4}.term{1}.loc = [s1,1+tau];     PDE.BC{4}.term{2}.loc = s1;          
                                        PDE.BC{4}.term{2}.C = -1; 

                                       
%%%%%   Step 4: initialize the PDE, and display                       %%%%%
PDE = initialize(PDE);

display_PDE(PDE);

%%%%%%%%%%%%%%%%%% End Code Snippet %%%%%%%%%%%%%%%%%%
echo off