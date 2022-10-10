% % This script illustrates the process of declaring a 2D PDE in PIETOOLS
% % 2022a, and subsequently testing statbility of this system. 
% %
% % As an example, we implement an ODE coupled to a heat equation with
% % delay, from (Kang, 2017) (see reference at bottom):
% %
% % ODE         X_{t} = A*X(t) + A1*X(t-tau) + B*x(t,s1=0)
% % PDE         x_{t} = x_{s1s1} + a*x + a2*x(t-tau);     s1 in [0,1]
% % With BCs    x_{s1}(t,s1=0) = 0;   
% %             x(t,s1=1) = u(t);   or x_{s1}(t,s1=1) = u(t);
% %
% %---------------------------------------------------------------------% %
%%
clc; clear;
echo on
%%%%%%%%%%%%%%%%%% Start Code Snippet %%%%%%%%%%%%%%%%%%

%%%%%   Step 0: Initialize the PDE structure and spatial variables    %%%%%
PDE = pde_struct();
pvar s1 s2


%%%%%   Step 1: Declare the state variabels, inputs and outputs       %%%%%
% ODE state variable x1(t).
PDE.x{1}.vars = [];

% PDE state variable x2(t,s1)
PDE.x{2}.vars = s1;         PDE.x{2}.dom = [0,3];       % s1 in [0,3]

% PDE state variables x3(t,s1,s2)
PDE.x{3}.vars = [s1,s2];    PDE.x{3}.dom = [0,3; -1,1]; % (s1,s2) in [0,3] x [-1,1];

% Inputs w1(t) and w2(t)
PDE.w{1}.vars = [];         PDE.w{2}.vars = [];

% Output z(t)
PDE.z{1}.vars = [];


%%%%%   Step 2: Declare the PDE equations                             %%%%%
% % Declare the ODE \dot{x1}(t) = ...
% First term \dot{x1}(t) = -x(t)
PDE.x{1}.term{1}.x = 1;
PDE.x{1}.term{1}.C = -1;

% Second term \dot{x1}(t) = ... + 2*int_{0}^{3} x2(t,s1) ds1
PDE.x{1}.term{2}.x = 2;
PDE.x{1}.term{2}.C = 2;
PDE.x{1}.term{2}.I{1} = [0,3];

% Third term \dot{x1}(t) = ... - 3*int_{-1}^{1} x3(t,s1=3,s2) ds2
PDE.x{1}.term{3}.x = 3; 
PDE.x{1}.term{3}.C = -3;
PDE.x{1}.term{3}.I{2} = [-1,1];
PDE.x{1}.term{3}.loc = [3,s2];

% Fourth term \dot{x1}(t) = ... + w1(t)
PDE.x{1}.term{4}.w = 1;


% % Declare the 1D PDE \dot{x2}(t,s1) = ...
% First term \dot{x2}(t,s1) = s1*x1(t)
PDE.x{2}.term{1}.x = 1;
PDE.x{2}.term{1}.C = s1;

% Second term \dot{x2}(t,s1) = ... + \partial_{s1} x2(t,s1)
PDE.x{2}.term{2}.x = 2;
PDE.x{2}.term{2}.D = 1;

% Third term \dot{x2}(t,s1) = ... - int_{-1}^{1} x2(t,s1,s2) ds2
PDE.x{2}.term{3}.x = 3;
PDE.x{2}.term{3}.C = -1;
PDE.x{2}.term{3}.I{2} = [-1,1];

% Fourth term \dot{x2}(t,s1) = - s1*w2(t)
PDE.x{2}.term{4}.w = 2;
PDE.x{2}.term{4}.C = -s1;


% % Declare the 2D PDE \dot{x3}(t,s1,s2) = ...
% First term \dot{x3}(t,s1,s2) = s1*s2*x1(t)
PDE.x{3}.term{1}.x = 1;  
PDE.x{3}.term{1}.C = s1*s2;

% Second term \dot{x3}(t,s1,s2) = ... + x2(t,s1)
PDE.x{3}.term{2}.x = 2;

% Third term \dot{x3}(t,s1,s2) = ... + \partial_{s2}x3(t,s1,s2)
PDE.x{3}.term{3}.x = 3;
PDE.x{3}.term{3}.D = [0,1];


%%% Third step: Declare the output equation
% First term z(t) = int_{0}^{3} int_{-1}^{1} x3(t,s1,s2) ds2 ds1
PDE.z{1}.term{1}.x = 3;
PDE.z{1}.term{1}.I{1} = [0,3];
PDE.z{1}.term{1}.I{2} = [-1,1];



%%%%%   Step 4: Declare the boundary conditions                       %%%%%
% % Declare the first BC 0 = F1(t)
% First term F1(t) = x2(t,s1=3)
PDE.BC{1}.term{1}.x = 2;
PDE.BC{1}.term{1}.loc = 3;

% Second term F1(t) = ... - x1(t)
PDE.BC{1}.term{2}.x = 1;
PDE.BC{1}.term{2}.C = -1;


% % Declare the second BC 0 = F2(t)
% First term F2(t) = x3(t,s1=3,s2=-1)
PDE.BC{2}.term{1}.x = 3;
PDE.BC{2}.term{1}.loc = [3,-1];

% Second term F2(t) = ... - x1(t)
PDE.BC{2}.term{2}.x = 1;
PDE.BC{2}.term{2}.C = -1;


% % Declare the third BC 0 = F3(t,s1)
% First term F3(t,s1) = \partial_{s1} x3(t,s1,s2=-1)
PDE.BC{3}.term{1}.x = 3; 
PDE.BC{3}.term{1}.loc = [s1,-1];
PDE.BC{3}.term{1}.D = [1,0];

% Second term F3(t,s1) = ... - x2(t,s1)
PDE.BC{3}.term{2}.x = 2;
PDE.BC{3}.term{2}.loc = s1;
PDE.BC{3}.term{2}.C = -1;


%%%%%   Step 5: initialize the PDE, and display                       %%%%%
PDE = initialize(PDE);

display_PDE(PDE);

%%%%%%%%%%%%%%%%%% End Code Snippet %%%%%%%%%%%%%%%%%%
echo off