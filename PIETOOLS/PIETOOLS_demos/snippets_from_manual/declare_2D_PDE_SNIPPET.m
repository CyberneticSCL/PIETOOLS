% % Code Snippet illustrating how to declare 2D PDE in PIETOOLS
% % See Section 8.2.3 in Manual for a description.
%
% % This script illustrates the process of declaring a 2D PDE in PIETOOLS
% % 2022. 
% %
% % System:
% % PDE:        x_{t}(t,s1,s2) = x_{s1s2}(t,s1,s2) + (3-s1)(1-s2)(s2+1)*w(t)
% % Outputs:    z(t) = 10 int_{0}^{3} int_{-1}^{1} x(t,s1,s2) ds2 ds1
% %             y1(t,s1) = int_{-1}^{1} x(t,s2) ds2 + s1*w(t)
% %             y2(t,s2) = x(t,0,s2)
% % With BCs:   0 = x(t,0,-1)
% %             0 = x_{s1}(t,s1,-1) - u1(t-0.5)
% %             0 = x_{s2}(t,0,s2) - u2(t,s2)
% %
% % For (s1,s2) in [0,3]x[-1,1]
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
% PDE state variable x(t,s1,s2)
PDE.x{1}.vars = [s1,s2];    PDE.x{1}.dom = [0,3; -1,1]; % (s1,s2) in [0,3] x [-1,1];

% Exogenous input w(t)
PDE.w{1}.vars = [];

% Controlled inputs u1(t) and u2(t,s2)
PDE.u{1}.vars = [];
PDE.u{2}.vars = s2;     PDE.u{2}.dom = [-1,1];

% Regulated output z(t)
PDE.z{1}.vars = [];

% Observed outputs y1(t,s1) and y2(t,s2)
PDE.y{1}.vars = s1;     PDE.y{1}.dom = [0,3];
PDE.y{2}.vars = s2;     PDE.y{2}.dom = [-1,1]; 


%%%%%   Step 2: Declare the PDE equations                             %%%%%
% % Declare the PDE \dot{x}(t,s1,s2) = ...
% First term \dot{x}(t,s1,s2) = x_{s1s2}(t,s1,s2)
PDE.x{1}.term{1}.x = 1;
PDE.x{1}.term{1}.D = [1,1];     % Take first order derivative wrt s1 and s2

% Second term \dot{x}(t,s1,s2) = ... + (3-s1)(1-s2)(s2+1)*w(t)
PDE.x{1}.term{2}.w = 1;
PDE.x{1}.term{2}.C = (3-s1)*(1-s2)*(s2+1);  % Scale w with (3-s1)(1-s2)(s2+1)


%%%%%   Step 3: Declare the output equations                          %%%%%
% Regulated output z(t) = 10 * int_{0}^{3} int_{-1}^{1} x(t,s1,s2) ds2 ds1
PDE.z{1}.term{1}.x = 1;
PDE.z{1}.term{1}.C = 10;
PDE.z{1}.term{1}.I{1} = [0,3];
PDE.z{1}.term{1}.I{2} = [-1,1];

% Observed output y1(t) = int_{-1}^{1} x(t,s2) ds2 + s1*w(t)
PDE.y{1}.term{1}.x = 1;             PDE.y{1}.term{2}.w = 1;
PDE.y{1}.term{1}.I{2} = [-1,1];     PDE.y{1}.term{2}.C = s1;

% Observed output y2(t) = x(t,0,s2)
PDE.y{2}.term{1}.x = 1;
PDE.y{2}.term{1}.loc = [0,s2];


%%%%%   Step 4: Declare the boundary conditions                       %%%%%
% % Declare the first BC 0 = x(t,0,-1)
PDE.BC{1}.term{1}.x = 1;
PDE.BC{1}.term{1}.loc = [0,-1];


% % Declare the second BC 0 = x_{s1}(t,s1,-1) - u1(t-0.5)
PDE.BC{2}.term{1}.x = 1;            PDE.BC{2}.term{2}.u = 1;
PDE.BC{2}.term{1}.D = [1,0];        PDE.BC{2}.term{2}.delay = 0.5;
PDE.BC{2}.term{1}.loc = [s1,-1];    PDE.BC{2}.term{2}.C = -1;


% % Declare the third BC 0 = x_{s2}(t,0,s2) - u2(t,s2)
PDE.BC{3}.term{1}.x = 1;            PDE.BC{3}.term{2}.u = 2;
PDE.BC{3}.term{1}.D = [0,1];       
PDE.BC{3}.term{1}.loc = [0,s2];     PDE.BC{3}.term{2}.C = -1;


%%%%%   Step 5: initialize the PDE, and display                       %%%%%
PDE = initialize(PDE);
display_PDE(PDE);

%%%%%%%%%%%%%%%%%% End Code Snippet %%%%%%%%%%%%%%%%%%
echo off