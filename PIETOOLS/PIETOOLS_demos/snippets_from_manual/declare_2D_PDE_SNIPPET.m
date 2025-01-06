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

%%% Zeroth step: Define temporal variable t and spatial variable s in [a,b]
pvar s1 s2 t


%%%%%   Step 1: Declare the state variabels, inputs and outputs       %%%%%
% PDE state variable x(t,s1,s2)
x = pde_var('state',1,[s1,s2],[0,3;-1 1]); % (s1,s2) in [0,3] x [-1,1];

% Exogenous input w(t)
w = pde_var('in');    

% Controlled inputs u1(t) and u2(t,s2)
u1 = pde_var('control');
u2 = pde_var('control',1,s2,[-1,1]);

% Regulated output z(t)
z = pde_var('out'); 

% Observed outputs y1(t,s1) and y2(t,s2)
y1 = pde_var('sense',1,s1,[0,3]);
y2 = pde_var('sense',1,s2,[-1,1]);

%%%%%   Step 2: Declare the PDE equations                             %%%%%
% % Declare the PDE  x_{t}(t,s1,s2) = x_{s1s2}(t,s1,s2) + (3-s1)(1-s2)(s2+1)*w(t)
PDE_dyn=diff(x,t,1)==diff(diff(x,s1,1),s2,1)+(3-s1)*(1-s2)*(s2+1)*w;

%%%%%   Step 3: Declare the output equations                          %%%%%
% Regulated output z(t) = 10 * int_{0}^{3} int_{-1}^{1} x(t,s1,s2) ds2 ds1
PDE_z=z==10*int(int(x,s2,[-1,1]),s1,[0,3]);

% Observed output y1(t) = int_{-1}^{1} x(t,s2) ds2 + s1*w(t)
PDE_y1=y1==int(x,s2,[-1,1])+s1*w;

% Observed output y2(t) = x(t,0,s2)
PDE_y2=y2==subs(x,s1,0);

%%%%%   Step 4: Declare the boundary conditions                       %%%%%
% % Declare the first BC 0 = x(t,0,-1)
BC1=subs(subs(x,s1,0),s2,-1)==0;

% % Declare the second BC 0 = x_{s1}(t,s1,-1) - u1(t-0.5)
BC2=subs(x,s2,-1)==subs(u1,t,t-0.5);


% % Declare the third BC 0 = x_{s2}(t,0,s2) - u2(t,s2)
BC3=subs(x,s1,0)==u2;

%%% Fourth step: create the PDE system by combining the equations
PDE =[PDE_dyn;PDE_z;PDE_y1;PDE_y2;BC1;BC2;BC3];

%%%%%   Step 5: initialize the PDE, and display                       %%%%%
PDE = initialize(PDE);
display_PDE(PDE);

%%%%%%%%%%%%%%%%%% End Code Snippet %%%%%%%%%%%%%%%%%%
echo off