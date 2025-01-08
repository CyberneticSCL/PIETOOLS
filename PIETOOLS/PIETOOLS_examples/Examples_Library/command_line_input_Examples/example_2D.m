% % Code Snippet illustrating how to declare 2D PDE in PIETOOLS
% % See also Chapter 4 and Chapter 8 of the manual
%
% % This script illustrates the process of declaring a 2D PDE in PIETOOLS
% % 2024. 
% %
% % System:
% % ODE:        x1_{t}(t) = -x1(t) +5*x2(t,3,-1) +u1(t);
% % PDE:        x2_{t}(t,s1,s2) = x2_{s1s2}(t,s1,s2) + (3-s1)(1-s2)(s2+1)*w(t);
% % Outputs:    z(t) = [x1(t); 10*int_{0}^{3} int_{-1}^{1} x2(t,s1,s2) ds2 ds1];
% %             y1(t,s1) = int_{-1}^{1} x2(t,s2) ds2 + s1*w(t);
% %             y2(t,s2) = x2(t,0,s2);
% % With BCs:   0 = x2(t,s1,-1)
% %             0 = x2(t,0,s2) - u2(t-1,s2)
% %
% % For (s1,s2) in [0,3]x[-1,1]
% %
% %---------------------------------------------------------------------% %
%%
clc; clear; clear stateNameGenerator

%%% Zeroth step: Define temporal variable t and spatial variable s in [a,b]
pvar s1 s2 t

%%%%%   Step 1: Declare the state variabels, inputs and outputs       %%%%%
% % General format: obj = pde_var(type,size,vars,dom).
% PDE state variable x(t,s1,s2)
x1 = pde_var();
x2 = pde_var('state',[s1,s2],[0,3;-1 1]); % (s1,s2) in [0,3] x [-1,1];

% Exogenous input w(t), scalar-valued
w = pde_var('in');    

% Controlled inputs u1(t) and u2(t,s2)
u1 = pde_var('control');
u2 = pde_var('control',s2,[-1,1]);

% Regulated output z(t)
z = pde_var('out',2);       % vector-valued of size 2

% Observed outputs y1(t,s1) and y2(t,s2)
y1 = pde_var('sense',s1,[0,3]);
y2 = pde_var('sense',s2,[-1,1]);

%%%%%   Step 2: Declare the ODE and  PDE equations                    %%%%%
% % Declare the ODE  x1_{t}(t) = -x1(t)+x2(t,3,-1)+u1(t);
ODE_dyn = diff(x1,'t')==-x1 +5*subs(x2,[s1;s2],[3;-1]) +u1;

% % Declare the PDE  x_{t}(t,s1,s2) = x_{s1s2}(t,s1,s2) + (3-s1)(1-s2)(s2+1)*w(t)
PDE_dyn = diff(x2,t)==diff(diff(x2,s1,1),s2,1)+(3-s1)*(1-s2)*(s2+1)*w;

%%%%%   Step 3: Declare the output equations                          %%%%%
% Regulated output z(t) = 10 * int_{0}^{3} int_{-1}^{1} x(t,s1,s2) ds2 ds1
PDE_z = z==[x1; 10*int(x2,[s1;s2],[0,3;-1,1])];

% Observed output y1(t) = int_{-1}^{1} x(t,s2) ds2 + s1*w(t)
PDE_y1 = (y1==int(x2,s2,[-1,1]) + s1*w);

% Observed output y2(t) = x(t,0,s2)
PDE_y2 = y2==subs(x2,s1,0);

%%%%%   Step 4: Declare the boundary conditions                       %%%%%
% % Declare the first BC 0 = x2(t,s1,-1)
BC1 = subs(x2,s2,-1)==0;

% % Declare the second BC 0 = x2(t,0,s2) - u2(t-0.5,s2)
BC2 = subs(x2,s1,0)==subs(u2,t,t-0.5);

%%%%%   Step 5: initialize the PDE, and display                       %%%%%
PDE = [ODE_dyn;PDE_dyn;PDE_z;PDE_y1;PDE_y2;BC1;BC2];
PDE = initialize(PDE)



% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
% To perform e.g. analysis, the system must be converted to a PIE.
% For this, we first get rid of the temporal delay, introducing a state
% variable x3(t,s,s2)=u2(t-s*tau), modeled by a transport equation:
PDE2 = expand_delays(PDE)
% Note also that the PDE has been normalized, in that the spatial domain of
% all variables is now [-1,1]:
PDE2.dom

% Then we can convert to a PIE
PIE = convert(PDE2);
% Note that "expand_delays" is also called internally by convert, if not
% already performed manually.
% Note also that the state variables get reordered when converting to the
% PIE representation.
% Finally, note that the boundary input u2(t,s2) appears as a derivative
% d/ds2 u1(t,s2) in the PIE representation, so e.g. optimal control will be
% performed in terms of d/ds2 u1(t,s2).
