%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% examples_pde_library_PIESIM_2D.m     PIETOOLS 2024
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This file contains a library of 2D examples to use with PIESIM that includes: 
% PDE sys+++tems and coupled PDE/ODE systems. 
% Each example sets up the PDE structure of the system, 
% defines initial and boundary conditions, and, if applicable, 
% an exact PDE solution for verification (exact solution for PDE/ODE systems
% are not yet available

% Inputs:
% example - the index number of example from the library to be used (see
% below)
%
% Outputs:
% PDE - PDE structure of the problem
% uinput   - user-defined boundary inputs


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% If you modify this code, document all changes carefully and include date
% authorship, and a brief description of modifications
%
% Initial coding YP  - 5_29_2021
% YP 6/16/2022 - updated to include terms format and added more examples in
% terms format


function [PDE,uinput]=examples_pde_library_PIESIM_2D(example)
syms st sx sy real;
pvar s1 s2 
pde_struct PDE;
format long;

switch example
    
    % NOTE: uinput.ic.PDE signifies initial conditions
    %       uinput.w, uinput.u signifies inhomogeneous inputs (either through
    %       boundary conditions or distributed)
    
    % In the following examples: spatial variables are denoted as x, y
    % Temporal variable is denoted as t
    % Symbolic variabls for x is sx, for y is sy, and a symbolic variable for t is st
    % Domain boundaries are specified as x\in[a,b]x[c,d]
    
%----------------------------------------
% Parabolic equation examples
%----------------------------------------
%----------------------------------------
% 2D Heat equation u_t=visc*(u_xx+u_yy)
%----------------------------------------

%----------------------------------------
%% Example 1 - 2D Heat Equation with Dirichlet-Dirichlet boundary conditions
%----------------------------------------
    case 1

% u(x,y,t)=(sin(pi*x)+sin(pi*y))*exp(-2*visc*pi^2*t) - solution for
% u_t=visc*(uxx+uyy)

% Boundary conditions: u(a,y,t)=0; u(b,y,t)=0; u(x,c,t)=0; u(x,d,t)=0
            
   
             a=-1;b=1;
             c=-1;d=1;
             uinput.dom=[a b;c d]; 
             visc = 0.5;  

% Solving PDE in the form  x_{t}   = visc*(x_{s1s1} + x_{s2s2}) + f(s1,s2,t)

c1 = visc; c2 = visc;  ne = 1;

% % % Construct the PDE.
%%% Term-based input format
PDE.x{1}.vars = [s1;s2];   PDE.x{1}.dom = uinput.dom;

PDE.x{1}.term{1}.D = [2,0; 0,2];
PDE.x{1}.term{1}.C = [c1*eye(ne), c2*eye(ne)];

% BC1: 0 = x(s1,c)                     % BC3: 0 = x(s1,d)
PDE.BC{1}.term{1}.loc = [s1,c];      PDE.BC{3}.term{1}.loc = [s1,d];
% BC2: 0 = x(a,s2)                     % BC4: 0 = x(b,s2) 
PDE.BC{2}.term{1}.loc = [a,s2];      PDE.BC{4}.term{1}.loc = [b,s2];

% % % % 
%           Exact solution, initial conditions and inhomogeneous inputs   
             
              uinput.exact =  sin(pi*sx)*sin(pi*sy)*exp(-2*visc*pi^2*st);
             % Initial conditions for the primary states of the PDE
              uinput.ic.PDE=  subs(uinput.exact,st,0);
%----------------------------------------
%% Example 2 - 2D Heat Equation with Dirichlet-Dirichlet boundary conditions
% with non-polynomial in space forcing 
%----------------------------------------
    case 2

% u(x,y,t)=(sin(pi*x)+sin(pi*y))*t - solution for
% u_t=visc*(uxx+uyy)+f(x,y,t),
% with f(x,y,t)=sin(pi*x)*sin(pi*y)*(1+2*visc*pi^2*t);

% Boundary conditions: u(a,y,t)=0; u(b,y,t)=0; u(x,c,t)=0; u(x,d,t)=0
            
   
             a=0;b=1;
             c=0;d=1;
             uinput.dom=[a b;c d]; 
             visc = 0.5;  

% Solving PDE in the form  x_{t}   = visc*(x_{s1s1} + x_{s2s2}) + f(s1,s2,t)

c1 = visc; c2 = visc;  ne = 1;

% % % Construct the PDE.
%%% Term-based input format
PDE.x{1}.vars = [s1;s2];   PDE.x{1}.dom = uinput.dom;
PDE.w{1}.vars = [s1;s2];   PDE.w{1}.dom = uinput.dom;

PDE.x{1}.term{1}.D = [2,0; 0,2];
PDE.x{1}.term{1}.C = [c1*eye(ne), c2*eye(ne)];

PDE.x{1}.term{2}.w = 1;

% BC1: 0 = x(s1,c)                     % BC3: 0 = x(s1,d)
PDE.BC{1}.term{1}.loc = [s1,c];      PDE.BC{3}.term{1}.loc = [s1,d];
% BC2: 0 = x(a,s2)                     % BC4: 0 = x(b,s2) 
PDE.BC{2}.term{1}.loc = [a,s2];      PDE.BC{4}.term{1}.loc = [b,s2];

% % % % 
%           Exact solution, initial conditions and inhomogeneous inputs   
             
              uinput.exact =  sin(pi*sx)*sin(pi*sy)*st;
             % Initial conditions for the primary states of the PDE
              uinput.ic.PDE=  subs(uinput.exact,st,0);

 % When forcing term is non-polynomial in space, both spatial and temporal
 % content of disturbance must be added through uinput.w construct

             uinput.w(1)= sin(pi*sx)*sin(pi*sy)*(1+2*visc*pi^2*st);

%----------------------------------------
%% Example 3 - 2D Heat Equation with Dirichlet-Dirichlet boundary conditions
% with polynomial in space forcing 
%----------------------------------------
    case 3

% u(x,y,t)=100*(x-a)(x-b)(y-c)(y-d)*cos(t) - solution for
% u_t=nu*(uxx+uyy)+f(x,y,t), with
% f(x,y,t)=-200*visc*((x-a)(x-b)+(y-c)(y-d))*cos(t)-100*(x-a)(x-b)(y-c)(y-d)*sin(t)

% Boundary conditions: u(a,y,t)=0; u(b,y,t)=0; u(x,c,t)=0; u(x,d,t)=0
            
   
             a=0;b=1;
             c=0;d=1;
             uinput.dom=[a b;c d]; 
             visc = 1;  

% Solving PDE in the form  x_{t}   = visc*(x_{s1s1} + x_{s2s2}) + f(s1,s2,t)

c1 = visc; c2 = visc;  ne = 1;

% % % Construct the PDE.
%%% Term-based input format
PDE.x{1}.vars = [s1;s2];   PDE.x{1}.dom = uinput.dom;

PDE.x{1}.term{1}.D = [2,0; 0,2];
PDE.x{1}.term{1}.C = [c1*eye(ne), c2*eye(ne)];

% BC1: 0 = x(s1,c)                     % BC3: 0 = x(s1,d)
PDE.BC{1}.term{1}.loc = [s1,c];      PDE.BC{3}.term{1}.loc = [s1,d];
% BC2: 0 = x(a,s2)                     % BC4: 0 = x(b,s2) 
PDE.BC{2}.term{1}.loc = [a,s2];      PDE.BC{4}.term{1}.loc = [b,s2];

%  Specify forcing term

 % When forcing term is polynomial in space, spatial 
 % content of disturbance must be added through PDE.w construct while
 % temporal content through uinput.w construct (disturbance component index must
 % match)


PDE.w{2}.vars=[];

PDE.x{1}.term{2}.w=1;
PDE.x{1}.term{2}.C=200*visc*((s1-a)*(s1-b)+(s2-c)*(s2-d));

PDE.x{1}.term{3}.w=2;
PDE.x{1}.term{3}.C=100*(s1-a)*(s1-b)*(s2-c)*(s2-d);

% % % % 
%           Exact solution, initial conditions and inhomogeneous inputs  

               uinput.exact =  100*(sx-a)*(sx-b)*(sy-c)*(sy-d)*cos(st);
%              % Initial conditions for the primary states of the PDE
               uinput.ic.PDE= subs(uinput.exact,st,0);

       % Enter temporal content of disturbance here
               uinput.w(1)= -cos(st);
               uinput.w(2)= -sin(st);
%              

%----------------------------------------
%% Example 4 - 2D Heat Equation with Dirichlet-Neumann boundary conditions
%----------------------------------------

    case 4

% u(x,y,t)=(sin(pi*x)+sin(pi*y))*exp(-v2*isc*pi^2*t) - solution for u_t=visc*(uxx+uyy)

% Boundary conditions: u(a,y,t)=0; u_x(b,y,t)=0; u_y(x,c,t)=0; u(x,d,t)=0
            
   
             a=0;b=1.5;
             c=-0.5;d=1;
             uinput.dom=[a b;c d]; 
             visc = 0.5;  

% Solving PDE in the form  x_{t}   = visc*(x_{s1s1} + x_{s2s2})

c1 = visc; c2 = visc;  ne = 1;

% % % Construct the PDE.
%%% Term-based input format
PDE.x{1}.vars = [s1;s2];   PDE.x{1}.dom = uinput.dom;

PDE.x{1}.term{1}.D = [2,0; 0,2];
PDE.x{1}.term{1}.C = [c1*eye(ne), c2*eye(ne)];

% BC1: 0 = x_s2(s1,c)                   % BC3: 0 = x(s1,d)
PDE.BC{1}.term{1}.D = [0,1];         PDE.BC{3}.term{1}.loc = [s1,d];
PDE.BC{1}.term{1}.loc = [s1,c];      
% BC2: 0 = x(a,s2)                      % BC4: 0 = x_s1(b,s2)
PDE.BC{2}.term{1}.loc = [a,s2];      PDE.BC{4}.term{1}.D = [1,0];        
                                     PDE.BC{4}.term{1}.loc = [b,s2];

% % % % 
%           Exact solution, initial conditions and inhomogeneous inputs  
             
              uinput.exact =  sin(pi*sx)*sin(pi*sy)*exp(-2*visc*pi^2*st);

             % Initial conditions for the primary states of the PDE
              uinput.ic.PDE=  subs(uinput.exact,st,0);

%----------------------------------------
%% Example 5 - 2D KISS model
%----------------------------------------

    case 5

% u(x,y,t)=ampl*(x-a)(x-b)(y-c)(y-d) - solution for
% u_t=lam*u+c1*uxx+c2*uyy+f(x,y), with
% f(x,y)=-2*ampl*(c2*(sx-a)*(sx-b)+c1*(sy-c)*(sy-d))-lam*ampl*(x-a)(x-b)(y-c)(y-d) 


% Boundary conditions: u(a,y,t)=0; u(b,y,t)=0; u(x,c,t)=0; u(x,d,t)=0
            
   
             a=0;b=1;
             c=0;d=1;
             uinput.dom=[a b;c d]; 

% Solving PDE in the form  x_{t}   = lam*x + c1*x_{s1s1} + c2*x_{s2s2} + f(s1,s2,t)

c1 = s1; c2 = s2;  ne = 1; lam=19;

ampl=100; 

% % % Construct the PDE.
%%% Term-based input format
PDE.x{1}.vars = [s1;s2];   PDE.x{1}.dom = uinput.dom;


PDE.x{1}.term{1}.D = [0,0; 2,0; 0,2];
PDE.x{1}.term{1}.C = [lam*eye(ne), c1*eye(ne), c2*eye(ne)];

% BC1: 0 = x(s1,c)                     % BC3: 0 = x(s1,d)
PDE.BC{1}.term{1}.loc = [s1,c];      PDE.BC{3}.term{1}.loc = [s1,d];
% BC2: 0 = x(a,s2)                     % BC4: 0 = x(b,s2) 
PDE.BC{2}.term{1}.loc = [a,s2];      PDE.BC{4}.term{1}.loc = [b,s2];


%  Specify forcing term

 % When forcing term is polynomial in space, spatial 
 % content of disturbance must be added through PDE.w construct while
 % temporal content through uinput.w construct (disturbance component index must
 % match)
PDE.w{1}.vars=[];

PDE.x{1}.term{2}.w=1;
PDE.x{1}.term{2}.C=-2*ampl*(c2*(s1-a)*(s1-b)+c1*(s2-c)*(s2-d))-lam*ampl*(s1-a)*(s1-b)*(s2-c)*(s2-d);


% % % % 
%           Exact solution, initial conditions and inhomogeneous inputs  
              
               uinput.exact =  ampl*(sx-a)*(sx-b)*(sy-c)*(sy-d);
%              % Initial conditions for the primary states of the PDE
               uinput.ic.PDE= uinput.exact;

       % Enter temporal content of disturbance here
               uinput.w(1)=1;
%               uinput.w(1) = -2*ampl*(sy*(sx-a)*(sx-b)+sx*(sy-c)*(sy-d))-lam*ampl*(sx-a)*(sx-b)*(sy-c)*(sy-d);
%  

%----------------------------------------
%% Example 6 - 2D equation with different differentiability
%----------------------------------------
    case 6

% u(x,y,t)=100*(x-a)(x-b)(y-c)*(y-d)*g(t) - solution for
% u_t=lam*u+c1*uxx+c2*uy+f(x,y,t)

% Boundary conditions: u(a,y,t)=0; u(b,y,t)=0; u(x,c,t)=0
            
             a=0;b=1;
             c=0;d=1;
             uinput.dom=[a b;c d]; 

% Solving PDE in the form  x_{t}   = lam*x + c1*x_{s1s1} + c2*x_{s2} + f(s1,s2,t)

c1 = 1; c2 = 2;  ne = 1; lam=19;

ampl=100; 

% % % Construct the PDE.
%%% Term-based input format
PDE.x{1}.vars = [s1;s2];   PDE.x{1}.dom = uinput.dom;


PDE.x{1}.term{1}.D = [0,0; 2,0; 0,1];
PDE.x{1}.term{1}.C = [lam*eye(ne), c1*eye(ne), c2*eye(ne)];

% BC1: 0 = x(s1,c)                     
PDE.BC{1}.term{1}.loc = [s1,c];      
% BC2: 0 = x(a,s2)                     % BC3: 0 = x(b,s2) 
PDE.BC{2}.term{1}.loc = [a,s2];      PDE.BC{3}.term{1}.loc = [b,s2];

            
%  Specify forcing term

 % When forcing term is polynomial in space, spatial 
 % content of disturbance must be added through PDE.w construct while
 % temporal content through uinput.w construct (disturbance component index must
 % match)

PDE.w{1}.vars=[];

PDE.x{1}.term{2}.w=1;
PDE.x{1}.term{2}.C=-ampl*(c2*(s1-a)*(s1-b)*(2*s2-d-c)+2*c1*(s2-c)*(s2-d))-lam*ampl*(s1-a)*(s1-b)*(s2-c)*(s2-d);


% % % % 
%           Exact solution, initial conditions and inhomogeneous inputs  
              
               uinput.exact =  ampl*(sx-a)*(sx-b)*(sy-c)*(sy-d);
%              % Initial conditions for the primary states of the PDE

               uinput.ic.PDE = uinput.exact;

       % Enter temporal content of disturbance here
               uinput.w(1)=1;


%----------------------------------------
%% Example 7 - 2D heat equation with Neumann-Neumann BCs
%----------------------------------------
case 7
% % Standard 2D heat equation with Neumann BCs
% %   u_t=visc*(u_xx + u_yy)                        (x,y) in [a,b]x[c,d]
% %   u_{x}(a,y,t)=0; u_{x}(b,y,t)=0; u_{y}(x,c,t)=0;  u_{y}(x,d,t)=0; 
% % Starting with initial conditions
% %   u(x,y,0) = 1 + ampl*cos(pi*(x-a)/(b-a)) -ampl*cos(3*pi*(y-c)/(d-c));
% % We get exact solution
% %   u(x,y,t) = 1 + ampl*cos(pi*(x-a)/(b-a))*exp(-visc*pi^2*t/(b-a)^2) -ampl*cos(3*pi*(y-c)/(d-c))*exp(-visc*9*pi^2*t/(d-c)^2)

% % Set the spatial domain 
a=0;   b=1;
c=0;   d=3;
uinput.dom=[a b;c d]; 

% Set the viscosity and amplitude
visc = 0.5;       ampl = 2;

% % We can't represent the system with Neumann BCs directly as a PIE. 
% % Instead, introduce U0(t) = int_{a}^{b}int_{c}^{d} u(x,y,t) dy dx,
% % U1(x,t) = int_{c}^{d} u(x,y,t)dy and U2(y,t) = int_{a}^{b} u(x,y,t)dx. 
% % Then
% %   U0_{t}(t) = int_{a}^{b}int_{c}^{d}u_{t}(x,y,t) dy dx
% %             = visc*int_{c}^{d}int_{a}^{b}u_{xx}(x,y,t) dx dy + visc*int_{a}^{b}int_{c}^{d}u_{yy}(x,y,t) dy dx
% %             = 0;
% % and
% %   U1_{t}(x,t) = int_{c}^{d}u_{t}(x,y,t) dy
% %               = visc*int_{c}^{d}u_{xx}(x,y,t)dy +visc*int_{c}^{d}u_{yy}(x,y,t)dy
% %               = visc*U1_{xx} + 0;
% %   int_{a}^{b} U1(x,t)dx = U0(t);      U1_{x}(b,t) = 0;
% % and similarly
% %   U2_{t}(y,t) = visc*U2_{yy};
% %   int_{c}^{d} U2(y,t)dy = U0(t);      U2_{y}(d,t) = 0;
% % So, we get a coupled system
% %     U0_{t}(t) = 0;
% %     U1_{t}(x,t) = visc*U1_{xx}(x,t);
% %     U2_{t}(y,t) = visc*U2_{yy}(y,t);
% %     u_{t}(x,y,t) = visc*(u_xx(x,t) + u_yy(y,t))
% % with BCs
% %     int_{a}^{b} U1(x,t)dx = U0(t);      U1_{x}(b,t) = 0;
% %     int_{c}^{d} U2(y,t)dy = U0(t);      U2_{y}(d,t) = 0;
% %     int_{c}^{d} u(x,y,t)dy = U1(x,t);   int_{a}^{b} u(x,y,t)dx = U2(y,t); 
% %     u_{x}(b,y,t) = 0;   u_{y}(x,d,t) = 0;

PDE.x{1}.vars = [];
PDE.x{2}.vars = s1;         PDE.x{2}.dom = [a,b];
PDE.x{3}.vars = s2;         PDE.x{3}.dom = [c,d];
PDE.x{4}.vars = [s1;s2];    PDE.x{4}.dom = [a,b;c,d];

% U0_{t}(t) = 0;
PDE.x{1}.term{1}.x = 1;     PDE.x{1}.term{1}.C = 0;

% U1_{t}(x,t) = visc*U1_{xx}(x,t);
PDE.x{2}.term{1}.x = 2;     PDE.x{2}.term{1}.C = visc;
PDE.x{2}.term{1}.D = 2;

% U2_{t}(y,t) = visc*U2_{yy}(y,t);
PDE.x{3}.term{1}.x = 3;     PDE.x{3}.term{1}.C = visc;
PDE.x{3}.term{1}.D = 2;

% u_{t}(x,y,t) = visc*(u_xx(x,t) + u_yy(y,t))
PDE.x{4}.term{1}.x = 4;     PDE.x{4}.term{1}.C = [visc,visc];
PDE.x{4}.term{1}.D = [2,0;0,2];

% % BCs
% int_{a}^{b} U1(x,t)dx = U0(t);        U1_{x}(b,t) = 0;
PDE.BC{1}.term{1}.x = 1;                PDE.BC{2}.term{1}.x = 2;
PDE.BC{1}.term{1}.C = -1;
PDE.BC{1}.term{2}.x = 2;                PDE.BC{2}.term{1}.loc = b;
PDE.BC{1}.term{2}.I = {[a,b]};          PDE.BC{2}.term{1}.D = 1;

% int_{c}^{d} U2(y,t)dy = U0(t);        U2_{y}(d,t) = 0;
PDE.BC{3}.term{1}.x = 1;                PDE.BC{4}.term{1}.x = 3;
PDE.BC{3}.term{1}.C = -1;
PDE.BC{3}.term{2}.x = 3;                PDE.BC{4}.term{1}.loc = d;
PDE.BC{3}.term{2}.I = {[c,d]};          PDE.BC{4}.term{1}.D = 1;

% int_{c}^{d} u(x,y,t)dy = U1(x,t);     int_{a}^{b} u(x,y,t)dx = U2(y,t); 
PDE.BC{5}.term{1}.x = 2;                PDE.BC{6}.term{1}.x = 3;
PDE.BC{5}.term{1}.C = -1;               PDE.BC{6}.term{1}.C = -1;
PDE.BC{5}.term{2}.x = 4;                PDE.BC{6}.term{2}.x = 4;
PDE.BC{5}.term{2}.I = {[];[c,d]};       PDE.BC{6}.term{2}.I = {[a,b];[]};

% u_{x}(b,y,t) = 0;                     u_{y}(x,d,t) = 0;
PDE.BC{7}.term{1}.x = 4;                PDE.BC{8}.term{1}.x = 4;
PDE.BC{7}.term{1}.D = [1,0];            PDE.BC{8}.term{1}.D = [0,1];
PDE.BC{7}.term{1}.loc = [b,s2];         PDE.BC{8}.term{1}.loc = [s1,d];


% % Set the initial conditions
% %   u(x,y,0) = 1 + ampl*cos(pi*(x-a)/(b-a)) -ampl*cos(3*pi*(y-c)/(d-c));
u_ic = 1 +ampl*cos(sym(pi)*(sx-a)/(b-a)) -ampl*cos(3*sym(pi)*(sy-c)/(d-c));
uinput.ic.x(1) = int(int(u_ic,sy,c,d),sx,a,b);
uinput.ic.x(2) = int(u_ic,sy,c,d);
uinput.ic.x(3) = int(u_ic,sx,a,b);
uinput.ic.x(4) = u_ic;

% % Set the exact solution
% %   u(x,y,t) = 1 + ampl*cos(pi*(x-a)/(b-a))*exp(-visc*pi^2*t/(b-a)^2) -ampl*cos(3*pi*(y-c)/(d-c))*exp(-visc*9*pi^2*t/(d-c)^2)
u_ex = 1 +ampl*cos(sym(pi)*(sx-a)/(b-a))*exp(-visc*sym(pi)^2*st/(b-a)^2) -ampl*cos(3*sym(pi)*(sy-c)/(d-c))*exp(-visc*9*sym(pi)^2*st/(d-c)^2);
uinput.exact(1) = int(int(u_ex,sy,c,d),sx,a,b);
uinput.exact(2) = int(u_ex,sy,c,d);
uinput.exact(3) = int(u_ex,sx,a,b);
uinput.exact(4) = u_ex;



%----------------------------------------
%% Example 8 - 2D heat equation with Periodic BCs
%----------------------------------------
case 8
% % Standard 2D heat equation with Periodic BCs along x-direction
% %   u_t=visc*(u_xx + u_yy)                        (x,y) in [a,b]x[c,d]
% %   u(a,y,t) = u(b,y,t);      u_{x}(a,y,t) = u_{x}(b,y,t); 
% %   u(x,c,t) = 0;             u_{y}(x,d,t) = 0; 
% % Starting with initial conditions
% %   u(x,y,0) = ampl1*sin(1.5*pi*(y-c)/(d-c)) +ampl2*cos(6*pi*(x-0.5*(a+b))/(b-a))*sin(0.5*pi*(y-c)/(d-c)); 
% % We get exact solution
% %   u(x,y,t) = ampl1*sin(1.5*pi*(y-c)/(d-c))*exp(-visc*pi^2*t*(1.5^2/(d-c)^2)) +ampl2*cos(6*pi*(x-0.5*(a+b))/(b-a))*sin(0.5*pi*(y-c)/(d-c))*exp(-visc*pi^2*t*(36/(b-a)^2 +0.25/(d-c)^2)); 

% % Set the spatial domain 
a=-1;   b=1;
c=0;   d=1;
uinput.dom=[a b;c d]; 

% Set the viscosity and amplitude
visc = 0.1;       ampl1 = 5;   ampl2 = 10;

% % We can't represent the system with periodic BCs directly as a PIE. 
% % Instead, introduce U1(y,t) = int_{a}^{b} u(x,y,t)dx. Then
% %   U1_{t}(y,t) = int_{a}^{b}u_{t}(x,y,t) dx
% %               = visc*int_{a}^{b}u_{xx}(x,y,t)dx +visc*int_{a}^{b}u_{xx}(x,y,t)dx
% %               = 0 + visc*U1_{yy} + 0;
% %   U1(c,t) = 0;          U1_{x}(d,t) = 0;
% % So, we get a coupled system
% %     U1_{t}(x,t) = visc*U1_{yy}(y,t);
% %     u_{t}(x,y,t) = visc*(u_xx(x,t) + u_yy(y,t))
% % with BCs
% %     U1(c,t) = 0;        U1_{x}(d,t) = 0;
% %     int_{a}^{b} u(x,y,t)dy = U1(x,t);   u(b,y,t) = u(a,y,t) 
% %     u(x,c,t) = 0;       u_{y}(x,d,t) = 0;

PDE.x{1}.vars = s2;         PDE.x{1}.dom = [c,d];
PDE.x{2}.vars = [s1;s2];    PDE.x{2}.dom = [a,b;c,d];
% PDE.w{1}.vars = [];

% U1_{t}(y,t) = visc*U1_{yy}(y,t);
PDE.x{1}.term{1}.x = 1;     PDE.x{1}.term{1}.C = visc;
PDE.x{1}.term{1}.D = 2;

% u_{t}(x,y,t) = visc*(u_xx(x,t) + u_yy(y,t))
PDE.x{2}.term{1}.x = 2;     PDE.x{2}.term{1}.C = [visc,visc];
PDE.x{2}.term{1}.D = [2,0;0,2];

% % BCs
% U1(c,t) = (b-a)*w(t);                 U1_{y}(d,t) = 0;
PDE.BC{1}.term{1}.x = 1;                PDE.BC{2}.term{1}.x = 1;
PDE.BC{1}.term{1}.loc = c;              PDE.BC{2}.term{1}.loc = d;
                                        PDE.BC{2}.term{1}.D = 1;
% PDE.BC{1}.term{2}.w = 1;
% PDE.BC{1}.term{2}.C = -(b-a);

% int_{a}^{b} u(x,y,t)dx = U1(x,t);     u(a,y,t) = u(b,y,t); 
PDE.BC{3}.term{1}.x = 1;                PDE.BC{4}.term{1}.x = 2;
PDE.BC{3}.term{1}.C = -1;               PDE.BC{4}.term{1}.C = -1;
                                        PDE.BC{4}.term{1}.loc = [a,s2];
PDE.BC{3}.term{2}.x = 2;                PDE.BC{4}.term{2}.x = 2;
PDE.BC{3}.term{2}.I = {[a,b];[]};       PDE.BC{4}.term{2}.loc = [b,s2];

% u(x,c,t) = 0;                         u_{y}(x,d,t) = 0;
PDE.BC{5}.term{1}.x = 2;                PDE.BC{6}.term{1}.x = 2;
PDE.BC{5}.term{1}.loc = [s1,c];         PDE.BC{6}.term{1}.loc = [s1,d];
                                        PDE.BC{6}.term{1}.D = [0,1];
% PDE.BC{5}.term{2}.w = 1;
% PDE.BC{5}.term{2}.C = -1;

% % Set the initial conditions
% %   u(x,y,0) = ampl1*sin(1.5*pi*(y-c)/(d-c)) +ampl2*cos(6*pi*(x-0.5*(a+b))/(b-a))*sin(0.5*pi*(y-c)/(d-c)); 
u_ic = ampl1*sin(1.5*sym(pi)*(sy-c)/(d-c)) +ampl2*cos(6*sym(pi)*(sx-0.5*(a+b))/(b-a))*sin(0.5*sym(pi)*(sy-c)/(d-c));
uinput.ic.x(1) = int(u_ic,sx,a,b);
uinput.ic.x(2) = u_ic;

% % Set the exact solution
% %   u(x,y,t) = ampl1*sin(1.5*pi*(y-c)/(d-c))*exp(-visc*pi^2*t*(1.5^2/(d-c)^2)) +ampl2*cos(6*pi*(x-0.5*(a+b))/(b-a))*sin(0.5*pi*(y-c)/(d-c))*exp(-visc*pi^2*t*(36/(b-a)^2 +0.25/(d-c)^2)); 
u_ex = ampl1*sin(1.5*sym(pi)*(sy-c)/(d-c))*exp(-visc*sym(pi)^2*st*(1.5^2/(d-c)^2)) +ampl2*cos(6*sym(pi)*(sx-0.5*(a+b))/(b-a))*sin(0.5*sym(pi)*(sy-c)/(d-c))*exp(-visc*sym(pi)^2*st*(36/(b-a)^2 +0.25/(d-c)^2)); 
uinput.exact(1) = int(u_ex,sx,a,b);
uinput.exact(2) = u_ex;


% %----------------------------------------
% %% Example 8 - 2D heat equation with Periodic BCs in both directions
% %----------------------------------------
% case 8
% % % Standard 2D heat equation with Periodic BCs
% % %   u_t=visc*(u_xx + u_yy)                        (x,y) in [a,b]x[c,d]
% % %   u(a,y,t) = u(b,y,t);    u_{x}(a,y,t) = u_{x}(b,y,t); 
% % %   u(x,c,t) = u(d,y,t);    u_{y}(x,c,t) = u_{y}(x,d,t); 
% % % Starting with initial conditions
% % %   u(x,y,0) = -1 + ampl*cos(2*pi*(x-0.5*(a+b))/(b-a))*sin(6*pi*(y-0.5*(c+d))/(d-c));
% % % We get exact solution
% % %   u(x,y,t) = -1 + ampl*cos(2*pi*(x-0.5*(a+b))/(b-a))*sin(6*pi*(y-0.5*(c+d))/(d-c))*exp(-visc*4*pi^2*t/(b-a)^2 -visc*36*pi^2*t/(d-c)^2)
% 
% % % Set the spatial domain 
% a=-1;   b=1;
% c=-3;   d=3;
% uinput.dom=[a b;c d]; 
% 
% % Set the viscosity and amplitude
% visc = 0.5;       ampl = 2;
% 
% % % We can't represent the system with periodic BCs directly as a PIE. 
% % % Instead, introduce U0(t) = int_{a}^{b}int_{c}^{d} u(x,y,t) dy dx,
% % % U1(x,t) = int_{c}^{d} u(x,y,t)dy and U2(y,t) = int_{a}^{b} u(x,y,t)dx. 
% % % Then
% % %   U0_{t}(t) = int_{a}^{b}int_{c}^{d}u_{t}(x,y,t) dy dx
% % %             = visc*int_{c}^{d}int_{a}^{b}u_{xx}(x,y,t) dx dy + visc*int_{a}^{b}int_{c}^{d}u_{yy}(x,y,t) dy dx
% % %             = visc*int_{c}^{d}(u_{x}(b,y,t)-u_{x}(a,y,t)) dy + visc*int_{a}^{b}(u_{y}(x,d,t)-u_{y}(x,c,t)) dx
% % %             = 0;
% % % and
% % %   U1_{t}(x,t) = int_{c}^{d}u_{t}(x,y,t) dy
% % %               = visc*int_{c}^{d}u_{xx}(x,y,t)dy +visc*int_{c}^{d}u_{yy}(x,y,t)dy
% % %               = visc*U1_{xx} + 0;
% % %   int_{a}^{b} U1(x,t)dx = U0(t);      U1(b,t) = U1(a,t);
% % % and similarly
% % %   U2_{t}(y,t) = visc*U2_{yy};
% % %   int_{c}^{d} U2(y,t)dy = U0(t);      U2(d,t) = U2(c,t);
% % % So, we get a coupled system
% % %     U0_{t}(t) = 0;
% % %     U1_{t}(x,t) = visc*U1_{xx}(x,t);
% % %     U2_{t}(y,t) = visc*U2_{yy}(y,t);
% % %     u_{t}(x,y,t) = visc*(u_xx(x,t) + u_yy(y,t))
% % % with BCs
% % %     int_{a}^{b} U1(x,t)dx = U0(t);      U1(b,t) = U1(a,t);
% % %     int_{c}^{d} U2(y,t)dy = U0(t);      U2(d,t) = U2(c,t);
% % %     int_{c}^{d} u(x,y,t)dy = U1(x,t);   int_{a}^{b} u(x,y,t)dx = U2(y,t); 
% % %     u(b,y,t) = u(a,y,t);        u(x,d,t) = u(x,c,t);
% 
% PDE.x{1}.vars = [];
% PDE.x{2}.vars = s1;         PDE.x{2}.dom = [a,b];
% PDE.x{3}.vars = s2;         PDE.x{3}.dom = [c,d];
% PDE.x{4}.vars = [s1;s2];    PDE.x{4}.dom = [a,b;c,d];
% 
% % U0_{t}(t) = 0;
% PDE.x{1}.term{1}.x = 1;     PDE.x{1}.term{1}.C = 0;
% 
% % U1_{t}(x,t) = visc*U1_{xx}(x,t);
% PDE.x{2}.term{1}.x = 2;     PDE.x{2}.term{1}.C = visc;
% PDE.x{2}.term{1}.D = 2;
% 
% % U2_{t}(y,t) = visc*U2_{yy}(y,t);
% PDE.x{3}.term{1}.x = 3;     PDE.x{3}.term{1}.C = visc;
% PDE.x{3}.term{1}.D = 2;
% 
% % u_{t}(x,y,t) = visc*(u_xx(x,t) + u_yy(y,t))
% PDE.x{4}.term{1}.x = 4;     PDE.x{4}.term{1}.C = [visc,visc];
% PDE.x{4}.term{1}.D = [2,0;0,2];
% 
% % % BCs
% % int_{a}^{b} U1(x,t)dx = U0(t);        U1(b,t) = U1(a,t);
% PDE.BC{1}.term{1}.x = 1;                PDE.BC{2}.term{1}.x = 2;
% PDE.BC{1}.term{1}.C = -1;               PDE.BC{2}.term{1}.C = -1;
% PDE.BC{1}.term{2}.x = 2;                PDE.BC{2}.term{1}.loc = b;
% PDE.BC{1}.term{2}.I = {[a,b]};          PDE.BC{2}.term{2}.x = 2;
%                                         PDE.BC{2}.term{2}.loc = a;
% 
% % int_{c}^{d} U2(y,t)dy = U0(t);        U2(d,t) = U2(c,t);
% PDE.BC{3}.term{1}.x = 1;                PDE.BC{4}.term{1}.x = 3;
% PDE.BC{3}.term{1}.C = -1;               PDE.BC{4}.term{1}.C = -1;
% PDE.BC{3}.term{2}.x = 3;                PDE.BC{4}.term{1}.loc = d;
% PDE.BC{3}.term{2}.I = {[c,d]};          PDE.BC{4}.term{2}.x = 3;
%                                         PDE.BC{4}.term{2}.loc = c;
% 
% % int_{c}^{d} u(x,y,t)dy = U1(x,t);     int_{a}^{b} u(x,y,t)dx = U2(y,t); 
% PDE.BC{5}.term{1}.x = 2;                PDE.BC{6}.term{1}.x = 3;
% PDE.BC{5}.term{1}.C = -1;               PDE.BC{6}.term{1}.C = -1;
% PDE.BC{5}.term{2}.x = 4;                PDE.BC{6}.term{2}.x = 4;
% PDE.BC{5}.term{2}.I = {[];[c,d]};       PDE.BC{6}.term{2}.I = {[a,b];[]};
% 
% % u(b,y,t) = u(a,y,t);                  u(x,d,t) = u(x,c,t);
% PDE.BC{7}.term{1}.x = 4;                PDE.BC{8}.term{1}.x = 4;
% PDE.BC{7}.term{1}.C = -1;               PDE.BC{8}.term{1}.C = -1;
% PDE.BC{7}.term{1}.loc = [b,s2];         PDE.BC{8}.term{1}.loc = [s1,d];
% PDE.BC{7}.term{2}.x = 4;                PDE.BC{8}.term{2}.x = 4;
% PDE.BC{7}.term{2}.loc = [a,s2];         PDE.BC{8}.term{2}.loc = [s1,c];
% 
% 
% % % Set the initial conditions
% % %   u(x,y,0) = -1 + ampl*cos(2*pi*(x-0.5*(a+b))/(b-a))*sin(6*pi*(y-0.5*(c+d))/(d-c));
% u_ic = -1 + ampl*cos(2*sym(pi)*(sx-0.5*(a+b))/(b-a))*sin(6*sym(pi)*(sy-0.5*(c+d))/(d-c));
% uinput.ic.x(1) = int(int(u_ic,sy,c,d),sx,a,b);
% uinput.ic.x(2) = int(u_ic,sy,c,d);
% uinput.ic.x(3) = int(u_ic,sx,a,b);
% uinput.ic.x(4) = u_ic;
% 
% % % Set the exact solution
% % %   u(x,y,t) = -1 + ampl*cos(2*pi*(x-0.5*(a+b))/(b-a))*sin(6*pi*(y-0.5*(c+d))/(d-c))*exp(-visc*4*pi^2*t/(b-a)^2 -visc*36*pi^2*t/(d-c)^2)
% u_ex = -1 + ampl*cos(2*sym(pi)*(sx-0.5*(a+b))/(b-a))*sin(6*sym(pi)*(sy-0.5*(c+d))/(d-c))*exp(-visc*st*sym(pi)^2*(4/(b-a)^2 -36/(d-c)^2));
% uinput.exact(1) = int(int(u_ex,sy,c,d),sx,a,b);
% uinput.exact(2) = int(u_ex,sy,c,d);
% uinput.exact(3) = int(u_ex,sx,a,b);
% uinput.exact(4) = u_ex;



%----------------------------------------
%% Example 9 - 2D heat equation coupled to 1D heat equation at boundary
%----------------------------------------
case 9
% % Coupled ODE - 1D PDE - 2D PDE system
% %   X_{t}(t)      = -visc*X(t);
% %   u1_{t}(x,t)   = visc*(b-a)^2/pi^2*U1_{xx}
% %   u2_{t}(y,t)   = visc*(d-c)^2/(25*pi^2)*U2_{yy}
% %   u3_{t}(x,y,t) = 0.5*visc/pi^2*((b-a)^2*u_{xx}+(d-c)^2/25 *u_{yy})    (x,y) in [a,b]x[c,d]
% % with BCs
% %   u1(a,t) = X(t)    u1_{x}(b,t) = 0;
% %   u2(c,t) = X(t)    u2_{y}(d,t) = 0;
% %   u3(x,c,t) = u1(x,t);      u3_{y}(x,d,t) = 0;
% %   u3(a,y,t) = u2(y,t);      u3_{x}(b,y,t) = 0;
% % Starting with initial conditions
% %   X(0)      = ampl;
% %   u1(x,0)   = ampl*cos(pi*(x-a)/(b-a));
% %   u2(y,0)   = ampl*cos(5*pi*(y-c)/(d-c));
% %   u3(x,y,0) = ampl*cos(pi*(x-a)/(b-a))*cos(5*pi*(y-c)/(d-c));
% % We get exact solution
% %   X(t)      = ampl*exp(-visc*t);
% %   u1(x,t)   = ampl*cos(pi*(x-a)/(b-a))*exp(-visc*t);
% %   u2(y,t)   = ampl*cos(5*pi*(y-c)/(d-c))*exp(-visc*t);
% %   u3(x,y,t) = ampl*cos(pi*(x-a)/(b-a))*cos(5*pi*(y-c)/(d-c))*exp(-visc*t);
% 

% % Set the spatial domain 
a=0;   b=1;
c=0;   d=5;
uinput.dom=[a b;c d]; 

% Set the viscosity and amplitude
visc = 1;       ampl = 10;

% % Declare the PDE
PDE.x{1}.vars = [];
PDE.x{2}.vars = s1;         PDE.x{2}.dom = [a,b];
PDE.x{3}.vars = s2;         PDE.x{3}.dom = [c,d];
PDE.x{4}.vars = [s1;s2];    PDE.x{4}.dom = [a,b;c,d];

% X_{t}(t) = -visc*X(t);
PDE.x{1}.term{1}.x = 1;     PDE.x{1}.term{1}.C = -visc;

% u1_{t}(x,t)   = visc*(b-a)^2/pi^2*u1_{xx}
PDE.x{2}.term{1}.x = 2;     PDE.x{2}.term{1}.C = visc*(b-a)^2/pi^2;
PDE.x{2}.term{1}.D = 2;

% u2_{t}(y,t)   = visc*(d-c)^2/(25*pi^2)*u2_{yy}
PDE.x{3}.term{1}.x = 3;     PDE.x{3}.term{1}.C = visc*(d-c)^2/(25*pi^2);
PDE.x{3}.term{1}.D = 2;

% u3_{t}(x,y,t) = 0.5*visc/pi^2*((b-a)^2*u_{xx}+(d-c)^2/25 *u_{yy})
PDE.x{4}.term{1}.x = 4;     PDE.x{4}.term{1}.C = 0.5*visc/pi^2 *[(b-a)^2,(d-c)^2/25];
PDE.x{4}.term{1}.D = [2,0;0,2];

% % BCs
% u1(a,t) = X(t);                       u1_{x}(b,t) = 0;
PDE.BC{1}.term{1}.x = 1;                PDE.BC{2}.term{1}.x = 2;
PDE.BC{1}.term{1}.C = -1;
PDE.BC{1}.term{2}.x = 2;                PDE.BC{2}.term{1}.loc = b;
PDE.BC{1}.term{2}.loc = a;              PDE.BC{2}.term{1}.D = 1;

% u2(c,t) = X(t);                       u2_{y}(d,t) = 0;
PDE.BC{3}.term{1}.x = 1;                PDE.BC{4}.term{1}.x = 3;
PDE.BC{3}.term{1}.C = -1;
PDE.BC{3}.term{2}.x = 3;                PDE.BC{4}.term{1}.loc = d;
PDE.BC{3}.term{2}.loc = c;              PDE.BC{4}.term{1}.D = 1;

% u3(x,c,t) = u1(x,t);                  u3(a,y,t) = u2(y,t); 
PDE.BC{5}.term{1}.x = 2;                PDE.BC{6}.term{1}.x = 3;
PDE.BC{5}.term{1}.C = -1;               PDE.BC{6}.term{1}.C = -1;
PDE.BC{5}.term{2}.x = 4;                PDE.BC{6}.term{2}.x = 4;
PDE.BC{5}.term{2}.loc = [s1,c];         PDE.BC{6}.term{2}.loc = [a,s2];

% u3_{x}(b,y,t) = 0;                    u3_{y}(x,d,t) = 0;
PDE.BC{7}.term{1}.x = 4;                PDE.BC{8}.term{1}.x = 4;
PDE.BC{7}.term{1}.D = [1,0];            PDE.BC{8}.term{1}.D = [0,1];
PDE.BC{7}.term{1}.loc = [b,s2];         PDE.BC{8}.term{1}.loc = [s1,d];


% % Set the initial conditions
% %   u(x,y,0) = ampl*cos(pi*(x-a)/(b-a))*cos(5*pi*(y-c)/(d-c));
u_ic = ampl*cos(sym(pi)*(sx-a)/(b-a))*cos(5*sym(pi)*(sy-c)/(d-c));
uinput.ic.x(1) = subs(subs(u_ic,sy,c),sx,a);
uinput.ic.x(2) = subs(u_ic,sy,c);
uinput.ic.x(3) = subs(u_ic,sx,a);
uinput.ic.x(4) = u_ic;

% % Set the exact solution
% %   u(x,y,t) = ampl*cos(pi*(x-a)/(b-a))*cos(5*pi*(y-c)/(d-c))*exp(-visc*t)
u_ex = ampl*cos(sym(pi)*(sx-a)/(b-a))*cos(5*sym(pi)*(sy-c)/(d-c))*exp(-visc*st);
uinput.exact(1) = subs(subs(u_ex,sy,c),sx,a);
uinput.exact(2) = subs(u_ex,sy,c);
uinput.exact(3) = subs(u_ex,sx,a);
uinput.exact(4) = u_ex;



%% Example 10 - 2D heat equation with forcing throguh coupled 1D PDE
%----------------------------------------
case 10
% % Coupled 1D PDE - 2D PDE system
% %   u1_{t}(x,t)   = u1_{xx}(x,t) -u1(x,t) +w3(x,t)
% %   u2_{t}(y,t)   = u2_{yy}(y,t) +(2*pi/(d-c))^2*u2
% %   u3_{t}(y,t)   = 1/(2*pi/(d-c))^2 *u3_{yy}(y,t)
% %   u4_{t}(x,y,t) = u4_{xx} +u4_{yy} +(2*pi/(d-c))^2*u4 -w4(x,y,t)
% % with BCs
% %   u1(a,t) = w1(t);  u1(b,t) = w2(t);
% %   u2(c,t) = w1(t);  u2_{y}(d,t) = 0;
% %   u3(c,t) = w2(t);  u3_{y}(d,t) = 0;
% %   u4(x,c,t) = u1(x,t);      u4_{y}(x,d,t) = 0;
% %   u4(a,y,t) = u2(y,t);      u4(b,y,t) = u3(y,t);
% % Starting with initial conditions
% %   u1(x,0) = (x-b)/(a-b)*f1 +(x-a)/(b-a)*f2
% %   u2(y,0) = f1*cos(2*pi*(y-c)/(d-c))
% %   u3(y,0) = f2*cos(2*pi*(y-c)/(d-c))
% %   u4(x,y,0) = ((x-b)/(a-b)*f1 +(x-a)/(b-a)*f2)*cos(2*pi*(y-c)/(d-c))
% % Using forcing
% %   w1(t)   = f1;     w2(t) = f2*exp(-t);
% %   w3(x,t) = (x-b)/(a-b)*f1
% %   w4(x,y,t) = -f2*exp(-t)*(x-a)/(b-a)*cos(2*pi*(y-c)/(d-c))
% % We get exact solution
% %   u1(x,t) = (x-b)/(a-b)*f1 +(x-a)/(b-a)*f2*exp(-t)
% %   u2(y,t) = f1*cos(2*pi*(y-c)/(d-c))
% %   u3(y,t) = f2*exp(-t)*cos(2*pi*(y-c)/(d-c))
% %   u4(x,y,t) = ((x-b)/(a-b)*f1 +(x-a)/(b-a)*f2*exp(-t))*cos(2*pi*(y-c)/(d-c))


% % Set the spatial domain 
a=0;   b=1;
c=0;   d=3;
uinput.dom=[a b;c d]; 

% Set the forcing constants w1(t)=f1 and w2(t)=f2;
f1 = 5;    f2 = -5;

% % Declare the PDE
PDE.x{1}.vars = s1;         PDE.x{1}.dom = [a,b];
PDE.x{2}.vars = s2;         PDE.x{2}.dom = [c,d];
PDE.x{3}.vars = s2;         PDE.x{3}.dom = [c,d];
PDE.x{4}.vars = [s1;s2];    PDE.x{4}.dom = [a,b;c,d];
PDE.w{1}.vars = [];         PDE.w{2}.vars = [];
PDE.w{3}.vars = s1;         PDE.w{3}.dom = [a,b];
PDE.w{4}.vars = [s1;s2];    PDE.w{4}.dom = [a,b;c,d];

% u1_{t}(x,t)   = u1_{xx}(x,t) -u1(x,t) +w3(x,t)
PDE.x{1}.term{1}.x = 1;     PDE.x{1}.term{1}.C = [1,-1];
PDE.x{1}.term{1}.D = [2;0];
PDE.x{1}.term{2}.w = 3;

% %   u2_{t}(y,t)   = u2_{yy}(y,t) +(1.5*pi/(d-c))^2*u2
PDE.x{2}.term{1}.x = 2;
PDE.x{2}.term{1}.C = [1,(2*pi/(d-c))^2];
PDE.x{2}.term{1}.D = [2;0];

% %   u3_{t}(y,t)   = 1/(1.5*pi/(d-c))^2 *u3_{yy}(y,t)
PDE.x{3}.term{1}.x = 3;
PDE.x{3}.term{1}.C = 1/(2*pi/(d-c))^2;
PDE.x{3}.term{1}.D = 2;
% %   u4_{t}(x,y,t) = u4_{xx} +u4_{yy} +(1.5*pi/(d-c))^2*u4 +w4(x,y,t)

% u4_{t}(x,y,t) = u4_{xx} +u4_{yy} +(1.5*pi/(d-c))^2*u4 -w4(x,y,t)
PDE.x{4}.term{1}.x = 4;     
PDE.x{4}.term{1}.C = [1,1,(2*pi/(d-c))^2];
PDE.x{4}.term{1}.D = [2,0;0,2;0,0];
PDE.x{4}.term{2}.w = 4;
PDE.x{4}.term{2}.C = -1;

% % BCs
% u1(a,t) = w1(t);                      u1(b,t) = w2;
PDE.BC{1}.term{1}.x = 1;                PDE.BC{2}.term{1}.x = 1;
PDE.BC{1}.term{1}.loc = a;              PDE.BC{2}.term{1}.loc = b;
PDE.BC{1}.term{1}.C = -1;               PDE.BC{2}.term{1}.C = -1;
PDE.BC{1}.term{2}.w = 1;                PDE.BC{2}.term{2}.w = 2;

% u2(c,t) = w1(t);                      u2_{y}(d,t) = 0;
PDE.BC{3}.term{1}.x = 2;                PDE.BC{4}.term{1}.x = 2;
PDE.BC{3}.term{1}.loc = c;              PDE.BC{4}.term{1}.loc = d;
PDE.BC{3}.term{2}.w = 1;                PDE.BC{4}.term{1}.D = 1;
PDE.BC{3}.term{2}.C = -1;

% u3(c,t) = w2(t);                      u3_{y}(d,t) = 0;
PDE.BC{5}.term{1}.x = 3;                PDE.BC{6}.term{1}.x = 3;
PDE.BC{5}.term{1}.loc = c;              PDE.BC{6}.term{1}.loc = d;
PDE.BC{5}.term{2}.w = 2;                PDE.BC{6}.term{1}.D = 1;
PDE.BC{5}.term{2}.C = -1;

% u4(x,c,t) = u1(x,t);                  u4_{y}(x,d,t) = 0;
PDE.BC{7}.term{1}.x = 4;                PDE.BC{8}.term{1}.x = 4;
PDE.BC{7}.term{1}.loc = [s1,c];         PDE.BC{8}.term{1}.loc = [s1,d];
PDE.BC{7}.term{2}.x = 1;                PDE.BC{8}.term{1}.D = [0,1];
PDE.BC{7}.term{2}.C = -1; 

% u4(a,y,t) = u2(y,t);                  u4(b,y,t) = u3(y,t);
PDE.BC{9}.term{1}.x = 4;                PDE.BC{10}.term{1}.x = 4;
PDE.BC{9}.term{1}.loc = [a,s2];         PDE.BC{10}.term{1}.loc = [b,s2];
PDE.BC{9}.term{2}.x = 2;                PDE.BC{10}.term{2}.x = 3;
PDE.BC{9}.term{2}.C = -1;               PDE.BC{10}.term{2}.C = -1;


% % Set the initial conditions
% %   u1(x,0) = (x-b)/(a-b)*f1 +(x-a)/(b-a)*f2
% %   u2(y,0) = f1*cos(2*pi*(y-c)/(d-c))
% %   u3(y,0) = f2*cos(2*pi*(y-c)/(d-c))
% %   u4(x,y,0) = ((x-b)/(a-b)*f1 +(x-a)/(b-a)*f2)*cos(2*pi*(y-c)/(d-c))
u_ic = ((sx-b)/(a-b)*f1 +(sx-a)/(b-a)*f2)*cos(2*sym(pi)*(sy-c)/(d-c));
uinput.ic.x(1) = subs(u_ic,sy,c);
uinput.ic.x(2) = subs(u_ic,sx,a);
uinput.ic.x(3) = subs(u_ic,sx,b);
uinput.ic.x(4) = u_ic;

% % Set the forcing
% %   w1(t)   = f1;     w2(t) = f2*exp(-t);
% %   w3(x,t) = (x-b)/(a-b)*f1
% %   w4(x,y,t) = -f2*exp(-t)*(x-a)/(b-a)*cos(2*pi*(y-c)/(d-c))
uinput.w(1) = sym(f1);
uinput.w(2) = f2*exp(-st);
uinput.w(3) = f1*(sx-b)/(a-b);
uinput.w(4) = f2*exp(-st)*(sx-a)/(b-a)*cos(2*sym(pi)*(sy-c)/(d-c));

% % Set the exact solution
% %   u1(x,t) = (x-b)/(a-b)*f1 +(x-a)/(b-a)*f2*exp(-t)
% %   u2(y,t) = f1*cos(1.5*pi*(y-c)/(d-c))
% %   u3(y,t) = f2*exp(-t)*cos(1.5*pi*(y-c)/(d-c))
% %   u4(x,y,t) = ((x-b)/(a-b)*f1 +(x-a)/(b-a)*f2*exp(-t))*cos(2*pi*(y-c)/(d-c))
u_ex = (f1*(sx-b)/(a-b) +f2*(sx-a)/(b-a)*exp(-st))*cos(2*sym(pi)*(sy-c)/(d-c));
uinput.exact(1) = subs(u_ex,sy,c);
uinput.exact(2) = subs(u_ex,sx,a);
uinput.exact(3) = subs(u_ex,sx,b);
uinput.exact(4) = u_ex;


%----------------------------------------
%% Example 11 - 2D wave equation with Dirichlet-Dirichlet BCs
%----------------------------------------
case 11
% % Standard 2D wave equation with Dirichlet BCs
% %   u_tt=C^2*(u_xx + u_yy)                        (x,y) in [0,1]x[0,2]
% %   u(a,y,t)=0; u(b,y,t)=0; u(x,c,t)=0;  u(x,d,t)=0;   
% % Starting with initial condition
% %       u(x,y,0) = 2*sin(pi*x)*sin(pi*y) - sin(2*pi*x)*sin((pi*y)/2)
% %   u_{t}(x,y,0) = (5^(1/2)*pi*sin(2*pi*x)*sin(pi*y))/2
% % Yielding exact solution
% %   u(x,y,t) = 2*sin(pi*x)*sin(pi*y)*cos(pi*2^(1/2)*t) - sin(2*pi*x)*sin((pi*y)/2)*cos((pi*17^(1/2)*t)/2) + (sin(2*pi*x)*sin(pi*y)*sin(pi*5^(1/2)*st))/2

% % Set the spatial domain 
a=0;   b=1;
c=0;   d=2;
uinput.dom=[a b;c d]; 
% NOTE: a and c must be 0 for the exact solution to be valid.

% Set the wave velocity 
C = 1; 

% % % Construct the PDE.
% % We represent it in the form
% %     x1_{t} = x2;
% %     x2_{t} = C^2*(x1_{xx} +x1_{yy});
%%% Term-based input format
PDE.x{1}.vars = [s1;s2];   PDE.x{1}.dom = uinput.dom;   % x1 = u;
PDE.x{2}.vars = [s1;s2];   PDE.x{2}.dom = uinput.dom;   % x2 = u_{t}
PDE.x{2}.diff = [2,2];

% x1_{t} = x2;
PDE.x{1}.term{1}.x = 2;

% x2_{t} = C^2*(x1_{xx} + x1_{yy});
PDE.x{2}.term{1}.x = 1;
PDE.x{2}.term{1}.D = [2,0;0,2];
PDE.x{2}.term{1}.C = [C^2,C^2];

% % Set the boundary conditions
% BC1: 0 = x1(s1,c)                 % BC3: 0 = x1(s1,d)
PDE.BC{1}.term{1}.x = 1;            PDE.BC{3}.term{1}.x = 1;
PDE.BC{1}.term{1}.loc = [s1,c];     PDE.BC{3}.term{1}.loc = [s1,d];
% BC2: 0 = x1(a,s2)                 % BC4: 0 = x1(b,s2)
PDE.BC{2}.term{1}.x = 1;            PDE.BC{4}.term{1}.x = 1;
PDE.BC{2}.term{1}.loc = [a,s2];     PDE.BC{4}.term{1}.loc = [b,s2];

% BC5: 0 = x2(s1,c)                 % BC7: 0 = x2(s1,d)
PDE.BC{5}.term{1}.x = 2;            PDE.BC{7}.term{1}.x = 2;
PDE.BC{5}.term{1}.loc = [s1,c];     PDE.BC{7}.term{1}.loc = [s1,d];
% BC6: 0 = x2(a,s2)                 % BC8: 0 = x2(b,s2)
PDE.BC{6}.term{1}.x = 2;            PDE.BC{8}.term{1}.x = 2;
PDE.BC{6}.term{1}.loc = [a,s2];     PDE.BC{8}.term{1}.loc = [b,s2];

% % Set the initial conditions
% % We build an initial function of the form
% % 
% %     u(x,y) = sum_{m=1}^{infty} sum_{n=1}^{infty} B_{mn}*sin(m*pi/(b-a)*x)*sin(n*pi/(d-c)*y);
% %     u_{t}(x,t) = sum_{m=1}^{infty} sum_{n=1}^{infty} lam_{mn}*C_{mn}*sin(m*pi/(b-a)*x)*sin(n*pi/(d-c)*y);
% % 
% % where lam_{mn}=C*sqrt((m*pi/(b-a))^2 +(n*pi/(d-c))^2)
% Set the coefficients B_mn and C_mn, so that e.g. B_mn(m,n)=B_{mn};
B_mn = [0,2,0;-1,0,0];   C_mn = [0,0,0;0,0.5,0];
% Set the factors lam_{mn} accordingly.
lam_mn = C*sqrt(((1:size(C_mn,1))'*sym(pi)/(b-a)).^2 +((1:size(C_mn,2))*sym(pi)/(d-c)).^2);
% Set the spatial basis functions
basis_funs = sin(((1:size(B_mn,1))'*pi/(b-a))*sx)*sin(((1:size(B_mn,2))*pi/(d-c))*sy);
% Finally, set the initial function.
uinput.ic.PDE(1) = sum(sum(B_mn.*basis_funs));
uinput.ic.PDE(2) = sum(sum(C_mn.*lam_mn.*basis_funs));

% % Given this initial function, the exact solution is given by
% %     u(x,y) = sum_{m=1}^{infty} sum_{n=1}^{infty} sin(m*pi/(b-a)*x)*sin(n*pi/(d-c)*y)*(B_{mn}*cos(lam_{mn}*t)+C_{mn}*sin(lam_{mn}*t));
u_exact = sum(sum(basis_funs.*(B_mn.*cos(lam_mn*st)+C_mn.*sin(lam_mn*st))));
uinput.exact(1) = u_exact;
uinput.exact(2) = diff(u_exact,st);


%----------------------------------------
%% Example 12 - 2D wave equation with Dirichlet-Neumann BCs
%----------------------------------------
case 12
% % Standard 2D wave equation with Dirichlet and Neumann BCs
% %   u_tt=C^2*(u_xx + u_yy)                        (x,y) in [0,1]x[0,2]
% %   u_{x}(a,y,t)=0; u(b,y,t)=0; u_{y}(x,c,t)=0;  u(x,d,t)=0;   
% % Starting with initial condition
% %       u(x,y,0) = 2*sin(pi*x)*sin(pi*y) - sin(2*pi*x)*sin((pi*y)/2)
% %   u_{t}(x,y,0) = (5^(1/2)*pi*sin(2*pi*x)*sin(pi*y))/2
% % Yielding exact solution
% %   u(x,y,t) = 2*sin(pi*x)*sin(pi*y)*cos(pi*2^(1/2)*t) - sin(2*pi*x)*sin((pi*y)/2)*cos((pi*17^(1/2)*t)/2) + (sin(2*pi*x)*sin(pi*y)*sin(pi*5^(1/2)*st))/2

% % Set the spatial domain 
a=0;   b=1;
c=0;   d=0.5;
uinput.dom=[a b;c d]; 
% NOTE: a and c must be 0 for the exact solution to be valid.

% Set the wave velocity 
C = 1; 

% % % Construct the PDE.
% % We represent it in the form
% %     x1_{t} = x2;
% %     x2_{t} = C^2*(x1_{xx} +x1_{yy});
%%% Term-based input format
PDE.x{1}.vars = [s1;s2];   PDE.x{1}.dom = uinput.dom;   % x1 = u;
PDE.x{2}.vars = [s1;s2];   PDE.x{2}.dom = uinput.dom;   % x2 = u_{t}
PDE.x{2}.diff = [2,2];

% x1_{t} = x2;
PDE.x{1}.term{1}.x = 2;

% x2_{t} = C^2*(x1_{xx} + x1_{yy});
PDE.x{2}.term{1}.x = 1;
PDE.x{2}.term{1}.D = [2,0;0,2];
PDE.x{2}.term{1}.C = [C^2,C^2];

% % Set the boundary conditions
% BC1: 0 = x1(s1,c)                 % BC3: 0 = x1_{y}(s1,d)
PDE.BC{1}.term{1}.x = 1;            PDE.BC{3}.term{1}.x = 1;
PDE.BC{1}.term{1}.loc = [s1,c];     PDE.BC{3}.term{1}.loc = [s1,d];
                                    PDE.BC{3}.term{1}.D = [0,1];
% BC2: 0 = x1_{x}(a,s2)             % BC4: 0 = x1(b,s2)
PDE.BC{2}.term{1}.x = 1;            PDE.BC{4}.term{1}.x = 1;
PDE.BC{2}.term{1}.loc = [a,s2];     PDE.BC{4}.term{1}.loc = [b,s2];
PDE.BC{2}.term{1}.D = [1,0];

% BC5: 0 = x2(s1,c)                 % BC7: 0 = x2_{y}(s1,d)
PDE.BC{5}.term{1}.x = 2;            PDE.BC{7}.term{1}.x = 2;
PDE.BC{5}.term{1}.loc = [s1,c];     PDE.BC{7}.term{1}.loc = [s1,d];
                                    PDE.BC{7}.term{1}.D = [0,1];
% BC6: 0 = x2_{x}(a,s2)             % BC8: 0 = x2(b,s2)
PDE.BC{6}.term{1}.x = 2;            PDE.BC{8}.term{1}.x = 2;
PDE.BC{6}.term{1}.loc = [a,s2];     PDE.BC{8}.term{1}.loc = [b,s2];
PDE.BC{6}.term{1}.D = [1,0];

% % Set the initial conditions
% % We build an initial function of the form
% % 
% %     u(x,y) = sum_{m=1}^{infty} sum_{n=1}^{infty} B_{mn}*cos((m-0.5)*pi/(b-a)*x)*sin((n-0.5)*pi/(d-c)*y);
% %     u_{t}(x,t) = sum_{m=1}^{infty} sum_{n=1}^{infty} lam_{mn}*C_{mn}*cos((m-0.5)*pi/(b-a)*x)*sin((n-0.5)*pi/(d-c)*y);
% % 
% % where lam_{mn}=C*sqrt(((m-0.5)*pi/(b-a))^2 +((n-0.5)*pi/(d-c))^2)
% Set the coefficients B_mn and C_mn, so that e.g. B_mn(m,n)=B_{mn};
B_mn = [0,2,0;-1,0,0];   C_mn = [0,0,0;0,0.5,0];
% Set the factors lam_{mn} accordingly.
lam_mn = C*sqrt((((1:size(C_mn,1))-0.5)'*sym(pi)/(b-a)).^2 +(((1:size(C_mn,2))-0.5)*sym(pi)/(d-c)).^2);
% Set the spatial basis functions
basis_funs = cos((((1:size(B_mn,1))-0.5)'*pi/(b-a))*sx)*sin((((1:size(B_mn,2))-0.5)*pi/(d-c))*sy);
% Finally, set the initial function.
uinput.ic.PDE(1) = sum(sum(B_mn.*basis_funs));
uinput.ic.PDE(2) = sum(sum(C_mn.*lam_mn.*basis_funs));

% % Given this initial function, the exact solution is given by
% %     u(x,y) = sum_{m=1}^{infty} sum_{n=1}^{infty} sin(m*pi/(b-a)*x)*sin(n*pi/(d-c)*y)*(B_{mn}*cos(lam_{mn}*t)+C_{mn}*sin(lam_{mn}*t));
u_exact = sum(sum(basis_funs.*(B_mn.*cos(lam_mn*st)+C_mn.*sin(lam_mn*st))));
uinput.exact(1) = u_exact;
uinput.exact(2) = diff(u_exact,st);


%----------------------------------------
%% Example 13 - 2D wave equation with Dirichlet-Dirichlet BCs and forcing
%----------------------------------------
case 13
% % Standard 2D wave equation with Dirichlet BCs and forcing
% %   u_tt=C^2*(u_xx + u_yy) +f(x,y,t)                 (x,y) in [0,1]x[0,2]
% %   u(a,y,t)=0; u(b,y,t)=0; u(x,c,t)=0;  u(x,d,t)=0;
% % We want an exact solution
% %   u(x,y,t) = (sin(2*pi*x/(b-a))*sin(2*pi*y/(d-c)))*(1-t)
% % To get this solution, introduce a forcing
% %   f(x,y,t) = u_tt-C^2*(u_xx + u_yy)
% %            = (sin(2*pi*x/(b-a))*sin(2*pi*y/(d-c)))*exp(-t)
% %                 +4*C^2*pi^2*(1/(b-a)+1/(d-c))*(sin(2*pi*x/(b-a))*sin(2*pi*y/(d-c)))*(1-t)
% % And start with initial conditions
% %       u(x,y,0) = 0;
% %   u_{t}(x,y,0) = -(sin(2*pi*x/(b-a))*sin(2*pi*y/(d-c)));

% % Set the spatial domain 
a=-1;   b=1;
c=-1;   d=1;
uinput.dom=[a b;c d];

% Set the wave velocity 
C = 1; 

% % % Construct the PDE.
% % We represent it in the form
% %     x1_{t} = x2;
% %     x2_{t} = C^2*(x1_{xx} +x1_{yy}) +w;
%%% Term-based input format
PDE.x{1}.vars = [s1;s2];   PDE.x{1}.dom = uinput.dom;   % x1 = u;
PDE.x{2}.vars = [s1;s2];   PDE.x{2}.dom = uinput.dom;   % x2 = u_{t}
PDE.x{2}.diff = [2,2];
PDE.w{1}.vars = [s1;s2];   PDE.w{1}.dom = uinput.dom;   % w = f;

% x1_{t} = x2;
PDE.x{1}.term{1}.x = 2;

% x2_{t} = C^2*(x1_{xx} + x1_{yy}) +w;
PDE.x{2}.term{1}.x = 1;
PDE.x{2}.term{1}.D = [2,0;0,2];
PDE.x{2}.term{1}.C = [C^2,C^2];

PDE.x{2}.term{2}.w = 1;

% % Set the boundary conditions
% BC1: 0 = x1(s1,c)                 % BC3: 0 = x1(s1,d)
PDE.BC{1}.term{1}.x = 1;            PDE.BC{3}.term{1}.x = 1;
PDE.BC{1}.term{1}.loc = [s1,c];     PDE.BC{3}.term{1}.loc = [s1,d];
% BC2: 0 = x1(a,s2)                 % BC4: 0 = x1(b,s2)
PDE.BC{2}.term{1}.x = 1;            PDE.BC{4}.term{1}.x = 1;
PDE.BC{2}.term{1}.loc = [a,s2];     PDE.BC{4}.term{1}.loc = [b,s2];

% BC5: 0 = x2(s1,c)                 % BC7: 0 = x2(s1,d)
PDE.BC{5}.term{1}.x = 2;            PDE.BC{7}.term{1}.x = 2;
PDE.BC{5}.term{1}.loc = [s1,c];     PDE.BC{7}.term{1}.loc = [s1,d];
% BC6: 0 = x2(a,s2)                 % BC8: 0 = x2(b,s2)
PDE.BC{6}.term{1}.x = 2;            PDE.BC{8}.term{1}.x = 2;
PDE.BC{6}.term{1}.loc = [a,s2];     PDE.BC{8}.term{1}.loc = [b,s2];


% % Set the desired exact solution
% %     u(x,y,t) = (sin(2*pi*x/(b-a))*sin(2*pi*y/(d-c)))*(1-t)
u_exact = (sin(2*pi*(sx-0.5*(b+a))/(b-a))*sin(2*pi*(sy-0.5*(d+c))/(d-c)))*(1-st);
uinput.exact(1) = u_exact;
uinput.exact(2) = diff(u_exact,st);

% % Set the associated initial conditions
uinput.ic.PDE = subs(uinput.exact,st,0);

% % Set the associated forcing in the PDE
% %     f(x,y,t) = u_tt-C^2*(u_xx + u_yy)
uinput.w(1) = diff(u_exact,st,2) - C^2*(diff(u_exact,sx,2) +diff(u_exact,sy,2));



%----------------------------------------
%% Example 14 - 2D wave equation with Dirichlet-Neumann BCs and forcing
%----------------------------------------
case 14
% % Standard 2D wave equation with Dirichlet-Neumann BCs and forcing
% %   u_tt=C^2*(u_xx + u_yy) +f(x,y,t)                 (x,y) in [0,1]x[0,2]
% %   u(a,y,t)=0; u_{x}(b,y,t)=0; u(x,c,t)=0;  u_{y}(x,d,t)=0;
% % We want an exact solution
% %   u(x,y,t) = sin(1.5*pi*(sx-a)/(b-a))*(sy-c)*(sy-2*d+c);
% % To get this solution, introduce a forcing
% %   f(x,y,t) = u_tt-C^2*(u_xx + u_yy)
% % And start with initial conditions
% %       u(x,y,0) = sin(1.5*pi*(sx-a)/(b-a))*(sy-c)*(sy-2*d+c);
% %   u_{t}(x,y,0) = 0;

% % Set the spatial domain 
a=-1;   b=1;
c=0;   d=1;
uinput.dom=[a b;c d]; 
% NOTE: a, and c must be integer valued for the exact solution to 
% satisfy the boundary conditions.

% Set the wave velocity 
C = 1; 

% % % Construct the PDE.
% % We represent it in the form
% %     x(1)_{t} = x(2);
% %     x(2)_{t} = C^2*(x(1)_{xx} +x(1)_{yy}) +w;
%%% Term-based input format
PDE.x{1}.vars = [s1;s2];   PDE.x{1}.dom = uinput.dom;   % x1 = [u;u_{t}];
PDE.x{1}.size = 2;
PDE.w{1}.vars = [s1;s2];   PDE.w{1}.dom = uinput.dom;   % w = f;

% x_{t} = [0,1,0,0,0,0;0,0,C^2,0,C^2,0]*[x; x_{xx}; x_{yy}] +[0;1]*w;
PDE.x{1}.term{1}.x = 1;
PDE.x{1}.term{1}.C = [0,1,0,0,0,0;
                      0,0,C^2,0,C^2,0];
PDE.x{1}.term{1}.D = [0,0;2,0;0,2];

PDE.x{1}.term{2}.w = 1;
PDE.x{1}.term{2}.C = [0;1];

% % Set the boundary conditions
% BC1: 0 = x1(s1,c)                 % BC3: 0 = x1_{y}(s1,d)
PDE.BC{1}.term{1}.x = 1;            PDE.BC{3}.term{1}.x = 1;
PDE.BC{1}.term{1}.loc = [s1,c];     PDE.BC{3}.term{1}.loc = [s1,d];
                                    PDE.BC{3}.term{1}.D = [0,1];
% BC2: 0 = x1(a,s2)                 % BC4: 0 = x1_{x}(b,s2)
PDE.BC{2}.term{1}.x = 1;            PDE.BC{4}.term{1}.x = 1;
PDE.BC{2}.term{1}.loc = [a,s2];     PDE.BC{4}.term{1}.loc = [b,s2];
                                    PDE.BC{4}.term{1}.D = [1,0];


% % Set the desired exact solution
% %     u(x,y,t) = sin(1.5*pi*(sx-a)/(b-a))*(sy-c)*(sy-2*d+c)
u_exact = sin(1.5*pi*(sx-a)/(b-a))*(sy-c)*(sy-2*d+c);
uinput.exact(1) = u_exact;
uinput.exact(2) = diff(u_exact,st);

% % Set the associated initial conditions
uinput.ic.PDE = subs(uinput.exact,st,0);

% % Set the associated forcing in the PDE
% %     f(x,y,t) = u_tt-C^2*(u_xx + u_yy)
uinput.w(1) = diff(u_exact,st,2) - C^2*(diff(u_exact,sx,2) +diff(u_exact,sy,2));



%----------------------------------------
%% Example 15 - 2D wave equation with Neumann-Neumann BCs
%----------------------------------------
case 15
% % Standard 2D wave equation with Neumann-Neumann BCs
% %   u_tt=C^2*(u_xx + u_yy)                           (x,y) in [0,1]x[0,2]
% %   u_{x}(a,y,t)=0; u_{x}(b,y,t)=0; u_{y}(x,c,t)=0;  u_{y}(x,d,t)=0;

% % Set the spatial domain 
a=0;   b=1;
c=0;   d=2;
uinput.dom=[a b;c d]; 

% Set the wave velocity 
C = 1; 

% % % Construct the PDE.
% % We cannot represent the Neumann boundary conditions direclty as PIE.
% % Instead, introduce U0(t) = int_{a}^{b}int_{c}^{d}u(x,y,t)dy dx,
% % U1(x,t) = int_{c}^{d}u(x,y,t)dy, and U2(y,t) = int_{a}^{b}u(x,y,t)dx.
% % Then
% %     U0_{tt}(t) = int_{a}^{b}int_{c}^{d}u_{tt}(x,y,t)dy dx
% %                = 0;
% % and
% %     U1_{tt}(x,t) = int_{c}^{d}u_{tt}(x,y,t)dy dx
% %                  = C^2*int_{c}^{d}u_{xx}(x,y,t) dy = C^2*U1_{xx}(x,t);
% % and similarly.
% %     U2_{tt}(y,t) = C^2*U2_{yy}(y,t);
% % In addition, we get boundary conditions
% %     int_{a}^{b}U1(x,t)dx = U0(t);       U1_{x}(b,t) = 0;
% %     int_{c}^{d}U2(y,t)dy = U0(t);       U1_{y}(b,t) = 0;
% % Since U0_{tt}(t)=0, we have U0(t)=w1 +w2*t for some constants
% % w1 = U0(0)=  int_{a}^{b}int_{c}^{d}u(x,y,0)dy dx and
% % w2 = U0_{t}(0) = int_{a}^{b}int_{c}^{d}u_{t}(x,y,0)dy dx.
% %
% % Finally, to deal with the second order temporal derivative,
% % augment the state as u--> [u; u_{t}];

%%% Declare the PDE
PDE.x{1}.vars = s1;         PDE.x{1}.dom = [a,b];       % x1 = [U1(x,t);U1_{x}(x,t)]
PDE.x{1}.size = 2;
PDE.x{2}.vars = s2;         PDE.x{2}.dom = [c,d];       % x2 = [U2(y,t);U2_{y}(y,t)]
PDE.x{2}.size = 2;
PDE.x{3}.vars = [s1;s2];    PDE.x{3}.dom = [a,b;c,d];   % x3 = [u;u_{t}];
PDE.x{3}.size = 2;

PDE.w{1}.vars = [];
PDE.w{1}.size = 2;

% [U1_{t};U1_{tt}] = [x1(2); C^2*x1_{xx}(1)]
PDE.x{1}.term{1}.x = 1;
PDE.x{1}.term{1}.C = [0,1,0,0;0,0,C^2,0];
PDE.x{1}.term{1}.D = [0;2];

% [U2_{t};U2_{tt}] = [x2(2); C^2*x2_{yy}(1)]
PDE.x{2}.term{1}.x = 2;
PDE.x{2}.term{1}.C = [0,1,0,0;0,0,C^2,0];
PDE.x{2}.term{1}.D = [0;2];

% [u_{t};u_{tt}] = [0,1,0,0,0,0;0,0,C^2,0,C^2,0]*[x3(2); x3_{xx}(1); x3_{yy}(1)];
PDE.x{3}.term{1}.x = 3;
PDE.x{3}.term{1}.C = [0,1,0,0,0,0;
                      0,0,C^2,0,C^2,0];
PDE.x{3}.term{1}.D = [0,0;2,0;0,2];

% % Set the boundary conditions
% int_{a}^{b}U1(x,t)dx = w(t);      U1_{x}(b,t) = 0;
PDE.BC{1}.term{1}.x = 1;            PDE.BC{2}.term{1}.x = 1;
PDE.BC{1}.term{1}.I = {[a,b]};      PDE.BC{2}.term{1}.loc = b;
PDE.BC{1}.term{2}.w = 1;            PDE.BC{2}.term{1}.D = 1;
PDE.BC{1}.term{2}.C = -1*eye(2);

% int_{c}^{d}U2(y,t)dy = w(t);      U1_{y}(b,t) = 0;
PDE.BC{3}.term{1}.x = 2;            PDE.BC{4}.term{1}.x = 2;
PDE.BC{3}.term{1}.I = {[c,d]};      PDE.BC{4}.term{1}.loc = d;
PDE.BC{3}.term{2}.w = 1;            PDE.BC{4}.term{1}.D = 1;
PDE.BC{3}.term{2}.C = -1*eye(2);

% int_{c}^{d}u(x,y,t)dy = U1(t,x);  u_{y}(x,d,t) = 0;
PDE.BC{5}.term{1}.x = 3;            PDE.BC{6}.term{1}.x = 3;
PDE.BC{5}.term{1}.I = {[];[c,d]};   PDE.BC{6}.term{1}.loc = [s1,d];
PDE.BC{5}.term{2}.x = 1;            PDE.BC{6}.term{1}.D = [0,1];
PDE.BC{5}.term{2}.C = -1*eye(2);

% int_{a}^{b}u(x,y,t)dx = U2(t,y);  u_{x}(b,y,t) = 0;
PDE.BC{7}.term{1}.x = 3;            PDE.BC{8}.term{1}.x = 3;
PDE.BC{7}.term{1}.I = {[a,b];[]};   PDE.BC{8}.term{1}.loc = [b,s2];
PDE.BC{7}.term{2}.x = 2;            PDE.BC{8}.term{1}.D = [1,0];
PDE.BC{7}.term{2}.C = -1*eye(2);


% % Set the initial conditions
% % We build an initial function of the form
% % 
% %     u(x,y) = sum_{m=0}^{infty} sum_{n=0}^{infty} B_{mn}*cos(m*pi/(b-a)*x)*cos(n*pi/(d-c)*y);
% %     u_{t}(x,t) = sum_{m=0}^{infty} sum_{n=0}^{infty} lam_{mn}*C_{mn}*cos(m*pi/(b-a)*x)*cos(n*pi/(d-c)*y);
% % 
% % where lam_{mn}=C*sqrt((m*pi/(b-a))^2 +(n*pi/(d-c))^2)
% Set the coefficients B_mn and C_mn, so that e.g. B_mn(m,n)=B_{(m-1)(n-1)};
B_mn = [5,10,0;-10,0,0];   C_mn = [0,0,-2;0,0,3];
% Set the factors lam_{mn} accordingly.
lam_mn = C*sqrt(((0:size(C_mn,1)-1)'*sym(pi)/(b-a)).^2 +((0:size(C_mn,2)-1)*sym(pi)/(d-c)).^2);
% Set the spatial basis functions
basis_funs = cos(((0:size(B_mn,1)-1)'*sym(pi)/(b-a))*(sx-a))*cos(((0:size(B_mn,2)-1)*sym(pi)/(d-c))*(sy-c));
% Finally, set the initial function.
u_ic = [sum(sum(B_mn.*basis_funs)), sum(sum(C_mn.*lam_mn.*basis_funs))];
uinput.ic.PDE([1,2]) = int(u_ic,sy,c,d);
uinput.ic.PDE([3,4]) = int(u_ic,sx,a,b);
uinput.ic.PDE([5,6]) = u_ic;

% % Set the "forcing" U0(t) = [w1+w2*t; w2];
w_params = int(int(u_ic,sx,a,b),sy,c,d);
uinput.w([1,2]) = [w_params(1)+w_params(2)*st, w_params(2)];

% % Given this initial function, the exact solution is given by
% %     u(x,y) = sum_{m=0}^{infty} sum_{n=0}^{infty} sin(m*pi/(b-a)*x)*sin(n*pi/(d-c)*y)*(B_{mn}*cos(lam_{mn}*t)+C_{mn}*sin(lam_{mn}*t));
u_ex1 = sum(sum(basis_funs.*(B_mn.*cos(lam_mn*st)+C_mn.*sin(lam_mn*st))));
u_ex = [u_ex1, diff(u_ex1,st)];
uinput.exact([1,2]) = int(u_ex,sy,c,d);
uinput.exact([3,4]) = int(u_ex,sx,a,b);
uinput.exact([5,6]) = u_ex;

%----------------------------------------
%% Example 16 - 2D wave equation with Periodic BCs along x-direction
%----------------------------------------
case 16
% % Standard 2D wave equation with Periodic BCs along x-direction
% %   u_tt=C^2*(u_xx + u_yy)                          (x,y) in [-1,1]x[0,1]
% %   u(a,y,t)=u(b,y,t);    u_{x}(a,y,t)=u_{x}(b,y,t); 
% %   u(x,c,t)=0;           u_{y}(x,d,t)=0;
% % With initial conditions
% %     u(x,y) = ampl0*sin(0.5*pi/(d-c)*y) +ampl1*sin(freq*pi*x-shft)*sin(2.5*pi/(d-c)*y);
% %     u_{t}(x,y) = ampl2*C*sqrt((freq*pi)^2 +(2.5*pi/(d-c))^2)*sin(freq*pi*x-shft)*sin(2.5*pi/(d-c)*y);
% % We get an exact solution
% %     u(x,y) = ampl0*sin(0.5*pi/(d-c)*y)*cos(C*0.5*pi/(d-c)*t) +sin(freq*pi*x-shft)*sin(2.5*pi/(d-c)*y)*(ampl1*cos(lam*t)+ampl2*sin(lam*t));

% % Set the spatial domain 
a=-1;   b=1;
c=0;   d=1;
uinput.dom=[a b;c d]; 

% Set the wave velocity, and initial condition parameters.
%C = 1; 
%ampl0 = 3;      ampl1 = 2;     ampl2 = 1;
%freq = 3;       shft = 0.25*pi;

C = 1; 
ampl0 = 1;      ampl1 = 1;     ampl2 = 0;
freq = 1;       shft = 0;


% % % Construct the PDE.
% % We cannot represent the peridoic boundary conditions direclty as PIE.
% % Instead, introduce U1(y,t) = int_{a}^{b}u(x,y,t) dx, then
% %     U1_{tt}(y,t) = int_{a}^{b}u_{tt}(x,y,t) dx
% %                  = C^2*int_{a}^{b}u_{yy}(x,y,t) dx = C^2*U1_{yy}(y,t);
% % In addition, we get boundary conditions
% %     U1(c,t) = 0;        U1_{y}(d,t) = 0;
% % To deal with the second order temporal derivative,
% % augment the state as u--> [u; u_{t}];

%%% Declare the PDE
PDE.x{1}.vars = s2;         PDE.x{1}.dom = [c,d];       % x1 = [U1(y,t);U1_{y}(y,t)]
PDE.x{1}.size = 2;
PDE.x{2}.vars = [s1;s2];    PDE.x{2}.dom = [a,b;c,d];   % x2 = [u;u_{t}];
PDE.x{2}.size = 2;

% [U1_{t};U1_{tt}] = [x1(2); C^2*x1_{yy}(1)]
PDE.x{1}.term{1}.x = 1;
PDE.x{1}.term{1}.C = [0,1,0,0;0,0,C^2,0];
PDE.x{1}.term{1}.D = [0;2];

% [u_{t};u_{tt}] = [0,1,0,0,0,0;0,0,C^2,0,C^2,0]*[x3(2); x3_{xx}(1); x3_{yy}(1)];
PDE.x{2}.term{1}.x = 2;
PDE.x{2}.term{1}.C = [0,1,0,0,0,0;
                      0,0,C^2,0,C^2,0];
PDE.x{2}.term{1}.D = [0,0;2,0;0,2];

% % Set the boundary conditions
% U1(c,t) = 0;                      U1_{y}(d,t) = 0;
PDE.BC{1}.term{1}.x = 1;            PDE.BC{2}.term{1}.x = 1;
PDE.BC{1}.term{1}.loc = c;          PDE.BC{2}.term{1}.loc = d;
                                    PDE.BC{2}.term{1}.D = 1;

% int_{a}^{b}u(x,y,t)dx = U1(y,t);  u(a,y,t) = u(b,y,t);
PDE.BC{3}.term{1}.x = 2;            PDE.BC{4}.term{1}.x = 2;
PDE.BC{3}.term{1}.I = {[a,b];[]};   PDE.BC{4}.term{1}.loc = [a,s2];
PDE.BC{3}.term{2}.x = 1;            PDE.BC{4}.term{2}.x = 2;
PDE.BC{3}.term{2}.C = -1*eye(2);    PDE.BC{4}.term{2}.loc = [b,s2];
                                    PDE.BC{4}.term{2}.C = -1*eye(2);

% u(x,c,t) = 0;                     u_{y}(x,d,t) = 0;
PDE.BC{5}.term{1}.x = 2;            PDE.BC{6}.term{1}.x = 2;
PDE.BC{5}.term{1}.loc = [s1,c];     PDE.BC{6}.term{1}.loc = [s1,d];
                                    PDE.BC{6}.term{1}.D = [0,1];


% % Set the initial conditions
% %     u(x,y) = ampl0*sin(0.5*pi/(d-c)*y) +ampl1*sin(freq*pi*x-shft)*sin(2.5*pi/(d-c)*y);
% %     u_{t}(x,t) = ampl2*C*sqrt((freq*pi)^2 +(2.5*pi/(d-c))^2)*sin(freq*pi*x-shft)*sin(2.5*pi/(d-c)*y);
% % 
lam = C*sqrt((freq*pi)^2 +(2.5*pi/(d-c))^2);
u_ic1 = ampl0*sin(0.5*sym(pi)/(d-c)*sy) +ampl1*sin(freq*sym(pi)*sx-shft)*sin(2.5*sym(pi)/(d-c)*sy);
u_ic2 = ampl2*lam*sin(freq*sym(pi)*sx-shft)*sin(2.5*sym(pi)/(d-c)*sy);
u_ic = [u_ic1, u_ic2];
uinput.ic.PDE([1,2]) = int(u_ic,sx,a,b);
uinput.ic.PDE([3,4]) = u_ic;

% % Given this initial function, the exact solution is given by
% %     u(x,y) = ampl0*sin(0.5*pi/(d-c)*y)*cos(C*0.5*pi/(d-c)*t) +sin(freq*pi*x-shft)*sin(2.5*pi/(d-c)*y)*(ampl1*cos(lam*t)+ampl2*sin(lam*t));
u_ex1 = ampl0*sin(0.5*sym(pi)/(d-c)*sy)*cos(C*0.5*sym(pi)/(d-c)*st) +sin(freq*sym(pi)*sx-shft)*sin(2.5*sym(pi)/(d-c)*sy)*(ampl1*cos(lam*st)+ampl2*sin(lam*st));
u_ex = [u_ex1, diff(u_ex1,st)];
uinput.exact([1,2]) = int(u_ex,sx,a,b);
uinput.exact([3,4]) = u_ex;


%----------------------------------------
%% Example 17 - Coupled ODE - 1D PDE - 2D PDE system
%----------------------------------------
case 17
% % Coupled ODE - 1D PDE - 2D PDE system
% %   X_{tt}(t)      = -C^2*X(t);
% %   u1_{tt}(x,t)   = C^2*(b-a)^2/pi^2*U1_{xx}
% %   u2_{tt}(y,t)   = C^2*(d-c)^2/(25*pi^2)*U2_{yy}
% %   u3_{tt}(x,y,t) = 0.5*C^2/pi^2*((b-a)^2*u_{xx}+(d-c)^2/25 *u_{yy})    (x,y) in [a,b]x[c,d]
% % with BCs
% %   u1(a,t) = X(t)    u1_{x}(b,t) = 0;
% %   u2(c,t) = X(t)    u2_{y}(d,t) = 0;
% %   u3(x,c,t) = u1(x,t);      u3_{y}(x,d,t) = 0;
% %   u3(a,y,t) = u2(y,t);      u3_{x}(b,y,t) = 0;
% % Starting with initial conditions
% %   X(0)      = ampl;                             X_{t}(0)=0;
% %   u1(x,0)   = ampl*cos(pi*(x-a)/(b-a));         u1_{t}(x,0)=0;
% %   u2(y,0)   = ampl*cos(5*pi*(y-c)/(d-c));       u2_{t}(y,0)=0;
% %   u3(x,y,0) = ampl*cos(pi*(x-a)/(b-a))*cos(5*pi*(y-c)/(d-c));
% %     u3_{t}(x,y,0) = 0
% % We get exact solution
% %   X(t)      = ampl*cos(C*t);
% %   u1(x,t)   = ampl*cos(pi*(x-a)/(b-a))*cos(C*t);
% %   u2(y,t)   = ampl*cos(5*pi*(y-c)/(d-c))*cos(C*t);
% %   u3(x,y,t) = ampl*cos(pi*(x-a)/(b-a))*cos(5*pi*(y-c)/(d-c))*cos(C*t);
% 

% % Set the spatial domain 
a=0;   b=1;
c=0;   d=5;

%a=-1;   b=1;
%c=-1;   d=1;

uinput.dom=[a b;c d]; 

% Set the wave velocity and amplitude
C = 0.5;       ampl = 10;

% % Declare the PDE
PDE.x{1}.vars = [];
PDE.x{1}.size = 2;
PDE.x{2}.vars = s1;         PDE.x{2}.dom = [a,b];
PDE.x{2}.size = 2;
PDE.x{3}.vars = s2;         PDE.x{3}.dom = [c,d];
PDE.x{3}.size = 2;
PDE.x{4}.vars = [s1;s2];    PDE.x{4}.dom = [a,b;c,d];   % x{4}=[u3;u3_{t}]
PDE.x{4}.size = 2;

% X_{tt}(t) = -C^2*X(t);
PDE.x{1}.term{1}.x = 1;     PDE.x{1}.term{1}.C = [0,1;-C^2,0];

% u1_{tt}(x,t)   = C^2*(b-a)^2/pi^2*u1_{xx}
PDE.x{2}.term{1}.x = 2;     PDE.x{2}.term{1}.C = [0,1,0,0; 0,0,C^2*(b-a)^2/pi^2,0];
PDE.x{2}.term{1}.D = [0;2];

% u2_{t}(y,t)   = C^2*(d-c)^2/(25*pi^2)*u2_{yy}
PDE.x{3}.term{1}.x = 3;     PDE.x{3}.term{1}.C = [0,1,0,0; 0,0,C^2*(d-c)^2/(25*pi^2),0];
PDE.x{3}.term{1}.D = [0;2];

% u3_{t}(x,y,t) = 0.5*C^2/pi^2*((b-a)^2*u_{xx}+(d-c)^2/25 *u_{yy})
PDE.x{4}.term{1}.x = 4;     PDE.x{4}.term{1}.C = [0,1,0,0,0,0; 0,0,0.5*C^2/pi^2*(b-a)^2,0,0.5*C^2/pi^2*(d-c)^2/25,0];
PDE.x{4}.term{1}.D = [0,0;2,0;0,2];

% % BCs
% u1(a,t) = X(t);                       u1_{x}(b,t) = 0;
PDE.BC{1}.term{1}.x = 1;                PDE.BC{2}.term{1}.x = 2;
PDE.BC{1}.term{1}.C = -1*eye(2);
PDE.BC{1}.term{2}.x = 2;                PDE.BC{2}.term{1}.loc = b;
PDE.BC{1}.term{2}.loc = a;              PDE.BC{2}.term{1}.D = 1;

% u2(c,t) = X(t);                       u2_{y}(d,t) = 0;
PDE.BC{3}.term{1}.x = 1;                PDE.BC{4}.term{1}.x = 3;
PDE.BC{3}.term{1}.C = -1*eye(2);
PDE.BC{3}.term{2}.x = 3;                PDE.BC{4}.term{1}.loc = d;
PDE.BC{3}.term{2}.loc = c;              PDE.BC{4}.term{1}.D = 1;

% u3(x,c,t) = u1(x,t);                  u3(a,y,t) = u2(y,t); 
PDE.BC{5}.term{1}.x = 2;                PDE.BC{6}.term{1}.x = 3;
PDE.BC{5}.term{1}.C = -1*eye(2);        PDE.BC{6}.term{1}.C = -1*eye(2);
PDE.BC{5}.term{2}.x = 4;                PDE.BC{6}.term{2}.x = 4;
PDE.BC{5}.term{2}.loc = [s1,c];         PDE.BC{6}.term{2}.loc = [a,s2];

% u3_{x}(b,y,t) = 0;                    u3_{y}(x,d,t) = 0;
PDE.BC{7}.term{1}.x = 4;                PDE.BC{8}.term{1}.x = 4;
PDE.BC{7}.term{1}.D = [1,0];            PDE.BC{8}.term{1}.D = [0,1];
PDE.BC{7}.term{1}.loc = [b,s2];         PDE.BC{8}.term{1}.loc = [s1,d];


% % Set the initial conditions
% %   u(x,y,0) = ampl*cos(pi*(x-a)/(b-a))*cos(5*pi*(y-c)/(d-c));
u_ic1 = ampl*cos(sym(pi)*(sx-a)/(b-a))*cos(5*sym(pi)*(sy-c)/(d-c));
u_ic2 = 0;
u_ic = [u_ic1,u_ic2];
uinput.ic.x([1,2]) = subs(subs(u_ic,sy,c),sx,a);
uinput.ic.x([3,4]) = subs(u_ic,sy,c);
uinput.ic.x([5,6]) = subs(u_ic,sx,a);
uinput.ic.x([7,8]) = u_ic;

% % Set the exact solution
% %   u(x,y,t) = ampl*cos(pi*(x-a)/(b-a))*cos(5*pi*(y-c)/(d-c))*cos(C*t)
u_ex1 = ampl*cos(sym(pi)*(sx-a)/(b-a))*cos(5*sym(pi)*(sy-c)/(d-c))*cos(C*st);
u_ex2 = diff(u_ex1,st);
u_ex = [u_ex1,u_ex2];
uinput.exact([1,2]) = subs(subs(u_ex,sy,c),sx,a);
uinput.exact([3,4]) = subs(u_ex,sy,c);
uinput.exact([5,6]) = subs(u_ex,sx,a);
uinput.exact([7,8]) = u_ex;



%----------------------------------------
%% Example 18 - 2D transport
%----------------------------------------
case 18
% % Transport equation along x- and y-direction.
% %   u_{t}(x,y,t) = -v*(u_{x}(x,y,t)+u_{y}(x,y,t))      (x,y) in [0,b]x[0,d]
% %   u1_{t}(x,t) = -v*u1_{x}(x,y,t) -w2;    
% %   u2_{t}(y,t) = -v*u2_{y}(x,y,t) -w2;
% % with BCs
% %   u(0,y,t) = u2(y,t);     u(x,0,t) = u1(x,t);
% %   u1(0,t) = w1(t);         u2(0,t) = w1(t)
% % Starting with an initial condition
% %   u(x,y,0) = ampl*(sin(x)+sin(y));
% %   u1(x,0) = ampl*sin(x);        u2(y,0) = ampl*sin(y);
% % And using forcing
% %   w1(t) = 2*ampl*sin(-v*t);     w2 = ampl*v*cos(v*t)
% % we get an exact solution
% %   u(x,y,t) = ampl*(sin(x-v*t) +sin(y-v*t));
% %   u1(x,t) = ampl*(sin(x-v*t) +sin(-v*t));     
% %   u2(y,t) = ampl*(sin(y-v*t) +sin(-v*t));

% Set the spatial domain 
a=0;   b=2*pi;
c=0;   d=2*pi;


uinput.dom=[a b;c d]; 

% Set the wave velocity and amplitude
v = 0.5;       ampl = 5;

% % Declare the PDE
PDE.x{1}.vars = [s1;s2];    PDE.x{1}.dom = [a,b;c,d];
PDE.x{2}.vars = s1;         PDE.x{2}.dom = [a,b];
PDE.x{3}.vars = s2;         PDE.x{3}.dom = [c,d];
PDE.w{1}.vars = [];
PDE.w{2}.vars = [];

% u_{t}(t) = -v*(u_{x}(t) +u_{y}(t));
PDE.x{1}.term{1}.x = 1;     PDE.x{1}.term{1}.C = [-v,-v];
PDE.x{1}.term{1}.D = [1,0;0,1];

% u1_{t}(t) = -v*u1_{x}(x,y,t) -w2;  
PDE.x{2}.term{1}.x = 2;     PDE.x{2}.term{1}.C = -v;
PDE.x{2}.term{1}.D = 1;
PDE.x{2}.term{2}.w = 2;     PDE.x{2}.term{2}.C = -1;

% u2_{t}(t) = -v*u2_{y}(x,y,t) -w2;  
PDE.x{3}.term{1}.x = 3;     PDE.x{3}.term{1}.C = -v;
PDE.x{3}.term{1}.D = 1;
PDE.x{3}.term{2}.w = 2;     PDE.x{3}.term{2}.C = -1;

% % Declare the BCs
% u(0,y,t) = u2(y,t)                u(x,0,t) = u1(x,t);
PDE.BC{1}.term{1}.x = 1;            PDE.BC{2}.term{1}.x = 1;            
PDE.BC{1}.term{1}.loc = [a,s2];     PDE.BC{2}.term{1}.loc = [s1,c];
PDE.BC{1}.term{2}.x = 3;            PDE.BC{2}.term{2}.x = 2;
PDE.BC{1}.term{2}.C = -1;           PDE.BC{2}.term{2}.C = -1;

% u1(0,t) = w1(t)                   % u2(0,t) = w1(t)
PDE.BC{3}.term{1}.x = 2;            PDE.BC{4}.term{1}.x = 3;
PDE.BC{3}.term{1}.loc = a;          PDE.BC{4}.term{1}.loc = c;
PDE.BC{3}.term{2}.w = 1;            PDE.BC{4}.term{2}.w = 1;
PDE.BC{3}.term{2}.C = -1;           PDE.BC{4}.term{2}.C = -1;

% % Set the initial conditions
% %   u(x,y,t) = ampl*(sin(x)+sin(y));
% %   u1(x,t) = ampl*sin(x);        u2(y,t) = ampl*sin(y);

u_ic = ampl*(sin(sx) +sin(sy));
uinput.ic.x(1) = u_ic;
uinput.ic.x(2) = subs(u_ic,sy,c);
uinput.ic.x(3) = subs(u_ic,sx,a);

% % Set the forcing
% %   w1(t) = 2*ampl*sin(-v*t);     w2 = ampl*v*cos(-v*t)
uinput.w(1) = 2*ampl*sin(-v*st);    uinput.w(2) = ampl*v*cos(-v*st);

% % Set the exact solution
% %   u(x,y,t) = ampl*(sin(x-v*t) +sin(y-v*t));
% %   u1(x,t) = ampl*(sin(x-v*t) +sin(-v*t));     
% %   u2(y,t) = ampl*(sin(y-v*t) +sin(-v*t));

u_ex = ampl*(sin(sx-v*st) +sin(sy-v*st));

uinput.exact(1) = u_ex;
uinput.exact(2) = subs(u_ex,sy,c);
uinput.exact(3) = subs(u_ex,sx,a);

%----------------------------------------
%% Example 19 - Mixed diffusive equation
%----------------------------------------
case 19
% % Mixed diffusive equation.
% %   u_{t}(x,y,t) = lam*u_{xy}(x,y,t)                 (x,y) in [0,2.5*pi]x[0,2.5*pi]
% %    u1_{t}(x,t) = lam*u1_{xx}(x,t)
% %    u2_{t}(y,t) = lam*u2_{yy}(y,t)
% % with Dirichlet BCs
% %   u(0,y,t) = u2(y,t);       u(x,0,t) = u1(x,t);
% %   u1(0,t) = 0;              u1_{x}(2.5*pi,t) = 0;
% %   u2(0,t) = 0;              u2_{y}(2.5*pi,t) = 0;
% % Starting with an initial condition
% %   u(x,y,t) = ampl*sin(x+y);
% %   u1(x,t) = ampl*sin(x);        u2(y,t) = ampl*sin(y);
% % we get an exact solution
% %   u(x,y,t) = ampl*sin(x+y)*exp(-lam*t);
% %   u1(x,t) = ampl*sin(x)*exp(-lam*t);
% %   u2(y,t) = ampl*sin(y)*exp(-lam*t);

% % Set the spatial domain 
a=0;   b=2.5*pi;
c=0;   d=2.5*pi;
uinput.dom=[a b;c d]; 

% Set the wave velocity and amplitude
lam = 0.5;       ampl = 5;

% % Declare the PDE
PDE.x{1}.vars = [s1;s2];    PDE.x{1}.dom = [a,b;c,d];
PDE.x{2}.vars = s1;         PDE.x{2}.dom = [a,b];
PDE.x{3}.vars = s2;         PDE.x{3}.dom = [c,d];

% u_{t}(t) = lam*u_{xy}(t);
PDE.x{1}.term{1}.x = 1;     PDE.x{1}.term{1}.C = lam;
PDE.x{1}.term{1}.D = [1,1];

% u1_{t}(x,t) = lam*u1_{xx}(x,t)
PDE.x{2}.term{1}.x = 2;     PDE.x{2}.term{1}.C = lam;
PDE.x{2}.term{1}.D = 2;

% u2_{t}(y,t) = lam*u2_{yy}(y,t)
PDE.x{3}.term{1}.x = 3;     PDE.x{3}.term{1}.C = lam;
PDE.x{3}.term{1}.D = 2;

% % Declare the BCs
% u(0,y,t) = u2(y,t)                u(x,0,t) = u1(x,t);
PDE.BC{1}.term{1}.x = 1;            PDE.BC{2}.term{1}.x = 1;            
PDE.BC{1}.term{1}.loc = [a,s2];     PDE.BC{2}.term{1}.loc = [s1,c];
PDE.BC{1}.term{2}.x = 3;            PDE.BC{2}.term{2}.x = 2;
PDE.BC{1}.term{2}.C = -1;           PDE.BC{2}.term{2}.C = -1;

% u1(0,t) = 0;                      u1_{x}(2.5*pi,t) = 0;
PDE.BC{3}.term{1}.x = 2;            PDE.BC{4}.term{1}.x = 2;
PDE.BC{3}.term{1}.loc = a;          PDE.BC{4}.term{1}.loc = b;
                                    PDE.BC{4}.term{1}.D = 1;

% u2(0,t) = 0;                      u2_{y}(2.5*pi,t) = 0;
PDE.BC{5}.term{1}.x = 3;            PDE.BC{6}.term{1}.x = 3;
PDE.BC{5}.term{1}.loc = c;          PDE.BC{6}.term{1}.loc = d;
                                    PDE.BC{6}.term{1}.D = 1;

% % Set the initial conditions
% %   u(x,y,t) = ampl*sin(x+y);
% %   u1(x,t) = ampl*sin(x);        u2(y,t) = ampl*sin(y);
uinput.ic.x(1) = ampl*sin(sx+sy);
uinput.ic.x(2) = ampl*sin(sx);      uinput.ic.x(3) = ampl*sin(sy);

% % Set the exact solution
% %   u(x,y,t) = ampl*sin(x+y)*exp(-lam*t);
% %   u1(x,t) = ampl*sin(x)*exp(-lam*t);
% %   u2(y,t) = ampl*sin(y)*exp(-lam*t);
uinput.exact(1) = ampl*sin(sx+sy)*exp(-lam*st);
uinput.exact(2) = ampl*sin(sx)*exp(-lam*st);
uinput.exact(3) = ampl*sin(sy)*exp(-lam*st);



%----------------------------------------
%% Example 20 - Heat Equation with Dirichlet Boundary Forcing
%----------------------------------------
case 20
% % Heat equation
% %   u_{t}(x,y,t) = visc*(u_{xx}(x,y,t) +u_{yy}(x,y,t))    (x,y) in [0,1]^2
% % with Dirichlet BCs
% %   u(0,y,t) = w(y,t);        u(1,y,t) = 0;
% %   u(x,0,t) = 0;             u(x,1,t) = 0;
% % Using forcing
% %   w(y,t) = ampl*sin(2*pi*y)*exp(-visc*((2pi)^2+(1.5pi)^2)*t);
% % Starting with an initial condition
% %   u(x,y,0) = ampl*sin(2*pi*y)*cos(1.5*pi*x);
% % we get an exact solution
% %   u(x,y,t) = ampl*sin(2*pi*y)*cos(1.5*pi*x)*exp(-visc*((2pi)^2+(1.5pi)^2)*t);

% % Set the spatial domain (exact solution is only valid for [0,1]^2 !!!)
a=0;   b=1;
c=0;   d=1;
uinput.dom=[a b;c d]; 

% Set the viscosity and amplitude
visc = 0.5;       ampl = 2;

% Initialize the state and forcing
PDE.x{1}.vars = [s1;s2];    PDE.x{1}.dom = [a,b;c,d];
PDE.w{1}.vars = s2;         PDE.w{1}.dom = [c,d];

% % Declare the PDE dynamics
% u_{t}(x,y,t) = visc*(u_xx(x,t) + u_yy(y,t))
PDE.x{1}.term{1}.x = 1;     PDE.x{1}.term{1}.C = [visc,visc];
PDE.x{1}.term{1}.D = [2,0;0,2];

% % BCs
% u(0,y,t) = w(y,t);                    u(1,y,t) = 0;
PDE.BC{1}.term{1}.x = 1;                PDE.BC{2}.term{1}.x = 1;
PDE.BC{1}.term{1}.loc = [a,s2];         PDE.BC{2}.term{1}.loc = [b,s2];
PDE.BC{1}.term{2}.w = 1;
PDE.BC{1}.term{2}.C = -1;

% u(x,0,t) = 0;                         u(x,1,t) = 0;
PDE.BC{3}.term{1}.x = 1;                PDE.BC{4}.term{1}.x = 1;
PDE.BC{3}.term{1}.loc = [s1,c];         PDE.BC{4}.term{1}.loc = [s1,d];


% % Set the forcing
% %   w(y,t) = sin(2*pi*y)*exp(-visc*((2pi)^2+(1.5pi)^2)*t);
uinput.w(1) = ampl*sin(2*sym(pi)*sy)*exp(-visc*((2*sym(pi))^2+(1.5*sym(pi))^2)*st);

% % Set the initial conditions
% %   u(x,y,0) = ampl*sin(2*pi*y)*cos(1.5*pi*x);
u_ic = ampl*sin(2*sym(pi)*sy)*cos(1.5*sym(pi)*sx);
uinput.ic.x(1) = u_ic;

% % Set the exact solution
% %   u(x,y,t) = ampl*sin(2*pi*y)*cos(1.5*pi*x)*exp(-visc*((2pi)^2+(1.5pi)^2)*t);
u_ex = ampl*sin(2*sym(pi)*sy)*cos(1.5*sym(pi)*sx)*exp(-visc*((2*sym(pi))^2+(1.5*sym(pi))^2)*st);
uinput.exact(1) = u_ex;



%----------------------------------------
%% Example 21 - Heat Equation with Neumann Boundary Forcing
%----------------------------------------
case 21
% % Heat equation
% %   u_{t}(x,y,t) = visc*(u_{xx}(x,y,t) +u_{yy}(x,y,t))    (x,y) in [0,1]^2
% % with homogeneous Dirichlet BCs but Neumann forcing
% %   u_{x}(0,y,t) = w(y,t);    u(1,y,t) = 0;
% %   u(x,0,t) = 0;             u(x,1,t) = 0;
% % Using forcing
% %   w(y,t) = 2*pi*ampl*sin(2*pi*y)*exp(-visc*8*pi^2*t);
% % Starting with an initial condition
% %   u(x,y,0) = ampl*sin(2*pi*y)*sin(2*pi*x);
% % we get an exact solution
% %   u(x,y,t) = ampl*sin(2*pi*y)*sin(2*pi*x)*exp(-visc*8*pi^2*t);

% % Set the spatial domain (exact solution is only valid for [0,1]^2 !!!)
a=0;   b=1;
c=0;   d=1;
uinput.dom=[a b;c d]; 

% Set the viscosity and amplitude
visc = 0.5;       ampl = 10;

% Initialize the state and forcing
PDE.x{1}.vars = [s1;s2];    PDE.x{1}.dom = [a,b;c,d];
PDE.w{1}.vars = s2;         PDE.w{1}.dom = [c,d];

% % Declare the PDE dynamics
% u_{t}(x,y,t) = visc*(u_xx(x,t) + u_yy(y,t))
PDE.x{1}.term{1}.x = 1;     PDE.x{1}.term{1}.C = [visc,visc];
PDE.x{1}.term{1}.D = [2,0;0,2];

% % BCs
% u_{x}(0,y,t) = w(y,t);                u(1,y,t) = 0;
PDE.BC{1}.term{1}.x = 1;                PDE.BC{2}.term{1}.x = 1;
PDE.BC{1}.term{1}.loc = [a,s2];         PDE.BC{2}.term{1}.loc = [b,s2];
PDE.BC{1}.term{1}.D = [1,0];
PDE.BC{1}.term{2}.w = 1;
PDE.BC{1}.term{2}.C = -1;

% u(x,0,t) = 0;                         u(x,1,t) = 0;
PDE.BC{3}.term{1}.x = 1;                PDE.BC{4}.term{1}.x = 1;
PDE.BC{3}.term{1}.loc = [s1,c];         PDE.BC{4}.term{1}.loc = [s1,d];


% % Set the forcing
% %   w(y,t) = 2*pi*ampl*sin(2*pi*y)*exp(-visc*8*pi^2*t);
uinput.w(1) = 2*sym(pi)*ampl*sin(2*sym(pi)*sy)*exp(-visc*8*sym(pi)^2*st);

% % Set the initial conditions
% %   u(x,y,0) = ampl*sin(2*pi*y)*sin(2*pi*x);
u_ic = ampl*sin(2*sym(pi)*sy)*sin(2*sym(pi)*sx);
uinput.ic.x(1) = u_ic;

% % Set the exact solution
% %   u(x,y,t) = ampl*sin(2*pi*y)*sin(2*pi*x)*exp(-visc*8*pi^2*t);
u_ex = ampl*sin(2*sym(pi)*sy)*sin(2*sym(pi)*sx)*exp(-visc*8*sym(pi)^2*st);
uinput.exact(1) = u_ex;



%----------------------------------------
%% Example 22 - Heat Equation with Multiple Boundary Forcing
%----------------------------------------
case 22
% % Heat equation
% %   u_{t}(x,y,t) = u_{xx}(x,y,t) +u_{yy}(x,y,t) +(2.5*pi)^2*u(x,y,t)    (x,y) in [-1,1]x[0,1]
% % with forcing along all boundaries
% %   u(-1,y,t) = w1(y,t);      u(y,1,t) = w2(y,t);
% %   u(x,0,t) = w3(x,t);       u_{y}(x,1,t) = w4(y,t);
% % Using forcing
% %   w1(y,t) = ampl*cos(2.5*pi*y);
% %   w2(y,t) = 3*ampl*cos(2.5*pi*y);
% %   w3(x,t) = ampl*(x+2);     
% %   w4(x,t) = -2.5*pi*ampl*(x+2);
% % Starting with an initial condition
% %   u(x,y,0) = ampl*sin(2*pi*y)*sin(2*pi*x);
% % we get an exact solution
% %   u(x,y,t) = ampl*(x+2)*cos(2.5*pi*y);

% % Set the spatial domain (exact solution is only valid for [-1,1]x[0,1] !!!)
a=-1;   b=1;
c=0;   d=1;
uinput.dom=[a b;c d]; 

% Set the amplitude
ampl = 10;

% Initialize the state and forcing
PDE.x{1}.vars = [s1;s2];    PDE.x{1}.dom = [a,b;c,d];
PDE.w{1}.vars = s2;         PDE.w{1}.dom = [c,d];
PDE.w{2}.vars = s2;         PDE.w{2}.dom = [c,d];
PDE.w{3}.vars = s1;         PDE.w{3}.dom = [a,b];
PDE.w{4}.vars = s1;         PDE.w{4}.dom = [a,b];

% % Declare the PDE dynamics
% u_{t}(x,y,t) = u_xx(x,y,t) + u_yy(x,y,t) +u(x,y,t)
PDE.x{1}.term{1}.x = 1;     PDE.x{1}.term{1}.C = [1,1,(2.5*pi)^2];
PDE.x{1}.term{1}.D = [2,0;0,2;0,0];

% % BCs
% u(-1,y,t) = w1(y,t);                  u(1,y,t) = w2(y,t);
PDE.BC{1}.term{1}.x = 1;                PDE.BC{2}.term{1}.x = 1;
PDE.BC{1}.term{1}.loc = [a,s2];         PDE.BC{2}.term{1}.loc = [b,s2];
PDE.BC{1}.term{2}.w = 1;                PDE.BC{2}.term{2}.w = 2;
PDE.BC{1}.term{2}.C = -1;               PDE.BC{2}.term{2}.C = -1;

% u(x,0,t) = w3(x,t);                   u_{y}(x,1,t) = w4(x,t);
PDE.BC{3}.term{1}.x = 1;                PDE.BC{4}.term{1}.x = 1;
PDE.BC{3}.term{1}.loc = [s1,c];         PDE.BC{4}.term{1}.loc = [s1,d];
                                        PDE.BC{4}.term{1}.D = [0,1];
PDE.BC{3}.term{2}.w = 3;                PDE.BC{4}.term{2}.w = 4;
PDE.BC{3}.term{2}.C = -1;               PDE.BC{4}.term{2}.C = -1;


% % Set the forcing
% %   w1(y,t) = ampl*cos(2.5*pi*y);
% %   w2(y,t) = 3*ampl*cos(2.5*pi*y);
% %   w3(x,t) = ampl*(x+2);     
% %   w4(x,t) = -2.5*pi*ampl*(x+2);
uinput.w(1) = ampl*cos(2.5*sym(pi)*sy);
uinput.w(2) = 3*ampl*cos(2.5*sym(pi)*sy);
uinput.w(3) = ampl*(sx+2);
uinput.w(4) = -2.5*sym(pi)*ampl*(sx+2);

% % Set the initial conditions
% %   u(x,y,0) = ampl*(x+2)*cos(2.5*pi*y);
u_ic = ampl*(sx+2)*cos(2.5*sym(pi)*sy);
uinput.ic.x(1) = u_ic;

% % Set the exact solution
% %   u(x,y,t) = ampl*(x+2)*cos(2.5*pi*y);
u_ex = ampl*(sx+2)*cos(2.5*sym(pi)*sy);
uinput.exact(1) = u_ex;



%----------------------------------------
%% Example 23 - Heat Equation with Boundary and Interior Forcing
%----------------------------------------
case 23
% % Heat equation
% %   u_{t}(x,y,t) = u_{xx}(x,y,t) +u_{yy}(x,y,t) +w1(x,y,t)    (x,y) in [-1,1]x[-1,1]
% % with forcing along all boundaries
% %   u(-1,y,t) = w2(y,t);       u(y,1,t) = w3(y,t);
% %   u(x,-1,t) = w4(x,t);       u(x,1,t) = w5(y,t);
% % Using forcing
% %   w1(x,y,t) = ((1.5*pi)^2*t+1)*ampl*y*sin(1.5*pi*x)
% %   w2(y,t) = ampl*y*t;                 w3(y,t) = -ampl*y*t;
% %   w4(x,t) = -ampl*sin(1.5*pi*x)*t;    w5(x,t) = ampl*sin(1.5*pi*x)*t;
% % Starting with an initial condition
% %   u(x,y,0) = 0;
% % we get an exact solution
% %   u(x,y,t) = ampl*y*sin(1.5*pi*x)*t;

% % Set the spatial domain (exact solution is only valid for [-1,1]x[-1,1] !!!)
a=-1;   b=1;
c=-1;   d=1;
uinput.dom=[a b;c d]; 

% Set the amplitude
ampl = 1;

% Initialize the state and forcing
PDE.x{1}.vars = [s1;s2];    PDE.x{1}.dom = [a,b;c,d];
PDE.w{1}.vars = [s1;s2];    PDE.w{1}.dom = [a,b;c,d];
PDE.w{2}.vars = s2;         PDE.w{2}.dom = [c,d];
PDE.w{3}.vars = s2;         PDE.w{3}.dom = [c,d];
PDE.w{4}.vars = s1;         PDE.w{4}.dom = [a,b];
PDE.w{5}.vars = s1;         PDE.w{5}.dom = [a,b];

% % Declare the PDE dynamics
% u_{t}(x,y,t) = u_xx(x,y,t) + u_yy(x,y,t) +w1(x,y,t)
PDE.x{1}.term{1}.x = 1;     PDE.x{1}.term{1}.C = [1,1];
PDE.x{1}.term{1}.D = [2,0;0,2];

PDE.x{1}.term{2}.w = 1;

% % BCs
% u(-1,y,t) = w2(y,t);                  u(1,y,t) = w3(y,t);
PDE.BC{1}.term{1}.x = 1;                PDE.BC{2}.term{1}.x = 1;
PDE.BC{1}.term{1}.loc = [a,s2];         PDE.BC{2}.term{1}.loc = [b,s2];
PDE.BC{1}.term{2}.w = 2;                PDE.BC{2}.term{2}.w = 3;
PDE.BC{1}.term{2}.C = -1;               PDE.BC{2}.term{2}.C = -1;

% u(x,-1,t) = w4(x,t);                   u_{y}(x,1,t) = w5(x,t);
PDE.BC{3}.term{1}.x = 1;                PDE.BC{4}.term{1}.x = 1;
PDE.BC{3}.term{1}.loc = [s1,c];         PDE.BC{4}.term{1}.loc = [s1,d];
                                        PDE.BC{4}.term{1}.D = [0,1];
PDE.BC{3}.term{2}.w = 4;                PDE.BC{4}.term{2}.w = 5;
PDE.BC{3}.term{2}.C = -1;               PDE.BC{4}.term{2}.C = -1;


% % Set the forcing
% %   w1(x,y,t) = ((1.5*pi)^2*t+1)*ampl*y*sin(1.5*pi*x)
% %   w2(y,t) = ampl*y*t;               w3(y,t) = -ampl*y*t;
% %   w4(x,t) = -ampl*sin(1.5*pi*x)*t;  w5(x,t) = ampl*sin(1.5*pi*x)*t;
uinput.w(1) = ((1.5*sym(pi))^2*st+1)*ampl*sy*sin(1.5*sym(pi)*sx);
uinput.w(2) = ampl*sy*st;
uinput.w(3) = -ampl*sy*st;
uinput.w(4) = -ampl*sin(1.5*sym(pi)*sx)*st;
uinput.w(5) = ampl*sin(1.5*sym(pi)*sx)*st;

% % Set the initial conditions
% %   u(x,y,0) = 0;
uinput.ic.x(1) = sym(0);

% % Set the exact solution
% %   u(x,y,t) = ampl*y*sin(1.5*pi*x)*t;
u_ex = ampl*sy*sin(1.5*sym(pi)*sx)*st;
uinput.exact(1) = u_ex;



%----------------------------------------
%% Example 24 - 2D wave equation with Neumann boundary forcing
%----------------------------------------
case 24
% % Standard 2D wave equation
% %   u_tt=C^2*(u_xx + u_yy)                        (x,y) in [0,pi]x[-pi,pi]
% % with forcing at the boundaries
% %   u_{x}(0,y,t)=w(y,t); u(pi,y,t)=0; u(x,-pi,t)=0;  u(x,pi,t)=0;   
% % Using forcing
% %   w(y,t) = 3*pi*sin(3*pi*y)*sin(sqrt(18)*C*t)
% % and starting with initial condition
% %       u(x,y,0) = 0;
% %   u_{t}(x,y,0) = sqrt(18)*C*sin(3*x)*sin(3*y);
% % we obtain exact solution
% %   u(x,y,t) = sin(sqrt(18)*C*t)*sin(3*x)*sin(3*y);

% % Set the spatial domain 
a=0;        b=pi;
c=-pi;      d=pi;
uinput.dom=[a b;c d]; 
% NOTE: Exact solution is valid only on [0,pi]x[-pi,pi];

% Set the wave velocity 
C = 1/sqrt(18); 

% % % Construct the PDE.
% % We represent it in the form
% %     x1_{t} = x2;
% %     x2_{t} = C^2*(x1_{xx} +x1_{yy});
%%% Term-based input format
PDE.x{1}.vars = [s1;s2];   PDE.x{1}.dom = uinput.dom;   % x1 = u;
PDE.x{2}.vars = [s1;s2];   PDE.x{2}.dom = uinput.dom;   % x2 = u_{t}
PDE.x{2}.diff = [2,2];
PDE.w{1}.vars = s2;        PDE.w{1}.dom = [c,d];    % w1 = w
PDE.w{2}.vars = s2;        PDE.w{2}.dom = [c,d];    % w2 = w_{t}

% x1_{t} = x2;
PDE.x{1}.term{1}.x = 2;

% x2_{t} = C^2*(x1_{xx} + x1_{yy});
PDE.x{2}.term{1}.x = 1;
PDE.x{2}.term{1}.D = [2,0;0,2];
PDE.x{2}.term{1}.C = [C^2,C^2];

% % Set the boundary conditions
% x1_{x}(0,y,t) = w1(y,t);              x1(pi,y,t) = 0;
PDE.BC{1}.term{1}.x = 1;                PDE.BC{2}.term{1}.x = 1;
PDE.BC{1}.term{1}.loc = [a,s2];         PDE.BC{2}.term{1}.loc = [b,s2];
PDE.BC{1}.term{1}.D = [1,0];
PDE.BC{1}.term{2}.w = 1;
PDE.BC{1}.term{2}.C = -1;

% x1(x,-pi,t) = 0;                      x1(x,pi,t) = 0;
PDE.BC{3}.term{1}.x = 1;                PDE.BC{4}.term{1}.x = 1;
PDE.BC{3}.term{1}.loc = [s1,c];         PDE.BC{4}.term{1}.loc = [s1,d];

% x2_{x}(0,y,t) = w2(y,t);              x2(pi,y,t) = 0;
PDE.BC{5}.term{1}.x = 2;                PDE.BC{6}.term{1}.x = 2;
PDE.BC{5}.term{1}.loc = [a,s2];         PDE.BC{6}.term{1}.loc = [b,s2];
PDE.BC{5}.term{1}.D = [1,0];
PDE.BC{5}.term{2}.w = 2;
PDE.BC{5}.term{2}.C = -1;

% x2(x,-pi,t) = 0;                      x2(x,pi,t) = 0;
PDE.BC{7}.term{1}.x = 2;                PDE.BC{8}.term{1}.x = 2;
PDE.BC{7}.term{1}.loc = [s1,c];         PDE.BC{8}.term{1}.loc = [s1,d];


% % Set the forcing
% %   w(y,t) = 3*sin(3*pi*y)*sin(sqrt(18)*C*t);
uinput.w(1) = 3*sin(3*sy)*sin(sqrt(18)*C*st);
uinput.w(2) = diff(uinput.w(1),st);

% % Set the initial conditions
% %       u(x,y,0) = 0;
% %   u_{t}(x,y,0) = sqrt(18)*C*sin(3*x)*sin(3*y);
uinput.ic.x(1) = sym(0);
uinput.ic.x(2) = sqrt(18)*C*sin(3*sx)*sin(3*sy);

% % Set the exact solution
% %   u(x,y,t) = sin(sqrt(18)*C*t)*sin(3*x)*sin(3*y);
uinput.exact(1) = sin(sqrt(18)*C*st)*sin(3*sx)*sin(3*sy);
uinput.exact(2) = diff(uinput.exact(1),st);

end 