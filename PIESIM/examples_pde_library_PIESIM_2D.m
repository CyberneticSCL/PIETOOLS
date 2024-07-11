%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% examples_pde_library_PIESIM_2D.m     PIETOOLS 2021b
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This file contains a library of 2D examples to use with PIESIM that includes: 
% PDE systems and coupled PDE/ODE systems. 
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
syms st sx sy;
pvar s1 s2 
format long;

switch example
    
    % NOTE: uinput.ic.PDE signifies initial conditions
    %       uinput.w, uinput.u signifies inhomogeneous inputs (either through
    %       boundary conditions or disctibuted)
    
    % In the following examples: spatial variable is denoted as x
    % Temporal variable is denoted as t
    % Symbolic variabls for x is sx, and a symbolic variable for t is st
    % Domain boundaries are specified as x\in[a,b]
    
%----------------------------------------
% Parabolic equation examples
%----------------------------------------
%----------------------------------------
% 2D Heat equation u_t=visc*(u_xx+u_yy)
%----------------------------------------

%----------------------------------------
% Example 1 - 2D Heat Equation with Dirichlet-Dirichlet boundary conditions
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
% Example 2 - 2D Heat Equation with Dirichlet-Dirichlet boundary conditions
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

PDE.x{1}.term{1}.D = [2,0; 0,2];
PDE.x{1}.term{1}.C = [c1*eye(ne), c2*eye(ne)];

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
% Example 3 - 2D Heat Equation with Dirichlet-Dirichlet boundary conditions
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
%Example 4 - 2D Heat Equation with Dirichlet-Neumann boundary conditions
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
% Example 5 - 2D KISS model
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
%  

%----------------------------------------
% Example 6 - 2D equation with different differentiability
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
end 
