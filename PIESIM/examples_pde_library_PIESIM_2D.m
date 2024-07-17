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
pde_struct PDE;
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
%PDE.x{1}.term{2}.C = 1;


% % % % 
%           Exact solution, initial conditions and inhomogeneous inputs  
              
               uinput.exact =  ampl*(sx-a)*(sx-b)*(sy-c)*(sy-d);
%              % Initial conditions for the primary states of the PDE
               uinput.ic.PDE= uinput.exact;

       % Enter temporal content of disturbance here
               uinput.w(1)=1;
               %uinput.w(1) = -2*ampl*(sy*(sx-a)*(sx-b)+sx*(sy-c)*(sy-d))-lam*ampl*(sx-a)*(sx-b)*(sy-c)*(sy-d);
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


%----------------------------------------
% Example 7 - 2D wave equation with Dirichlet-Dirichlet BCs
%----------------------------------------
case 7
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
% %     u_{t}(x,t) = sum_{m=0}^{infty} sum_{n=0}^{infty} lam_{mn}*C_{mn}*sin(m*pi/(b-a)*x)*sin(n*pi/(d-c)*y);
% % 
% % where lam_{mn}=c*sqrt((m*pi/(b-a))^2 +(n*pi/(d-c))^2)
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
% Example 8 - 2D wave equation with Dirichlet-Neumann BCs
%----------------------------------------
case 8
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
% %     u_{t}(x,t) = sum_{m=0}^{infty} sum_{n=0}^{infty} lam_{mn}*C_{mn}*cos((m-0.5)*pi/(b-a)*x)*sin((n-0.5)*pi/(d-c)*y);
% % 
% % where lam_{mn}=c*sqrt(((m-0.5)*pi/(b-a))^2 +((n-0.5)*pi/(d-c))^2)
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
% Example 9 - 2D wave equation with Dirichlet-Dirichlet BCs and forcing
%----------------------------------------
case 9
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
% Example 10 - 2D wave equation with Dirichlet-Neumann BCs and forcing
%----------------------------------------
case 10
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

end 