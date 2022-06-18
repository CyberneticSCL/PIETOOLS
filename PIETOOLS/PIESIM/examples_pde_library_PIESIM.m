
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% examples_pde_library_PIESIM.m     PIETOOLS 2021b
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This file contains a library of examples to use with PIESIM that includes: 
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
% authorship, and a brief description of modifuinput.ic.PDEations
%
% Initial coding YP  - 5_29_2021


function [PDE,uinput]=examples_pde_library_PIESIM(example);
syms st sx;
pvar s;

switch example
    
    % NOTE: uinput.ic.PDE signifies initial conditions
    %       uinput.w, uinput.u signifies inhomogeneous inputs (either through
    %       boundary conditions or disctibuted)
    %       input.uinput.wdot, uinput.udot signifies temporal derivatives of the inputs uinput.w,
    %       uinput.u
    
    % In the following examples: spatial variable is denoted as x
    % Temporal variable is denoted as t
    % Symbolic variabls for x is sx, and a symbolic variable for t is st
    % Domain boundaries are specified as x\in[a,b]
    
%----------------------------------------
% Parabolic equation examples
%----------------------------------------
%----------------------------------------
% Heat equation u_t=visc*u_xx
%----------------------------------------
%----------------------------------------
%   Dirichlet-Dirichlet examples
%----------------------------------------

%----------------------------------------
%Example 1 - Dirichlet-Dirichlet boundary conditions
%----------------------------------------

    case 1

% u(x,t)=sin(pi*x)*exp(-visc*pi^2*t) - solution for u_t=visc*uxx
% Function sin(pi*x); 2nd derivative is ic for p=2. 
% This is for Dirichlet-Dirichlet test

% Boundary conditions: u(a,t)=sin(pi*a)*exp(-visc*pi^2*t), u(b,t)=sin(pi*b)*exp(-visc*pi^2*t)
% Non-zero boundary conditions are treated as disturbances in PIE framework
            
   
             %a=-1;b=1; 
             a=2.2; b=3.3;
             PDE.dom=[a b];  
             PDE.n0=0; PDE.n1=0; PDE.n2=1; PDE.nf=0; PDE.nw=2; PDE.nu=0; PDE.nx=0;  
             visc = 0.5;  
             PDE.A0=0; PDE.A1=0; PDE.A2=visc;
             PDE.B=[1 0 0 0;0 1 0 0];
             PDE.Bw=[1 0;0 1];
             
             % Input initial conditions on the fundamental states:
             % n0+n1+n2 initial conditons
             % 2nd spatial derivative of u (=uxx) is ic for p=2;
             
             uinput.exact =  sin(pi*sx)*exp(-visc*pi^2*st);
             uinput.ic.PDE=  -pi^2*sin(pi*sx);
             
             % Input n1+2*n2 boundary conditions in the order corresponding
             % to the rows of matrix B
             uinput.w(1)= sin(pi*a)*exp(-visc*pi^2*st);
             uinput.w(2)= sin(pi*b)*exp(-visc*pi^2*st);
                   
 
%---------------------------------------------------        
% Example 2 - second example with Dirichlet-Dirichlet boundary conditions
%---------------------------------------------------
 
    case 2
        
% u(x,t)=sin(5*pi/4*x+alpha)*exp(-visc*pi^2*t) - solution for u_t=visc*uxx
% Function sin(5*pi/4*x+alpha); 2nd derivative is ic for p=2. 
% This is for Dirichlet-Dirichlet test

% Boundary conditions: u(a,t)=sin(5*pi/4*a+alpha)*exp(-visc*(5*pi/4)^2*t), u(b,t)=sin(5*pi/4*b+alpha)*exp(-visc*(5*pi/4)^2*t)
% Non-zero boundary conditions are treated as disturbances in PIE framework
% This would give a zero boundary condition on the left, and time-dependent
% % % boundary condition on the right for [-1,1] domain

             a=2.2;b=3.3; PDE.dom=[a b]; 
             PDE.n0=0; PDE.n1=0; PDE.n2=1; PDE.nf=0; 
             PDE.nw=2; PDE.nu=0; PDE.nx=0;
             visc=0.5;  
             PDE.A0=0; PDE.A1=0; PDE.A2=visc;
             PDE.B=[1 0 0 0;0 1 0 0];
             PDE.Bw=[1 0;0 1];
          alpha=pi/8;
          uinput.exact =  sin(5*pi/4*sx+alpha)*exp(-visc*(5*pi/4)^2*st);
% % % % %  2nd spatial derivative of u (=uxx) is ic for p=2;
          uinput.ic.PDE= -(5*pi/4)^2*sin(5*pi/4*sx+alpha);
%          funxxder= @(x,t) -(5*pi/4)^2*sin(5*pi/4*x+alpha)*exp(-visc*(5*pi/4)^2*t);
          uinput.w(1)=sin(5*pi/4*a+alpha)*exp(-visc*(5*pi/4)^2*st);
          uinput.w(2)=sin(5*pi/4*b+alpha)*exp(-visc*(5*pi/4)^2*st);     



%------------------------------------------------------------
% Example 3 - Third example with Dirichlet-Dirichlet boundary conditions 
%-------------------------------------------------------------

    case 3

% u(x,t)=sin(pi*(x-alpha))*exp(-visc*pi^2*t) - solution for u_t=visc*uxx
% Function sin(pi*(x-alpha)); 2nd derivative is ic for p=2. 
% This is for Dirichlet-Dirichlet
% test

% Boundary conditions: u(a,t)=sin(pi*(a-alpha))*exp(-visc*pi^2*t), u(b,t)=sin(pi*(b-alpha))*exp(-visc*pi^2*t)


% This would give time-dependence on both ends for [-1,1] domain (unless alpha=0 == Test case
% 1)

             alpha=pi/8;
             a=2;b=3; PDE.dom=[a b]; 
             PDE.n0=0; PDE.n1=0; PDE.n2=1; PDE.nf=0;  
             PDE.nw=2;
             visc = 0.5;  
             PDE.A0=0; PDE.A1=0; PDE.A2=visc;
             PDE.B=[1 0 0 0;0 1 0 0];
             PDE.Bw=[1 0;0 1];
% %             
           uinput.exact = sin(pi*(sx-alpha))*exp(-visc*pi^2*st);
% % % % %  2nd spatial derivative of u (=uxx) is ic for p=2;
           uinput.ic.PDE= -pi^2*sin(pi*(sx-alpha));
           uinput.w(1)= sin(pi*(a-alpha))*exp(-visc*pi^2*st);
           uinput.w(2)= sin(pi*(b-alpha))*exp(-visc*pi^2*st);
           
%-------------------------------------------------- 
% Example 4 -  Dirichlet-Dirichlet with forcing
%----------------------------------------------------

    case 4

% u(x,t)=sin(pi*x)*t - solution for u_t=visc*uxx+f(x,t), f(x,t)=sin(pi*x)+visc*pi^2*sin(pi*x)*t
% uinput.ic.PDE=0 
% This is for Dirichlet-Dirichlet test

% Boundary conditions: u(a,t)=sin(pi*a)t, u(b,t)=sin(pi*b)*t
% % %             
            a=-1;b=1; PDE.dom=[a b];   
            PDE.n0=0; PDE.n1=0; PDE.n2=1; PDE.nf=2;  
            PDE.nw=2;
            visc = 0.5;  
            PDE.A0=0; PDE.A1=0; PDE.A2=visc;
            PDE.B=[1 0 0 0;0 1 0 0];
            PDE.Bw=[1 0;0 1];
            
            uinput.exact(1) = sin(pi*sx)*st;
            
% % % % % %  2nd spatial derivative of u (=uxx) is ic for p=2;
            uinput.ic.PDE=  0;
            uinput.w(1)= sin(pi*a)*st;
            uinput.w(2)=sin(pi*b)*st;
                                    
%           Temporal dependence of force components in separated form
%           uinput.tforce(1:nf) - same temporal dependence for every state
            uinput.tforce(1)=1+0*st;
            uinput.tforce(2)=st;

%           
%           % Spatial dependence of force components in separated form
%           uinput.force(1:n0+n1+n2,1:nf) - can be different spatial dependence for
%           different states
            uinput.force(1,1) = sin(pi*sx); 
            uinput.force(1,2) = visc*pi^2*sin(pi*sx);

%----------------------------------------------------------------------
%            Example 5 -  Dirichlet-Dirichlet with nonlienar forcing in
%            time
%----------------------------------------------------------------------

    case 5

% u(x,t)=sin(pi*x)*sqrt(t+1) - solution for u_t=visc*uxx+f(x,t),
% f(x,t)=0.5*(t+1)^{-1/2}*sin(pi*x)+visc*pi^2*sin(pi*x)*sqrt(t+1) 
% This is for Dirichlet-Dirichlet test

% Boundary conditions: u(a,t)=sin(pi*a)sqrt(t+1),
% u(b,t)=sin(pi*b)*sqrt(t+1)
% % %             
            a=-1;b=1;  PDE.dom=[a b];  
            PDE.n0=0; PDE.n1=0; PDE.n2=1; PDE.nf=2; 
            PDE.nw=2;
            visc = 0.5;  
            PDE.A0=0; PDE.A1=0; PDE.A2=visc;
            PDE.B=[1 0 0 0;0 1 0 0];
            PDE.Bw=[1 0;0 1];
            
         uinput.exact(1) = sin(pi*sx)*sqrt(st+1);
% % % % % %  2nd spatial derivative of u (=uxx) is ic for p=2;
           uinput.ic.PDE(1)=  -pi^2*sin(pi*sx);
           uinput.w(1)= sin(pi*a)*sqrt(st+1);
           uinput.w(2)= sin(pi*b)*sqrt(st+1);
           
%           % Temporal dependence of force components in separated form
            uinput.tforce(1)=0.5/sqrt(st+1);
            uinput.tforce(2)=sqrt(st+1);
%           
%           % Spatial dependence of force components in separated form
            uinput.force(1,1) =  sin(pi*sx); 
            uinput.force(1,2) =  visc*pi^2*sin(pi*sx);
%            

%----------------------------------------
%   Dirichlet-Neumann examples
%----------------------------------------

%-------------------------------------------------------      
% Example 6 - Dirichlet-Neumann boundary conditions
%-------------------------------------------------------

case 6

% u(x,t)=sin(pi*x)*exp(-visc*pi^2*t) - solution for u_t=visc*uxx
% Function sin(pi*x); 2nd derivative is ic for p=2. 
% This is for Dirichlet-Neumann test (Dirichlet on the
% right, Neumann on the left)
             a=-1;b=1; PDE.dom=[a b];   
             PDE.n0=0; PDE.n1=0; PDE.n2=1; PDE.nf=0; 
             PDE.nw=2;
             visc = 0.5;  
             PDE.A0=0; PDE.A1=0; PDE.A2=visc;
             PDE.B=[1 0 0 0;0 0 0 1];
             PDE.Bw=[1 0;0 1];
 
         uinput.exact(1) = sin(pi*sx)*exp(-visc*pi^2*st);
% % % % %  2nd spatial derivative of u (=uxx) is ic for p=2;
          uinput.ic.PDE(1)= -pi^2*sin(pi*sx);
          uinput.w(1)= sin(pi*a)*exp(-visc*pi^2*st);
          uinput.w(2)= pi*cos(pi*b)*exp(-visc*pi^2*st)*(b-a)/2;


%-------------------------------------------------------      
% Example 7 - inhomogeneous Dirichlet-Neumann boundary conditions
%-------------------------------------------------------

case 7

% u(x,t)=sin(5*pi/4*x+alpha)*exp(-visc*pi^2*t) - solution for u_t=visc*uxx
% Function sin(5*pi/4*x+alpha); 2nd derivative is ic for p=2. 
% This is for Dirichlet-Neumann test (Dirichlet on the
% left, Neumann on the right)
  
          alpha=pi/8; 
          a=1.7;b=1.8; PDE.dom=[a b];   
          PDE.n0=0; PDE.n1=0; PDE.n2=1; PDE.nf=0;   
          PDE.nw=2;
          visc = 0.5;  
          PDE.A0=0; PDE.A1=0; PDE.A2=visc;
          PDE.B=[1 0 0 0;0 0 0 1];
          PDE.Bw=[1 0;0 1];
          uinput.exact(1) = sin(5*pi/4*sx+alpha)*exp(-visc*(5*pi/4)^2*st);
% % % % %  2nd spatial derivative of u (=uxx) is ic for p=2;
          uinput.ic.PDE(1)=  -(5*pi/4)^2*sin(5*pi/4*sx+alpha);
          uinput.w(1)= sin(5*pi/4*a+alpha)*exp(-visc*(5*pi/4)^2*st);
          uinput.w(2)= 5*pi/4*cos(5*pi/4*b+alpha)*exp(-visc*(5*pi/4)^2*st);


%----------------------------------------------------------------------------
%             Example 8 -  Dirichlet-Neumann with nonlienar forcing in
%            time
%----------------------------------------------------------------------------

    case 8

% u(x,t)=sin(pi*x)*sqrt(t+1) - solution for u_t=visc*uxx+f(x),
% f(x)=0.5*(t+1)^{-1/2}*sin(pi*x)+visc*pi^2*sin(pi*x)*sqrt(t+1)
% This is for Dirichlet-Neumann test

% Boundary conditions: u(a,t)=sin(pi*a)t, u(b,t)=pi*cos(pi*b)*t
% % % %             
          a=2;b=3;   PDE.dom=[a b]; 
          PDE.n0=0; PDE.n1=0; PDE.n2=1; PDE.nf=2; 
          PDE.nw=2;
          visc = 0.5;  
          PDE.A0=0; PDE.A1=0; PDE.A2=visc;
          PDE.B=[1 0 0 0;0 0 0 1];
          PDE.Bw=[1 0;0 1];
          uinput.exact(1) = sin(pi*sx)*sqrt(st+1);
% % % % % %  2nd spatial derivative of u (=uxx) is ic for p=2;
           uinput.ic.PDE(1)=  -pi^2*sin(pi*sx);
           uinput.w(1)= sin(pi*a)*sqrt(st+1);
           uinput.w(2)= pi*cos(pi*b)*sqrt(st+1);
           
%           % Temporal dependence of force components in separated form
            uinput.tforce(1)=0.5/sqrt(st+1);
            uinput.tforce(2)=sqrt(st+1);
%           
%           % Spatial dependence of force components in separated form
            uinput.force(1) =  sin(pi*sx); 
            uinput.force(2) =  visc*pi^2*sin(pi*sx);
%            


%-----------------------------------------------------------------------------
%  Example 9 - Convection-diffusion equation with Dirichlet-Dirichlet boundary
%  conditions
%-----------------------------------------------------------------------------

  case 9
        
% u(x,t)=sin(pi*(x-ct))*exp(-visc*pi^2*t) - solution for u_t+c*ux=visc*uxx
% Function sin(pi*x); 2nd derivative is ic for p=2. 
% This is for Dirichlet-Dirichlet test

            a=0;b=1; PDE.dom=[a b];   
            PDE.n0=0; PDE.n1=0; PDE.n2=1; PDE.nf=0; PDE.nw=2;            
            visc = 0.5;  
            c = 2;
            PDE.A0=0; PDE.A1=-c; PDE.A2=visc;
            PDE.B=[1 0 0 0;0 1 0 0];
            PDE.Bw=[1 0;0 1];

         uinput.exact(1) = sin(pi*(sx-c*st))*exp(-visc*pi^2*st);
% % % % %  2nd spatial derivative of u (=uxx) is ic for p=2;
          uinput.ic.PDE(1)= -pi^2*sin(pi*sx);
          uinput.w(1)= sin(pi*(a-c*st))*exp(-visc*pi^2*st);
          uinput.w(2)= sin(pi*(b-c*st))*exp(-visc*pi^2*st);


%--------------------------------------------
%  Parabolic equation with all A0, A1, A2 terms and forcing present
%--------------------------------------------
%--------------------------------------------
% Example 10 - Dirichlet-Dirichlet with linear forcing in time
%--------------------------------------------

    case 10

%  u(x,t)=t sin(pi*x) is a solution for u_t=A0_pde*u + A1_pde*u_x+A2_pde*u_xx+f(x,t), where
%  f(x,t)=sin(pi*x) - A0_pde*t*sin(pi*x)-A1_pde*t*pi*cos(pi*x)+t*A2_pde*pi^2*sin(pi*x) 
         
% Diriclet-Dirichlet boundary conditions

            a=1.25;b=2.5;PDE.dom=[a b];    
            PDE.n0=0; PDE.n1=0; PDE.n2=1; PDE.nf=2; PDE.nw=2;
            PDE.A0=4; PDE.A1=2; PDE.A2=0.5;
            PDE.B=[1 0 0 0;0 1 0 0];
            PDE.Bw=[1 0;0 1];
         
           uinput.exact(1) = st*sin(pi*sx);
% % % % % %  2nd spatial derivative of u (=uxx) is ic for p=2;
             uinput.ic.PDE(1)=  0;
             uinput.w(1)= st*sin(pi*a);
             uinput.w(2)= st*sin(pi*b);
                                   
%           % Temporal dependence of force components in separated form
            uinput.tforce(1)=1+0*st;
            uinput.tforce(2)=st;
%           
%           % Spatial dependence of force components in separated form
            uinput.force(1,1) =  sin(pi*sx); 
            uinput.force(1,2) =  -PDE.A0*sin(pi*sx)-PDE.A1*pi*cos(pi*sx)+PDE.A2*pi^2*sin(pi*sx);
             
%-------------------------------------------
% Example 11 - Neumann (left) - Dirichlet (right) boundary conditions with
% linear forcing in time
%--------------------------------------------

    case 11

%  u(x,t)=t sin(pi*x) is a solution for u_t=A0_pde*u + A1_pde*u_x+A2_pde*u_xx+f(x,t), where
%  f(x,t)=sin(pi*x) - A0_pde*t*sin(pi*x)-A1_pde*t*pi*cos(pi*x)+t*A2_pde*pi^2*sin(pi*x) 
         
% Neumann (left) - Dirichlet (right) boundary conditions
% %          
            a=1.7;b=1.8;  PDE.dom=[a b];  
            PDE.n0=0; PDE.n1=0; PDE.n2=1; PDE.nf=2; PDE.nw=2;
            PDE.A0=4; PDE.A1=2; PDE.A2=0.5;
            PDE.B=[0 0 1 0;0 1 0 0];
            PDE.Bw=[1 0;0 1];
            
           uinput.exact(1) = st*sin(pi*sx);
% % % % % % %  2nd spatial derivative of u (=uxx) is ic for p=2;
             uinput.ic.PDE(1)= 0;
%             funxxder=@(x,t) -t*pi^2*sin(pi*x);
             uinput.w(1)= st*pi*cos(pi*a);
             uinput.w(2)= st*sin(pi*b);
             
              % Temporal dependence of force components in separated form
            uinput.tforce(1)=1+0*st;
            uinput.tforce(2)=st;
%           
%           % Spatial dependence of force components in separated form
            uinput.force(1,1) =  sin(pi*sx); 
            uinput.force(1,2) =  -PDE.A0*sin(pi*sx)-PDE.A1*pi*cos(pi*sx)+PDE.A2*pi^2*sin(pi*sx);
 

%-------------------------------------------
% Example 12 - Neumann (left) -Dirichlet (right) boundary conditions, with
% solution and forcing nonlinear in time
%--------------------------------------------

    case 12
             
%  u(x,t)=sqrt(t+1) sin(pi*x) is a solution for u_t=A0_pde*u + A1_pde*u_x+A2_pde*u_xx+f(x,t), where
%  f(x,t)=1/(2*sqrt(t+1))*sin(pi*x) - A0*sqrt(t+1)*sin(pi*x)-A1*sqrt(t+1)*pi*cos(pi*x)+sqrt(t+1)*A2*pi^2*sin(pi*x) 
         
% Neumann (left) - Dirichlet (right) boundary conditions
            a=1.7;b=1.8; PDE.dom=[a b];   
            PDE.n0=0; PDE.n1=0; PDE.n2=1; PDE.nf=2; PDE.nw=2;
            PDE.A0=4; PDE.A1=2; PDE.A2=0.5;
            PDE.B=[0 0 1 0;0 1 0 0];
            PDE.Bw=[1 0;0 1];
            uinput.exact(1) = (st+1)^0.5*sin(pi*sx);
% % % % % % %  2nd spatial derivative of primary state u (=uxx) is ic for p=2;
             uinput.ic.PDE(1)= @(x) -pi^2*sin(pi*sx);
%             funxxder=@(x,t) -(t+1)^0.5*pi^2*sin(pi*x);
             uinput.w(1)= (st+1)^0.5*pi*cos(pi*a);
             uinput.w(2)= (st+1)^0.5*sin(pi*b);
             
             
            % Temporal dependence of force components in separated form
            uinput.tforce(1)=0.5*(st+1)^(-0.5);
            uinput.tforce(2)=(st+1)^0.5;
%           
%           % Spatial dependence of force components in separated form
            uinput.force(1,1) =  sin(pi*sx); 
            uinput.force(1,2) =  -PDE.A0*sin(pi*sx)-PDE.A1*pi*cos(pi*sx)+PDE.A2*pi^2*sin(pi*sx);

%-------------------------------------------------------------
% Euler-Bernoulli beam equation examples u_tt=-c u_xxxx
%-------------------------------------------------------------
%-------------------------------------------------------------
% Example 13 - exact eigenmodes of a cantilever beam vibration
%-------------------------------------------------------------
             case 13
                 
 % NOTE: don't change the domain limits for this test case. Beta values
 % below correspond to the beam of length L=2 only
                 
            a=0;b=2; PDE.dom=[a b]; 
            PDE.n0=0; PDE.n1=0; PDE.n2=2; PDE.nf=0; PDE.nw=0; 
            L=2;
            c=2;
            PDE.A0=0; PDE.A1=0; PDE.A2=[0 -c;1 0];
            PDE.B=[1 0 0 0 0 0 0 0;0 0 0 1 0 0 0 0; 0 0 0 0 1 0 0 0; 0 0 0 0 0 0 0 1];
% %           
% %             
% % % %       exact solution for EB-beam vibration;
% % %         Primary states are state(1)=u_t; state(2)=u_xx; Fundamental states are state(1)_xx=u_txx, 
%             state(2)_xx=u_xxxx.

%               beta=0.937552034355981;  % First eigenmode
%               beta=2.347045566487088;  % Second eigenmode
%               beta=3.927378719118806;  % Third eigenmode
                beta=5.497770367437734;  % Fourth eigenmode
              omega=beta^2*sqrt(c);
              coeff=(cos(beta*L)+cosh(beta*L))/(sin(beta*L)+sinh(beta*L));
              
%             exact solution is as follows:
             mode = cosh(beta*sx)-cos(beta*sx)+coeff*(sin(beta*sx)-sinh(beta*sx));
             mode_x= beta*(sinh(beta*sx)+sin(beta*sx)+coeff*(cos(beta*sx)-cosh(beta*sx)));
             mode_xx= beta^2*(cosh(beta*sx)+cos(beta*sx)+coeff*(-sin(beta*sx)-sinh(beta*sx)));
             mode_xxx= beta^3*(sinh(beta*sx)-sin(beta*sx)+coeff*(-cos(beta*sx)-cosh(beta*sx)));
             mode_xxxx= beta^4*(cosh(beta*sx)-cos(beta*sx)+coeff*(sin(beta*sx)-sinh(beta*sx)));
%            u = mode*cos(omega*st);
             u_t = -omega*mode*sin(omega*st);
             u_x =  mode_x*cos(omega*st);
             u_xx=  mode_xx*cos(omega*st);
             u_tx =  -omega*mode_x*sin(omega*st);
%  %  %        exact solution for fundamental states u_txx and u_xxxx
             u_txx = -omega*mode_xx*sin(omega*st);
             u_xxxx =  mode_xxxx*cos(omega*st);
% % %           
% %  %        exact solution for primary states u_t(x,t), u_xx(x,t)
             uinput.exact(1) =  u_t;
             uinput.exact(2) =  u_xx;

% % % % % % % Initial conditions on fundamental states are exact solutions
% % % % % % % for fundamental states u_txx and u_xxxx evaluated at time t=0;

               uinput.ic.PDE(1) =  subs(u_txx,st,0);
               uinput.ic.PDE(2) =  subs(u_xxxx,st,0);
% % %           
% %             


%-------------------------------------------------------------
% Forced Euler-Bernoulli beam equation examples u_tt=-c u_xxxx+f(x,t)
%-------------------------------------------------------------
%-------------------------------------------------------------
% Example 14 - exact solution with forcing, solution is polynomial in time
%-------------------------------------------------------------
    case 14

% % %       u(x,t)=x^3*t^2 solution for u_tt=-c u_xxxx+2x^3;
% % Primary states are state(1)=u_t; state(2)=u_xx; Fundamental states are 
%   state(1)_xx=u_txx, state(2)_xx=u_xxxx.
% 
%            u =     x.^3*t^2;
%            u_t =   2*x.^3*t;
%            u_tt=   2*x.^3;
%            u_tx=   6*x.^2*t;
%            u_txx=  12*x.*t;
%            u_x =   3*x.^2*t^2;
%            u_xx =  6*x*t^2; 
%            u_xxx = 6*t^2; 
%            u_xxxx=0;
            
            a=1;b=2; PDE.dom=[a b];   
            PDE.n0=0; PDE.n1=0; PDE.n2=2; PDE.nf=1; PDE.nw=4;
            c=2;
            PDE.A0=0; PDE.A1=0; PDE.A2=[0 -c;1 0];
            PDE.B=[1 0 0 0 0 0 0 0;0 0 0 1 0 0 0 0; 0 0 0 0 1 0 0 0; 0 0 0 0 0 0 0 1];
            PDE.Bw=eye(4);
            
            % exact solution of the PDE
            
            u = @(x,t) x.^3*t^2;   
            ux = @(x,t) 3*x.^2*t^2;
% %           
% %  %        exact solution for primary states u_t(x,t), u_xx(x,t)
             uinput.exact(1) =  2*sx.^3*st;
             uinput.exact(2) =  6*sx*st^2;
% %           
% % % % % % % Initial conditions on fundamental states are exact solutions
% % % % % % % for fundamental states u_txx and u_xxxx evaluated at time t=0;

               uinput.ic.PDE(1)= 0;
               uinput.ic.PDE(2)= 0;
               
               % Boundary conditions: need four boundary conditions on primary states 
% %           
               uinput.w(1)= 2*a^3*st;       % Dirichlet: state(1): u_t at x=a
               uinput.w(2)= 6*b*st^2;       % Dirichlet: state(2): u_xx at x=b
               uinput.w(3)= 6*a^2*st;       % Neumann:   state(1)_x: u_tx at x=a
               uinput.w(4)= 6*st^2;         % Neumann:   state(2)_x: u_xxx at x=b
               

% %           Temporal dependence of force components in separated form
              uinput.tforce(1)= 1;
             
%             Spatial dependence of force components in separated form.
%             Needs to be input for every state for every force component.
%             Index numbder of state is first, index number of force
%             component (1:nf) is second
              uinput.force(1,1) = 2*sx.^3; 
              uinput.force(2,1) = 0;
% %            
%--------------------------------------------
% Example 15 - exact solution with forcing, solution is exponentially decaying in time
%--------------------------------------------
    case 15

% % %       u(x,t)=sin x e^(-t) solution for u_tt=-c u_xxxx+(c+1)sinx e^(-t);
%   Primary states are state(1)=u_t; state(2)=u_xx; Fundamental states are 
%   state(1)_xx=u_txx, state(2)_xx=u_xxxx.

            a=2;b=3;PDE.dom=[a b];    
            PDE.n0=0; PDE.n1=0; PDE.n2=2; PDE.nf=1; PDE.nw=4;
            c=2;
            PDE.A0=0; PDE.A1=0; PDE.A2=[0 -c;1 0];
            PDE.B=[1 0 0 0 0 0 0 0;0 0 0 1 0 0 0 0; 0 0 0 0 1 0 0 0; 0 0 0 0 0 0 0 1];
            PDE.Bw=eye(4);
            
% 
%            u = sin(sx)*exp(-st);
%            ux = cos(sx)*exp(-st);
%            ut = -sin(sx)*exp(-st);
%            utx = -cos(sx)*exp(-st);

% exact solution of the PDE
            
            u = @(x,t) sin(x)*exp(-t);   
            ux = @(x,t) cos(x)*exp(-t);
            
% %  %        exact solution for primary states u_t(x,t), u_xx(x,t)
             uinput.exact(1) =  -sin(sx)*exp(-st);
             uinput.exact(2) =  -sin(sx)*exp(-st);
% %           

% % % % % % % Initial conditions on fundamental states are exact solutions
% % % % % % % for fundamental states u_txx and u_xxxx evaluated at time t=0;

               uinput.ic.PDE(1)= sin(sx);
               uinput.ic.PDE(2)= sin(sx);
               
               
               % Boundary conditions: need four boundary conditions on primary states 
% %           
               uinput.w(1)= -sin(a)*exp(-st);  % Dirichlet: state(1): u_t at x=a
               uinput.w(2)= -sin(b)*exp(-st);  % Dirichlet: state(2): u_xx at x=b
               uinput.w(3)= -cos(a)*exp(-st); % Neumann:   state(1)_x: u_tx at x=a
               uinput.w(4)= -cos(b)*exp(-st);  % Neumann:   state(2)_x: u_xxx at x=b
           
               
% %           Temporal dependence of force components in separated form
              uinput.tforce(1)= exp(-st);
             
%             Spatial dependence of force components in separated form.
%             Needs to be input for every state for every force component.
%             Index numbder of state is first, index number of force
%             component (1:nf) is second
              uinput.force(1,1) = (c+1)*sin(sx); 
              uinput.force(2,1) = 0;
              
%---------------------------------------
% Hyperbolic.PDE equation examples
%----------------------------------------
%---------------------------------------
% Transport equation u_t+c u_x=0
%----------------------------------------
%--------------------------------------------
% Example 16 - traveling wave solution
%--------------------------------------------
    case 16 
% u(x,t)=sin(x-ct) solution for u_t+c u_x=0;
% % % % % 
            c=-2;
            a=1;b=2;   PDE.dom=[a b]; 
            PDE.n0=0; PDE.n1=1; PDE.n2=0; PDE.nf=0; PDE.nw=1;
            PDE.A0=0; PDE.A1=-c; PDE.A2=0;
            if (c>=0)
            xb=a;
            PDE.B=[1 0];
            else
            xb=b;
            PDE.B=[0 1];
            end
            PDE.Bw=1;
                
         uinput.exact(1) = sin(sx-c*st);
% % % % % % %  1st spatial derivative is ic for p=1;
         uinput.ic.PDE(1)=  cos(sx);
         uinput.w(1)= sin(xb-c*st);
% %       
%    
%--------------------------------------------
% Example 17 - Gauss bump propagating solution
%--------------------------------------------
    case 17
% Gauss exponential 1/(sigma*sqrt(2*pi))*exp(-0.5*(x-ct+1)/sigma)^2) solution for u_t+cu_x=0;
% % 
 sigma=0.2;
 mu=0;
 c=-2;
            a=-1;b=1;  PDE.dom=[a b];  
            PDE.n0=0; PDE.n1=1; PDE.n2=0; PDE.nf=0; PDE.nw=1;
            PDE.A0=0; PDE.A1=-c; PDE.A2=0;
            if (c>=0)
            xb=a;
            PDE.B=[1 0];
            else
            xb=b;
            PDE.B=[0 1];
            end
            PDE.Bw=1;
                
% 
  uinput.exact(1) = 1/(sigma*sqrt(2*pi))*exp(-0.5*((sx-c*st-mu)/sigma).^2);
% %  1st derivative is ic for p=1;
  uinput.ic.PDE(1)= 1/(sigma*sqrt(2*pi))*exp(-0.5*((sx-mu)/sigma).^2).*(-(sx-mu)/sigma^2);
% % % %  
% % 

%---------------------------------------
% Wave equation u_tt=c^2 u_xx
%----------------------------------------
%--------------------------------------------
% Example 18 - standing wave solution
% Boundary conditions 
% B=[1 0 0 0;0 0 0 1], u_t(a,t), u_x(b,t)
%--------------------------------------------

    case 18

% u(x,t)=sin(pi*x)*cos(c*pi*t) solution for u_tt=c^2 u_xx;
% u_t is primary first state, u_x is second primary state
% Boundary conditions: u_t(a,t)=-c*pi*sin(pi*a)*sin(c*pi*t),
% u_x(b,t)=pi*cos(pi*b)*cos(c*pi*st);
% % dirside - side on whuinput.ic.PDEh Dirichlet condition is assigned: 1 for left
% boundary, 0 for right boundary
            a=-1;b=1; PDE.dom=[a b];   
            PDE.n0=0; PDE.n1=2; PDE.n2=0; PDE.nf=0; PDE.nw=2;
            c=2;
            PDE.A0=0; PDE.A1=[0 c^2; 1 0]; PDE.A2=0;
            PDE.B=[1 0 0 0;0 0 0 1];
            PDE.Bw=eye(2);

%         u_t and u_x are the primary states 
          uinput.exact(1) = -c*pi*sin(pi*sx)*sin(c*pi*st);
          uinput.exact(2) = pi*cos(pi*sx)*cos(c*pi*st);
          
%       Initial conditions for first fundamental state:  u_tx 
%       Initial conditions for second fundamental state:  u_xx 
          uinput.ic.PDE(1) =  0*sx;
          uinput.ic.PDE(2)=  -pi^2*sin(pi*sx);
          uinput.w(1)= -c*pi*sin(pi*a)*sin(c*pi*st);
          uinput.w(2)= pi*cos(pi*b)*cos(c*pi*st);

%------------------------------------------

% Example 19 - another standing wave solution
% Boundary conditions

% B=[1 0 0 0;0 0 0 1], u_t(a,t), u_x(b,t)
%--------------------------------------------

    case 19

% % u(x,t)=sin(pi*x)*sin(c*pi*t) solution for u_tt=c^2 u_xx;
% % Boundary conditions: ut(a,t)=c*pi*sin(pi*a)*cos(c*pi*t), ux(b,t)=pi*cos(pi*b)*sin(c*pi*t)
% % % % 
            a=-1;b=1;  PDE.dom=[a b];  
            PDE.n0=0; PDE.n1=2; PDE.n2=0; PDE.nf=0; PDE.nw=2;
            c=-2;
            PDE.A0=0; PDE.A1=[0 c^2; 1 0]; PDE.A2=0;
            PDE.B=[1 0 0 0;0 0 0 1];
            PDE.Bw=eye(2);

%         ut and ux are the primary states 
          uinput.exact(1) = c*pi*sin(pi*sx)*cos(c*pi*st);
          uinput.exact(2) = pi*cos(pi*sx)*sin(c*pi*st);

%         Initial conditions for the first fundamental state:  utx 
%         Initial conditions for second fundamental state:  uxx 
          uinput.ic.PDE(1) = c*pi^2*cos(pi*sx);
          uinput.ic.PDE(2)=  0;
          uinput.w(1)= c*pi*sin(pi*a)*cos(c*pi*st);
          uinput.w(2)= pi*cos(pi*b)*sin(c*pi*st);

%-----------------------------------------
% Example 20 - traveling wave solution
% Boundary conditions

% B=[1 0 0 0;0 0 0 1], u_t(a,t), u_x(b,t)
%--------------------------------------------

    case 20

% %       u(x,t)=sin(x-ct) solution for u_tt=c^2 u_xx;
% % % Boundary conditions: u_t(a,t)=-c*cos(a-ct), u_x(b,t)=cos(b-ct)
            a=-1;b=1;   PDE.dom=[a b]; 
            PDE.n0=0; PDE.n1=2; PDE.n2=0; PDE.nf=0; PDE.nw=2;
            c=2;
            PDE.A0=0; PDE.A1=[0 c^2; 1 0]; PDE.A2=0;
            PDE.B=[1 0 0 0;0 0 0 1];
            PDE.Bw=eye(2);

%         u_t and u_x are the primary states 
          uinput.exact(1) = -c*cos(sx-c*st);
          uinput.exact(2) = cos(sx-c*st);
%           
%         Initial conditions for the first fundamental state:  u_tx 
%         Initial conditions for second fundamental state:  u_xx 
          uinput.ic.PDE(1) = c*sin(sx);
          uinput.ic.PDE(2)=  -sin(sx);
          uinput.w(1)= -c*cos(a-c*st);
          uinput.w(2)=  cos(b-c*st);

%-----------------------------------------
% Example 21 - traveling wave solution with characteristuinput.ic.PDE boundary
% conditions at the right end
% Use these boundary conditions with c>0 for enhanced stability 
% Boundary conditions

% B=[1 0 0 0;0 0 1 c], u_t(a,t), u_t(b,t)+cu_x(b,t)=0
%--------------------------------------------

    case 21

% %       u(x,t)=sin(x-ct) solution for u_tt=c^2 u_xx;
% % % Boundary conditions: u_t(a,t)=-c*cos(a-ct), u_t(b,t)+c u_x(b,t)=0;
            a=-1;b=1;  PDE.dom=[a b];  
            PDE.n0=0; PDE.n1=2; PDE.n2=0; PDE.nf=0; PDE.nw=1;
            c=2;
            PDE.A0=0; PDE.A1=[0 c^2; 1 0]; PDE.A2=0;
            PDE.B=[1 0 0 0;0 0 1 c];
            PDE.Bw=[1;0];
% % % % % 
%         u_t and u_x are the primary states 
          uinput.exact(1) = -c*cos(sx-c*st);
          uinput.exact(2) = cos(sx-c*st);
%           
%         Initial conditions for the first fundamental state:  u_tx 
%         Initial conditions for second fundamental state:  u_xx 
          uinput.ic.PDE(1) = c*sin(sx);
          uinput.ic.PDE(2)=  -sin(sx);
          uinput.w(1)= -c*cos(a-c*st);


%-----------------------------------------
% Example 22 - propagating Gauss bump
% Boundary conditions 

% B=[1 0 0 0;0 0 0 1], u_t(a,t), u_x(b,t)
%--------------------------------------------
    case 22

% Propagating Gauss bump 1/(sigma*sqrt(2*pi))*exp(-0.5*(x-ct-mu)/sigma)^2) solution for u_tt=c^2u_xx;
% % % 

            a=-1;b=1;  PDE.dom=[a b];  
            PDE.n0=0; PDE.n1=2; PDE.n2=0; PDE.nf=0; PDE.nw=2;
            c=2;
            PDE.A0=0; PDE.A1=[0 c^2; 1 0]; PDE.A2=0;
            PDE.B=[1 0 0 0;0 0 0 1];
            PDE.Bw=eye(2);
            sigma=0.2;
            mu=0;
% % 
    u = 1/(sigma*sqrt(2*pi))*exp(-0.5*((sx-c*st-mu)/sigma).^2);

%         u_t and u_x are the primary states 
    uinput.exact(1) =c*(sx-c*st-mu)/sigma^2.*u;   
    uinput.exact(2) =-1*(sx-c*st-mu)/sigma^2.*u;

%         Initial conditions for the first fundamental state:  u_tx 
%         Initial conditions for second fundamental state:  u_xx 

    uinput.ic.PDE(1)=(c/sigma^2-c*((sx-mu)/sigma^2).^2).*subs(u,st,0);
    uinput.ic.PDE(2)=(-1/sigma^2+((sx-mu)/sigma^2).^2).*subs(u,st,0);
    
    uinput.w(1)=subs(uinput.exact(1),sx,a);
    uinput.w(2)=subs(uinput.exact(2),sx,b);
    
  %-----------------------------------------
% Example 23 - propagating Gauss bump with characteristuinput.ic.PDEs boundary
% conditions at the right end
% Use these boundary conditions with c>0 for enhanced stability 
% Boundary conditions 
% B=[1 0 0 0;0 0 1 c], u_t(a,t), u_t(b,t)+cu_x(b,t)=0

%--------------------------------------------
    case 23

% Propagating Gauss bump 1/(sigma*sqrt(2*pi))*exp(-0.5*(x-ct-mu)/sigma)^2) solution for u_tt=c^2u_xx;
% % % 

            a=-1;b=1;  PDE.dom=[a b];  
            PDE.n0=0; PDE.n1=2; PDE.n2=0; PDE.nf=0; PDE.nw=1;
            c=2;
            PDE.A0=0; PDE.A1=[0 c^2; 1 0]; PDE.A2=0;
            PDE.B=[1 0 0 0;0 0 1 c];
            PDE.Bw=[1;0];
            sigma=0.2;
            mu=0;
% % 
    u = 1/(sigma*sqrt(2*pi))*exp(-0.5*((sx-c*st-mu)/sigma).^2);

%         u_t and u_x are the primary states 
    uinput.exact(1) =c*(sx-c*st-mu)/sigma^2.*u;   
    uinput.exact(2) =-1*(sx-c*st-mu)/sigma^2.*u;

%         Initial conditions for the first fundamental state:  u_tx 
%         Initial conditions for second fundamental state:  u_xx 

    uinput.ic.PDE(1)=(c/sigma^2-c*((sx-mu)/sigma^2).^2).*subs(u,st,0);
    uinput.ic.PDE(2)=(-1/sigma^2+((sx-mu)/sigma^2).^2).*subs(u,st,0);
    
    uinput.w(1)=subs(uinput.exact(1),sx,a);
  
%-----------------------------------------
%  Example 24 - splitting Gauss bump
%  Boundary conditions
%  B=[1 0 0 0;0 0 0 1], u_t(a,t), u_x(b,t)
%--------------------------------------------
    case 24

% Splitting Gauss bump 1/2(sigma*sqrt(2*pi))*exp(-0.5*(x-ct-mu)/sigma)^2)+1/2(sigma*sqrt(2*pi))*exp(-0.5*(x+ct-mu)/sigma)^2) solution for u_tt=c^2u_xx;

            a=-1;b=1;  PDE.dom=[a b];  
            PDE.n0=0; PDE.n1=2; PDE.n2=0; PDE.nf=0; PDE.nw=2;
            c=-2;
            PDE.A0=0; PDE.A1=[0 c^2; 1 0]; PDE.A2=0;
            PDE.B=[1 0 0 0;0 0 0 1];
            PDE.Bw=eye(2);
            sigma=0.2;
            mu=0;
% % 
    uminus =  0.5/(sigma*sqrt(2*pi))*exp(-0.5*((sx-c*st-mu)/sigma).^2);
    uplus=  0.5/(sigma*sqrt(2*pi))*exp(-0.5*((sx+c*st-mu)/sigma).^2);
    u = uminus+uplus;
    
%         u_t and u_x are the primary states 
    uinput.exact(1) = c*(sx-c*st-mu)/sigma^2.*uminus-c*(sx+c*st-mu)/sigma^2.*uplus;   
    uinput.exact(2) = -(sx-c*st-mu)/sigma^2.*uminus-(sx+c*st-mu)/sigma^2.*uplus;
    
%         Initial conditions for the first fundamental state:  u_tx 
%         Initial conditions for second fundamental state:  u_xx 

    uinput.ic.PDE(1)=0*sx;
    uinput.ic.PDE(2)=(-1/sigma^2+((sx-mu)/sigma^2).^2).*subs(u,st,0);
    uinput.w(1)=subs(uinput.exact(1),sx,a);
    uinput.w(2)=subs(uinput.exact(2),sx,b);    
    

%-----------------------------------------
%  Example 25 - This example is added to document a solution of PDEs in a system form. 
%  NOTE: PDEs in this example are uncoupled and are simply solved together to test a matrix implementation.  
%  This example combines u_t=-lambda u equation, transport equation and heat equation
%  into one system. 
%--------------------------------------------
    case 25
        
        % State 1 u(x,t)=u(x,0)exp(-lambda t) solution for u_t=-lambda u;
        
        % State 2 u(x,t)=sin(x-ct) solution for u_t+c u_x=0;
        
        % State 3 u(x,t)=sin(5*pi/4*x+alpha)*exp(-visc*pi^2*t) - solution for u_t=visc*uxx
        
% For state 3: solutoin function is sin(5*pi/4*x+alpha)
% Dirichlet-Dirichlet boundary conditions

% Boundary conditions: u(a,t)=sin(5*pi/4*a+alpha)*exp(-visc*(5*pi/4)^2*t), u(b,t)=sin(5*pi/4*b+alpha)*exp(-visc*(5*pi/4)^2*t)
% % 
% This would give a zero boundary condition on the left, and time-dependent
% % % boundary condition on the right

            c=2; visc=0.5; lambda=1;
            a=2;b=3;   PDE.dom=[a b]; 
            PDE.n0=1; PDE.n1=1; PDE.n2=1; PDE.nf=0; PDE.nw=3;
            PDE.A0=[-lambda 0 0;0 0 0;0 0 0]; PDE.A1=[0 0;-c 0;0 0]; PDE.A2=[0;0;visc];
            PDE.Bw=eye(3);
           
  % State 1: n0 state. No boundary conditions required on this state
            uinput.ic.PDE(1)=cos(sx);
            uinput.exact(1)=uinput.ic.PDE(1)*exp(-lambda*st);
              
  % State 2: n1 state. One boundary condition required on this state
            if (c>=0)
            xb=a;
            PDE.B(1,:)=[1 0 0 0 0 0];
            else
            xb=b;
            PDE.B(1,:)=[0 1 0 0 0 0];
            end
         uinput.exact(2) = sin(sx-c*st);
% % % % % % %  1st spatial derivative is ic for p=1;
         uinput.ic.PDE(2)=  cos(sx);
         uinput.w(1)= sin(xb-c*st);
% %       

 % State 3: n2 state. Two boundary conditions required on this state
             PDE.B(2:3,:)=[0 0 1 0 0 0;0 0 0 1 0 0];
          alpha=pi/8;
          uinput.exact(3) =  sin(5*pi/4*sx+alpha)*exp(-visc*(5*pi/4)^2*st);
% % % % % %  2nd spatial derivative of u (=uxx) is ic for p=2;
           uinput.ic.PDE(3)= -(5*pi/4)^2*sin(5*pi/4*sx+alpha);
           uinput.w(2)=sin(5*pi/4*a+alpha)*exp(-visc*(5*pi/4)^2*st);
           uinput.w(3)=sin(5*pi/4*b+alpha)*exp(-visc*(5*pi/4)^2*st);
         
%-------------------------------------------------------------------------------------
%  Example 26 - This example combines u_t=-lambda u equation, 
%  transport equation and heat equation
%  into one system and has forcing
%-------------------------------------------------------------------------------------
    case 26
        
        % State 1 u(x,t)=cos(x)cos(t) solution for u_t=-lambda u+lambda cos(x)cos(t)-cos(x)sin(t);
        
        % State 2 u(x,t)=x^10*t^5 solution for u_t+c u_x=5*x^10*t^4+10*c*x^9*t^5;
        
        % State 3 u(x,t)=sin(pi*x)*sqrt(t+1) - solution for
        % u_t=visc*uxx+f(x,t), where f(x,t)=0.5*(t+t0)^{-1/2}*sin(pi*x)+visc*pi^2*sin(pi*x)*sqrt(t+t0)  
        
            c=10; visc=[1 3]; lambda=1;
            a=3;b=4;   PDE.dom=[a b]; 
            PDE.n0=1; PDE.n1=1; PDE.n2=2; PDE.nf=6; PDE.nw=5; PDE.nx=0; PDE.nu=0;
            PDE.A0=[-lambda 0 0 0;0 0 0 0;0 0 0 0;0 0 0 0]; PDE.A1=[0 0 0; -c 0 0;0 0 0;0 0 0]; PDE.A2=[0 0; 0 0; visc(1) 0;0 visc(2)];
            t0=1;
            PDE.Bw=eye(5);
            
  % State 1: n0 state. No boundary conditions required on this state
            uinput.ic.PDE(1)=cos(sx);
            uinput.exact(1)=uinput.ic.PDE(1)*cos(st);
            
  % State 2: n1 state. One boundary condition required on this state  
            
            xb=a;
            PDE.B(1,:)=[1 0 0 0 0 0 0 0 0 0];
           
 
         uinput.exact(2) = sx^10*st^5;
% % % % % % %  1st spatial derivative is ic for p=1;
         uinput.ic.PDE(2)=  0;
         
         uinput.w(1)= xb^10*st^5;
         
% %       

 % State 3: n2 state. Two boundary conditions required on this state
             PDE.B(2:5,:)=[0 0 1 0 0 0 0 0 0 0;0 0 0 0 1 0 0 0 0 0;0 0 0 1 0 0 0 0 0 0;0 0 0 0 0 1 0 0 0 0];
             
           uinput.exact(3) =  sin(pi*sx)*sqrt(st+t0);
           uinput.ic.PDE(3)= -pi^2*sin(pi*sx)*sqrt(t0);
           uinput.exact(4) =  sin(pi*sx)*sqrt(st+t0);
           uinput.ic.PDE(4)= -pi^2*sin(pi*sx)*sqrt(t0);
           
           uinput.w(2)=sin(pi*a)*sqrt(st+t0);
           uinput.w(3)=sin(pi*b)*sqrt(st+t0);
           
           uinput.w(4)=sin(pi*a)*sqrt(st+t0);
           uinput.w(5)=sin(pi*b)*sqrt(st+t0);
                  
           
%           % Temporal dependence of force components in separated form
            uinput.tforce(1)=0.5/sqrt(st+t0);
            uinput.tforce(2)=sqrt(st+t0);
            uinput.tforce(3)=st^4;
            uinput.tforce(4)=st^5;
            uinput.tforce(5)=cos(st);
            uinput.tforce(6)=sin(st);
%           
%           % Spatial dependence of force components in separated form
%           uinput.force(1:n0+n1+n2,1:nf) - can be different spatial dependence for
%           different states
            uinput.force(1,1) =  0*sx; 
            uinput.force(2,1) =  0*sx; 
            uinput.force(3,1) =  sin(pi*sx); 
            uinput.force(4,1) =  sin(pi*sx); 
            uinput.force(1,1) =  0*sx;            
            uinput.force(2,2) =  0*sx;
            uinput.force(3,2) =  visc(1)*pi^2*sin(pi*sx);
            uinput.force(4,2) =  visc(2)*pi^2*sin(pi*sx);
            uinput.force(1,3) =  0*sx;            
            uinput.force(2,3) =  5*sx^10;
            uinput.force(3,3) =  0;
            uinput.force(4,3) =  0;
            uinput.force(1,4) =  0*sx;            
            uinput.force(2,4) =  10*c*sx^9;
            uinput.force(3,4) =  0;
            uinput.force(4,4) =  0;
            uinput.force(1,5) =  lambda*cos(sx);            
            uinput.force(2,5) =  0;
            uinput.force(3,5) =  0;
            uinput.force(4,5) =  0;
            uinput.force(1,6) =  -cos(sx);            
            uinput.force(2,6) =  0;
            uinput.force(3,6) =  0;
            uinput.force(4,6) =  0;
            
%----------------------------------------------------------------------------
% Examples of PDE/ODE coupling
%----------------------------------------------------------------------------

%----------------------------------------------------------------------------
%Example 27 - 2 PDE states + 2 ODE states. Dynamics uncoupled.
% ODE here is not very interesting. Seems like nothing happens to ODE -
% check
%----------------------------------------------------------------------------

    case 27
            
% %       u(x,t)=x^3*t^2 solution for u_tt=-c u_xxxx+2x^3;
% % Primary states are state(1)=u_t; state(2)=u_xx; Fundamental states are 
%   state(1)_xx=u_txx, state(2)_xx=u_xxxx.
% 
%            u =     x.^3*t^2;
%            u_t =   2*x.^3*t;
%            u_tt=   2*x.^3;
%            u_tx=   6*x.^2*t;
%            u_txx=  12*x.*t;
%            u_x =   3*x.^2*t^2;
%            u_xx =  6*x*t^2; 
%            u_xxx = 6*t^2; 
%            u_xxxx=0;

% Forcing is done through B matrix
            
            a=-1;b=1; PDE.dom=[a b]; 
            PDE.n0=0; PDE.n1=0; PDE.n2=2; PDE.nf=0; PDE.nw=5; PDE.nx=2; PDE.nu=0;
            c=2;
            PDE.A0=0; PDE.A1=0; PDE.A2=[0 -c;1 0];
            PDE.B=[1 0 0 0 0 0 0 0;0 0 0 1 0 0 0 0; 0 0 0 0 1 0 0 0; 0 0 0 0 0 0 0 1];
            PDE.Bw=[1 0 0 0 0;0 1 0 0 0; 0 0 1 0 0;0 0 0 1 0];
            
            uinput.ic.ODE=[1 2];
            
            % exact solution of the PDE
            
            u = @(x,t) x.^3*t^2;   
            ux = @(x,t) 3*x.^2*t^2;
% %           
% %  %        exact solution for primary states u_t(x,t), u_xx(x,t)
             uinput.exact(1) =  2*sx.^3*st;
             uinput.exact(2) =  6*sx*st^2;
% %           
% % % % % % % Initial conditions on fundamental states are exact solutions
% % % % % % % for fundamental states u_txx and u_xxxx evaluated at time t=0;

               uinput.ic.PDE(1)= 0;
               uinput.ic.PDE(2)= 0;
               
               % Boundary conditions: need four boundary conditions on primary states 
% %           
               uinput.w(1)= 2*a^3*st;       % Dirichlet: state(1): u_t at x=a
               uinput.w(2)= 6*b*st^2;       % Dirichlet: state(2): u_xx at x=b
               uinput.w(3)= 6*a^2*st;       % Neumann:   state(1)_x: u_tx at x=a
               uinput.w(4)= 6*st^2;         % Neumann:   state(2)_x: u_xxx at x=b
               
               % 5th input is distributed forcing
               uinput.w(5)=1;
               
               
               PDE.B21=s*zeros(PDE.n0+PDE.n1+PDE.n2,PDE.nw);  
               PDE.B21(1,5)=2*s^3;
               
% %           Temporal dependence of force components in separated form
%               uinput.tforce(1)= 1;
%              
% %             Spatial dependence of force components in separated form.
% %             Needs to be input for every state for every force component.
% %             Index numbder of state is first, index number of force
% %             component (1:nf) is second
%               uinput.force(1,1) = 2*sx.^3; 
%               uinput.force(2,1) = 0;


%----------------------------------------------------------------------------
% PDE/ODE with coupled dynamics
%----------------------------------------------------------------------------

%----------------------------------------
%Example 28 - Heat Equation coupled with unstable ODE at the boundary 
%----------------------------------------

    case 28

 PDE.n0=0;PDE.n1=0;PDE.n2=1; PDE.nx =1;
 a=0; b=1; PDE.dom=[a b]; 
%
 PDE.A0=0; PDE.A1=0; PDE.A2=1; PDE.A = 1; PDE.E =0;
%
%%% % x(a)=0, x(b)=0   - Unstable always (this line should remain commented)
 PDE.B = [1 0 0 0; 0 1 0 0]; PDE.Bx = [0;1];
 
 uinput.ic.ODE=2;
 uinput.ic.PDE=1;
 
 
%----------------------------------------
% Example 29 - x_t = -x(t-tau) %stable if tau<=pi/2 
%----------------------------------------
 

  case 29
  PDE.nx = 1; PDE.n0=0; PDE.n1 =1; PDE.n2 =0; 
  PDE.A = 0; PDE.E = 0;
  tau = pi/2 -5e-1; 
  PDE.A0 = 0; PDE.A1 = 1/tau; PDE.A2 = zeros(PDE.n2);
  PDE.E0 = -[1 0]; PDE.Ea = 0; PDE.Eb = 0;
  PDE.B=[0 1];
  PDE.Bx = [1];
  a = -1; b =0; PDE.dom = [-1 0];
  uinput.ic.ODE=[1];
  uinput.ic.PDE=[1];
  
  
%----------------------------------------
% Example 30 - Example 1 from Amritam's paper: stable
%----------------------------------------
    case 30
  % % 
 a=0; b= 1; PDE.dom = [a b];
 PDE.nx = 1; PDE.n0=0; PDE.n1=0; PDE.n2 =1;
 PDE.A = -3; 
 PDE.A2 = 1; 
 PDE.B = [0, 0, 1, 0; 0, 1, 0, 0]; 
 PDE.Bx = [0; 0];
 uinput.ic.ODE=[1];
 uinput.ic.PDE=[1];
 
 
%----------------------------------------
% Example 31 - Example 2 from Amritam's paper: I get unstable, Amritam says
% stable
%----------------------------------------
 
    case 31

 lambda = pi^2-0.01; %stable for lambda<pi^2
 a=0; b= 1; PDE.dom = [a b];
 PDE.nx = 4; PDE.n0=0; PDE.n1=0; PDE.n2 =2;
 PDE.A = [-1.2142, 1.9649, 0.2232, 0.5616;
            -1.8042, -0.7260, -0.3479, 5.4355;
            -0.2898, 0.7381, -1.7606, 0.8294;
            -0.9417, -5.3399, -1.0704, -0.7590]; 
 PDE.E0 = [-1.5368 0;0 0.8871;1.0656 0;1.1882 0]*[zeros(PDE.n2) zeros(PDE.n2) eye(PDE.n2) zeros(PDE.n2)];
 PDE.A2 = eye(PDE.n2); PDE.A0 = lambda*eye(PDE.n2); 
 PDE.B = [eye(PDE.n2), zeros(PDE.n2), zeros(PDE.n2),zeros(PDE.n2); 
           zeros(PDE.n2), eye(PDE.n2), zeros(PDE.n2), zeros(PDE.n2)]; 
 PDE.E = [-2.5575 0 1.0368 0;-1.8067 0.4630 1.3621 0];
 PDE.Bx = zeros(4,4);
 
 uinput.ic.ODE=[1 1 1 1];
 uinput.ic.PDE=[1 1];
 
  
%----------------------------------------
% Example 32 
%----------------------------------------
 
    case 32
 
 p = 1; q = 0; qc = 3;

A = [0, -1/4, -1/5, 1/5, 1/6;
1/2, 1, -4, 9/2, 7/2;
-9/4, -1/2, -14, 23, 16;
-1/5, -1/2, -11/4, 1/10, 5/4;
-4/3, -4/3, -9, 9, 5/2];

B = [-7/2, -3/2, -1/10, 1/2, 1]';
C = [1/10, -1/3, -4, 7/8, 7/8];


 PDE.nw = 0;   PDE.ny = 0;   PDE.nz = 0;   PDE.nx = 5;   PDE.nu = 0;
 PDE.n0 = 0;   PDE.n1 = 0;   PDE.n2 = 1;
 PDE.dom = [0,1];

 PDE.A2 = p;   PDE.A0 = qc-q;
 PDE.B22 = 0;    
 
 PDE.A = A; PDE.E0 = [B, zeros(5,3)];
 
 
 PDE.B = [0 0 1 0;
            0 1 0 0];
 PDE.Bx = [zeros(1,5);C];
 
 uinput.ic.ODE=[1 1 1 1 1];
 uinput.ic.PDE=[1];

%----------------------------------------
% Example 33 - example 2 from  Lhachemi & Prieur 2021 -
% Stable for both Dirichlet-Dirichlet and Dirichlet-Neumann
%----------------------------------------
 
    case 33

p = 1; q = 0; qc = 3;

A = [-1/4, -1/6, 2, 1, 1/12;
-3/2, -3/2, 5, 5, 1/6;
3/2, -4, -15/2, -5, -1/3;
-13/2, 22, 22, -14, -1/2;
1/7, -1/2, -1/2, 1/5, -5/2];

B = [-5/4, 2/3, 1/6, -1/6, 0]';
C = [-2/5, -5/4, 3/2, 1/3, 1/40];
    

 PDE.nw = 0;   PDE.ny = 0;   PDE.nz = 0;   PDE.nx = 5;   PDE.nu = 0;
 PDE.n0 = 0;   PDE.n1 = 0;   PDE.n2 = 1;
 PDE.dom = [0,1];

 PDE.A2 = p;   PDE.A0 = qc-q;
 PDE.B22 = 0;    
 
 PDE.A = A; PDE.E0 = [zeros(5,2),B, zeros(5,1)];
 
 
 PDE.B = [1 0 0 0;
            0 0 0 1];
 PDE.Bx = [zeros(1,5);C];
 
  
 uinput.ic.ODE=[-1 1 -2 2 -1];
 uinput.ic.PDE=5*sx*(1-sx)^2*cos(3*pi*sx);


end % cases
% %            