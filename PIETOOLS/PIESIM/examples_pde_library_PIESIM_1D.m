
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% examples_pde_library_PIESIM_1D.m     PIETOOLS 2024
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This file contains a library of 1D examples to use with PIESIM that includes: 
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
% DJ 12/19/2024: Change variable s --> s1;
% VP 11/22/2025: added cylindrical coordinates examples


function [PDE,uinput]=examples_pde_library_PIESIM_1D(example)
syms st sx;
pvar s1 t;

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
% Heat equation u_t=visc*u_xx
%----------------------------------------
%----------------------------------------
%   Dirichlet-Dirichlet examples
%----------------------------------------

%----------------------------------------
%Example 1 - Heat Equation with Dirichlet-Dirichlet boundary conditions
%----------------------------------------

    case 1

% u(x,t)=sin(pi*x)*exp(-visc*pi^2*t) - solution for u_t=visc*uxx

% Boundary conditions: u(a,t)=sin(pi*a)*exp(-visc*pi^2*t), u(b,t)=sin(pi*b)*exp(-visc*pi^2*t)
% Non-zero boundary conditions are treated as disturbances in PIE framework
            
   
             a=-1;b=1; 
             PDE.dom=[a b]; 
             visc = 0.5;  

             % Batch format

              PDE.n0=0; PDE.n1=0; PDE.n2=1; PDE.nw=2; PDE.nu=0; PDE.no=0; 
              PDE.A0=0; PDE.A1=0; PDE.A2=visc;
              PDE.B=[1 0 0 0;0 1 0 0];
              PDE.Bw=[1 0;0 1];

% %    %          Terms format
%                   PDE.n.n_pde=[0,0,1];
%                   PDE.n.nw=2;
%                   PDE.n.nv=2;
%                   PDE.PDE.A{1}.Lstate=2;
%                   PDE.PDE.A{1}.Rstate=2;
%                   PDE.PDE.A{1}.D=2;
%                   PDE.PDE.A{1}.I=0;
%                   PDE.PDE.A{1}.coeff=visc;
% % % % % % 
% % % % % %        % BCs: u(x=a)       
%                  PDE.BC.Ebb{1}.coeff = [eye(1);zeros(1)];
%                  PDE.BC.Ebb{1}.Rstate = 2; PDE.BC.Ebb{1}.delta = 0;
% % % % % %  
% % % % % %        % BCs: u(x=b)
%                  PDE.BC.Ebb{2}.coeff = [zeros(1);eye(1)];
%                  PDE.BC.Ebb{2}.Rstate = 2; PDE.BC.Ebb{2}.delta = 1; 
% % % % % 
% % % % %          % BCs: inhomogeneous inputs (connect through ODE interconnected signals)
%                  PDE.BC.Ebv=[1 0;0 1];
%                  PDE.ODE.Dvw=[1 0;0 1];

%           Exact solution, initial conditions and inhomogeneous inputs   
             
              uinput.exact =  sin(pi*sx)*exp(-visc*pi^2*st);
             % Initial conditions for the primary states of the PDE
              uinput.ic.PDE=  sin(pi*sx);

             % Input n1+2*n2 boundary conditions in the order corresponding
             % to the rows of matrix B
             uinput.w(1)= sin(pi*a)*exp(-visc*pi^2*st);
             uinput.w(2)= sin(pi*b)*exp(-visc*pi^2*st);                    
%  
%---------------------------------------------------        
% Example 2 - second example of heat equation with Dirichlet-Dirichlet boundary conditions
%---------------------------------------------------
 
    case 2
        
% u(x,t)=sin(5*pi/4*x+alpha)*exp(-visc*pi^2*t) - solution for u_t=visc*uxx

% Boundary conditions: u(a,t)=sin(5*pi/4*a+alpha)*exp(-visc*(5*pi/4)^2*t), u(b,t)=sin(5*pi/4*b+alpha)*exp(-visc*(5*pi/4)^2*t)
% Non-zero boundary conditions are treated as disturbances in PIE framework
% This would give a zero boundary condition on the left, and time-dependent
% % % boundary condition on the right for [-1,1] domain

             a=0;b=1; PDE.dom=[a b]; 
             visc=1;  

             % Batch format
%              PDE.n0=0; PDE.n1=0; PDE.n2=1; 
%              PDE.nw=2; PDE.nu=0; PDE.no=0;
%              PDE.A0=0; PDE.A1=0; PDE.A2=visc;
%              PDE.B=[1 0 0 0;0 1 0 0];
%              PDE.Bw=[1 0;0 1];

             % Terms format

             PDE.n.n_pde=[0,0,1];
             PDE.n.nw=2;
             PDE.n.nv=2;
             PDE.PDE.A{1}.Lstate=2;
             PDE.PDE.A{1}.Rstate=2;
             PDE.PDE.A{1}.D=2;
             PDE.PDE.A{1}.I=0;
             PDE.PDE.A{1}.coeff=visc;

              % BCs: u(x=a)       
             PDE.BC.Ebb{1}.coeff = [eye(1);zeros(1)];
             PDE.BC.Ebb{1}.Rstate = 2; PDE.BC.Ebb{1}.delta = 0;
% % % % %  
% % % % %    % BCs: u(x=b)
             PDE.BC.Ebb{2}.coeff = [zeros(1);eye(1)];
             PDE.BC.Ebb{2}.Rstate = 2; PDE.BC.Ebb{2}.delta = 1; 
% % % % 
% % % %      % BCs: inhomogeneous inputs (connect through ODE interconnected signals)
             PDE.BC.Ebv=[1 0;0 1];
             PDE.ODE.Dvw=[1 0;0 1];


%           Exact solution, initial conditions and inhomogeneous inputs      
          alpha=pi/8;
          uinput.exact =  sin(5*pi/4*sx+alpha)*exp(-visc*(5*pi/4)^2*st);
% % % % %  2nd spatial derivative of u (=uxx) is ic for p=2;

             % Initial conditions for the primary states of the PDE
          uinput.ic.PDE= sin(5*pi/4*sx+alpha);

          uinput.w(1)=sin(5*pi/4*a+alpha)*exp(-visc*(5*pi/4)^2*st);
          uinput.w(2)=sin(5*pi/4*b+alpha)*exp(-visc*(5*pi/4)^2*st);     



%------------------------------------------------------------
% Example 3 - Third example with Dirichlet-Dirichlet boundary conditions 
%-------------------------------------------------------------

    case 3

% u(x,t)=sin(pi*(x-alpha))*exp(-visc*pi^2*t) - solution for u_t=visc*uxx

% Boundary conditions: u(a,t)=sin(pi*(a-alpha))*exp(-visc*pi^2*t), u(b,t)=sin(pi*(b-alpha))*exp(-visc*pi^2*t)


% This would give time-dependence on both ends for [-1,1] domain (unless alpha=0 == Test case
% 1)

             alpha=pi/8;
             a=2;b=3; PDE.dom=[a b]; 
             visc = 0.5;  

             % Batch format

%              PDE.n0=0; PDE.n1=0; PDE.n2=1; 
%              PDE.nw=2;
%              PDE.A0=0; PDE.A1=0; PDE.A2=visc;
%              PDE.B=[1 0 0 0;0 1 0 0];
%              PDE.Bw=[1 0;0 1];

             % Terms format

             PDE.n.n_pde=[0,0,1];
             PDE.n.nw=2;
             PDE.n.nv=2;
             PDE.PDE.A{1}.Lstate=2;
             PDE.PDE.A{1}.Rstate=2;
             PDE.PDE.A{1}.D=2;
             PDE.PDE.A{1}.I=0;
             PDE.PDE.A{1}.coeff=visc;

              % BCs: u(x=a)       
             PDE.BC.Ebb{1}.coeff = [eye(1);zeros(1)];
             PDE.BC.Ebb{1}.Rstate = 2; PDE.BC.Ebb{1}.delta = 0;
% % % % %  
% % % % %    % BCs: u(x=b)
             PDE.BC.Ebb{2}.coeff = [zeros(1);eye(1)];
             PDE.BC.Ebb{2}.Rstate = 2; PDE.BC.Ebb{2}.delta = 1; 
% % % % 
% % % %      % BCs: inhomogeneous inputs (connect through ODE interconnected signals)
             PDE.BC.Ebv=[1 0;0 1];
             PDE.ODE.Dvw=[1 0;0 1];

%           Exact solution, initial conditions and inhomogeneous inputs  
% %             
           uinput.exact = sin(pi*(sx-alpha))*exp(-visc*pi^2*st);

           % Initial conditions for the primary states of the PDE
           uinput.ic.PDE= sin(pi*(sx-alpha));
           uinput.w(1)= sin(pi*(a-alpha))*exp(-visc*pi^2*st);
           uinput.w(2)= sin(pi*(b-alpha))*exp(-visc*pi^2*st);

%----------------------------------------
% Example 4 - Heat Equation with Dirichlet-Dirichlet and non-constant coefficient
%----------------------------------------

    case 4

% u(x,t)=-2xt-x^2 - solution for u_t=x*uxx

% Boundary conditions: u(a,t)=-2at-a^2, u(b,t)=-2bt-b^2
% Non-zero boundary conditions are treated as disturbances in PIE framework
% Non-constant coefficient which is a polynomial in x should be specified
% throught the polynomial variable s
            
   
             a=0;b=2; 
             PDE.dom=[a b]; 

             % Batch format
%              PDE.n0=0; PDE.n1=0; PDE.n2=1; PDE.nw=2; PDE.nu=0; PDE.no=0; 
%              PDE.A0=0; PDE.A1=0; PDE.A2=s;
%              PDE.B=[1 0 0 0;0 1 0 0];
%              PDE.Bw=[1 0;0 1];

             % Terms format

             PDE.n.n_pde=[0,0,1];
             PDE.n.nw=2;
             PDE.n.nv=2;
             PDE.PDE.A{1}.Lstate=2;
             PDE.PDE.A{1}.Rstate=2;
             PDE.PDE.A{1}.D=2;
             PDE.PDE.A{1}.I=0;
             PDE.PDE.A{1}.coeff=s1;

              % BCs: u(x=a)       
             PDE.BC.Ebb{1}.coeff = [eye(1);zeros(1)];
             PDE.BC.Ebb{1}.Rstate = 2; PDE.BC.Ebb{1}.delta = 0;
% % % % %  
% % % % %    % BCs: u(x=b)
             PDE.BC.Ebb{2}.coeff = [zeros(1);eye(1)];
             PDE.BC.Ebb{2}.Rstate = 2; PDE.BC.Ebb{2}.delta = 1; 
% % % % 
% % % %      % BCs: inhomogeneous inputs (connect through ODE interconnected signals)
             PDE.BC.Ebv=[1 0;0 1];
             PDE.ODE.Dvw=[1 0;0 1];
            
%           Exact solution, initial conditions and inhomogeneous inputs  
             
              uinput.exact =  -2*sx*st-sx^2;
             % Initial conditions for the primary states of the PDE
              uinput.ic.PDE=  -sx^2;

             % Input n1+2*n2 boundary conditions in the order corresponding
             % to the rows of matrix B
             uinput.w(1)= -2*a*st-a^2;
             uinput.w(2)= -2*b*st-b^2;
                   

           
%-------------------------------------------------- 
% Example 5 -  Dirichlet-Dirichlet with forcing
%----------------------------------------------------

    case 5

% u(x,t)=sin(pi*x)*t - solution for u_t=visc*uxx+f(x,t), f(x,t)=sin(pi*x)+visc*pi^2*sin(pi*x)*t
% uinput.ic.PDE=0 
% This is for Dirichlet-Dirichlet test

% Boundary conditions: u(a,t)=sin(pi*a)t, u(b,t)=sin(pi*b)*t
% % %             
            a=2;b=3; PDE.dom=[a b];
            visc = 0.5; 

%             % Batch format
             PDE.n0=0; PDE.n1=0; PDE.n2=1; 
             PDE.nw=4; 
             PDE.A0=0; PDE.A1=0; PDE.A2=visc;
             PDE.B=[1 0 0 0;0 1 0 0];
             PDE.Bw=[1 0 0 0;0 1 0 0]; 

%            Non-polynomial in space forcing needs to be entered through
%            symbolic PDE.B21_nonpol matrix
            PDE.B21_nonpol=sx*zeros(PDE.n0+PDE.n1+PDE.n2,PDE.nw);  
            PDE.B21_nonpol(1,3:4)=[sin(pi*sx) visc*pi^2*sin(pi*sx)];

            % Terms format

%              PDE.n.n_pde=[0,0,1];
%              PDE.n.nw=4;
%              PDE.n.nv=4;
%              PDE.PDE.A{1}.Lstate=2;
%              PDE.PDE.A{1}.Rstate=2;
%              PDE.PDE.A{1}.D=2;
%              PDE.PDE.A{1}.I=0;
%              PDE.PDE.A{1}.coeff=visc;
% 
%              % BCs: u(x=a)       
%              PDE.BC.Ebb{1}.coeff = [eye(1);zeros(1)];
%              PDE.BC.Ebb{1}.Rstate = 2; PDE.BC.Ebb{1}.delta = 0;
% % % % % %  
% % % % % %    % BCs: u(x=b)
%              PDE.BC.Ebb{2}.coeff = [zeros(1);eye(1)];
%              PDE.BC.Ebb{2}.Rstate = 2; PDE.BC.Ebb{2}.delta = 1; 
% % % % % 
% % % % %      % BCs: inhomogeneous inputs (connect through ODE interconnected signals)
%              PDE.BC.Ebv=[1 0 0 0;0 1 0 0];
%              PDE.ODE.Dvw=eye(4);
% 
%             % Non-polynomial in space forcing needs to be entered through
%             % symbolic PDE.PDE.Bpw_nonpol matrix 
%             PDE.PDE.Bpw_nonpol=sx*zeros(sum(PDE.n.n_pde),PDE.n.nw);  
%             PDE.PDE.Bpw_nonpol(1,3:4)=[sin(pi*sx) visc*pi^2*sin(pi*sx)];
  

%           Exact solution, initial conditions and inhomogeneous inputs  
            
            uinput.exact(1) = sin(pi*sx)*st;
    
           % Initial conditions for the primary states of the PDE
            uinput.ic.PDE=  0;

            % Disturbances - both boundary disturbances (1,2) and forcing
            % (3,4)
            uinput.w(1)= sin(pi*a)*st;
            uinput.w(2)=sin(pi*b)*st;
            uinput.w(3)=1;
            uinput.w(4)=st;
         
%----------------------------------------------------------------------
%            Example 6 -  Dirichlet-Dirichlet with nonlinear forcing in
%            time
%----------------------------------------------------------------------

    case 6

% u(x,t)=sin(pi*x)*sqrt(t+1) - solution for u_t=visc*uxx+f(x,t),
% f(x,t)=0.5*(t+1)^{-1/2}*sin(pi*x)+visc*pi^2*sin(pi*x)*sqrt(t+1) 
% This is for Dirichlet-Dirichlet test

% Boundary conditions: u(a,t)=sin(pi*a)sqrt(t+1),
% u(b,t)=sin(pi*b)*sqrt(t+1)
% % %             
            a=-1;b=1;  PDE.dom=[a b];  
            visc = 0.5;  

            % Batch format
%             PDE.n0=0; PDE.n1=0; PDE.n2=1; 
%             PDE.nw=4;
%             PDE.A0=0; PDE.A1=0; PDE.A2=visc;
%             PDE.B=[1 0 0 0;0 1 0 0];
%             PDE.Bw=[1 0 0 0;0 1 0 0];
% 
%             % Non-polynomial in space forcing needs to be entered through
%             % symbolic PDE.B21_nonpol matrix
%             PDE.B21_nonpol=sx*zeros(PDE.n0+PDE.n1+PDE.n2,PDE.nw);  
%             PDE.B21_nonpol(1,3:4)=[sin(pi*sx) visc*pi^2*sin(pi*sx)];

            % Terms format

             PDE.n.n_pde=[0,0,1];
             PDE.n.nw=4;
             PDE.n.nv=4;
             PDE.PDE.A{1}.Lstate=2;
             PDE.PDE.A{1}.Rstate=2;
             PDE.PDE.A{1}.D=2;
             PDE.PDE.A{1}.I=0;
             PDE.PDE.A{1}.coeff=visc;

             % BCs: u(x=a)       
             PDE.BC.Ebb{1}.coeff = [eye(1);zeros(1)];
             PDE.BC.Ebb{1}.Rstate = 2; PDE.BC.Ebb{1}.delta = 0;
% % % % %  
% % % % %    % BCs: u(x=b)
             PDE.BC.Ebb{2}.coeff = [zeros(1);eye(1)];
             PDE.BC.Ebb{2}.Rstate = 2; PDE.BC.Ebb{2}.delta = 1; 
% % % % 
% % % %      % BCs: inhomogeneous inputs (connect through ODE interconnected signals)
             PDE.BC.Ebv=[1 0 0 0;0 1 0 0];
             PDE.ODE.Dvw=eye(4);

            % Non-polynomial in space forcing needs to be entered through
            % symbolic PDE.PDE.Bpw_nonpol matrix 
            PDE.PDE.Bpw_nonpol=sx*zeros(sum(PDE.n.n_pde),PDE.n.nw);  
            PDE.PDE.Bpw_nonpol(1,3:4)=[sin(pi*sx) visc*pi^2*sin(pi*sx)];


%           Exact solution, initial conditions and inhomogeneous inputs  
            
         uinput.exact(1) = sin(pi*sx)*sqrt(st+1);
           % Initial conditions for the primary states of the PDE
           uinput.ic.PDE(1)=  sin(pi*sx);
           uinput.w(1)= sin(pi*a)*sqrt(st+1);
           uinput.w(2)= sin(pi*b)*sqrt(st+1);
           uinput.w(3)=0.5/sqrt(st+1);
           uinput.w(4)=sqrt(st+1);
           
%            

%----------------------------------------
%   Dirichlet-Neumann examples
%----------------------------------------

%-------------------------------------------------------      
% Example 7 - Dirichlet-Neumann boundary conditions
%-------------------------------------------------------

case 7

% u(x,t)=sin(pi*x)*exp(-visc*pi^2*t) - solution for u_t=visc*uxx
% Function sin(pi*x); 2nd derivative is ic for p=2. 
% This is for Dirichlet-Neumann test (Dirichlet on the
% right, Neumann on the left)
             a=-1;b=2; PDE.dom=[a b]; 
             visc = 0.5; 

             % Batch format
%              PDE.n0=0; PDE.n1=0; PDE.n2=1; 
%              PDE.nw=2; 
%              PDE.A0=0; PDE.A1=0; PDE.A2=visc;
%              PDE.B=[1 0 0 0;0 0 0 1];
%              PDE.Bw=[1 0;0 1];

             % Terms format

             PDE.n.n_pde=[0,0,1];
             PDE.n.nw=2;
             PDE.n.nv=2;
             PDE.PDE.A{1}.Lstate=2;
             PDE.PDE.A{1}.Rstate=2;
             PDE.PDE.A{1}.D=2;
             PDE.PDE.A{1}.I=0;
             PDE.PDE.A{1}.coeff=visc;

              % BCs: u(x=a)       
             PDE.BC.Ebb{1}.coeff = [eye(1);zeros(1)];
             PDE.BC.Ebb{1}.Rstate = 2; PDE.BC.Ebb{1}.delta = 0;
% % % % %  
% % % % %    % BCs: u_x(x=b)
             PDE.BC.Ebb{2}.coeff = [zeros(1);eye(1)];
             PDE.BC.Ebb{2}.Rstate = 2; PDE.BC.Ebb{2}.delta = 1; 
             PDE.BC.Ebb{2}.D = 1;
% % % % 
% % % %      % BCs: inhomogeneous inputs (connect through ODE interconnected signals)
             PDE.BC.Ebv=[1 0;0 1];
             PDE.ODE.Dvw=[1 0;0 1];

%           Exact solution, initial conditions and inhomogeneous inputs  
 
         uinput.exact(1) = sin(pi*sx)*exp(-visc*pi^2*st);
           % Initial conditions for the primary states of the PDE
          uinput.ic.PDE(1)= sin(pi*sx);
          uinput.w(1)= sin(pi*a)*exp(-visc*pi^2*st);
          uinput.w(2)= pi*cos(pi*b)*exp(-visc*pi^2*st);


%-------------------------------------------------------      
% Example 8 - inhomogeneous Dirichlet-Neumann boundary conditions
%-------------------------------------------------------

case 8

% u(x,t)=sin(5*pi/4*x+alpha)*exp(-visc*pi^2*t) - solution for u_t=visc*uxx
% Function sin(5*pi/4*x+alpha); 
% This is for Dirichlet-Neumann test (Dirichlet on the
% left, Neumann on the right)
  
          alpha=pi/8;
          a=-1;b=1; PDE.dom=[a b]; 
          visc = 0.5;    

          % Batch format
           PDE.n0=0; PDE.n1=0; PDE.n2=1; 
           PDE.nw=2; 
           PDE.A0=0; PDE.A1=0; PDE.A2=visc;
           PDE.B=[1 0 0 0;0 0 0 1];
           PDE.Bw=[1 0;0 1];

%           % Terms format
%              PDE.n.n_pde=[0,0,1];
%              PDE.n.nw=2;
%              PDE.n.nv=2;
%              PDE.PDE.A{1}.Lstate=2;
%              PDE.PDE.A{1}.Rstate=2;
%              PDE.PDE.A{1}.D=2;
%              PDE.PDE.A{1}.I=0;
%              PDE.PDE.A{1}.coeff=visc;
% 
%               % BCs: u(x=a)       
%              PDE.BC.Ebb{1}.coeff = [eye(1);zeros(1)];
%              PDE.BC.Ebb{1}.Rstate = 2; PDE.BC.Ebb{1}.delta = 0;
% % % % % % %  
% % % % % %    % BCs: u_x(x=b)
%              PDE.BC.Ebb{2}.coeff = [zeros(1);eye(1)];
%              PDE.BC.Ebb{2}.Rstate = 2; PDE.BC.Ebb{2}.delta = 1; 
%              PDE.BC.Ebb{2}.D = 1;
% % % % % % 
% % % % %      % BCs: inhomogeneous inputs (connect through ODE interconnected signals)
%              PDE.BC.Ebv=[1 0;0 1];
%              PDE.ODE.Dvw=[1 0;0 1];

%           Exact solution, initial conditions and inhomogeneous inputs  

          uinput.exact(1) = sin(5*pi/4*sx+alpha)*exp(-visc*(5*pi/4)^2*st);
           % Initial conditions for the primary states of the PDE
          uinput.ic.PDE(1)=  sin(5*pi/4*sx+alpha);
          uinput.w(1)= sin(5*pi/4*a+alpha)*exp(-visc*(5*pi/4)^2*st);
          uinput.w(2)= 5*pi/4*cos(5*pi/4*b+alpha)*exp(-visc*(5*pi/4)^2*st);


%----------------------------------------------------------------------------
%             Example 9 -  Dirichlet-Neumann with nonlienar forcing in
%            time
%----------------------------------------------------------------------------

    case 9

% u(x,t)=sin(pi*x)*sqrt(t+1) - solution for u_t=visc*uxx+f(x),
% f(x)=0.5*(t+1)^{-1/2}*sin(pi*x)+visc*pi^2*sin(pi*x)*sqrt(t+1)
% This is for Dirichlet-Neumann test

% Boundary conditions: u(a,t)=sin(pi*a)t, u_x(b,t)=pi*cos(pi*b)*t
% % % %             
          a=1;b=3;   PDE.dom=[a b]; 
          visc = 0.5; 

          % Batch format
% 
%           PDE.n0=0; PDE.n1=0; PDE.n2=1; 
%           PDE.nw=4; 
%           PDE.A0=0; PDE.A1=0; PDE.A2=visc;
%           PDE.B=[1 0 0 0;0 0 0 1];
%           PDE.Bw=[1 0 0 0;0 1 0 0];
% 
%           % Non-polynomial in space forcing needs to be entered through
%             % symbolic PDE.B21_nonpol matrix
%             PDE.B21_nonpol=sx*zeros(PDE.n0+PDE.n1+PDE.n2,PDE.nw);  
%             PDE.B21_nonpol(1,3:4)=[sin(pi*sx) visc*pi^2*sin(pi*sx)];

          % Terms format

           PDE.n.n_pde=[0,0,1];
             PDE.n.nw=4;
             PDE.n.nv=4;
             PDE.PDE.A{1}.Lstate=2;
             PDE.PDE.A{1}.Rstate=2;
             PDE.PDE.A{1}.D=2;
             PDE.PDE.A{1}.I=0;
             PDE.PDE.A{1}.coeff=visc;

             % BCs: u(x=a)       
             PDE.BC.Ebb{1}.coeff = [eye(1);zeros(1)];
             PDE.BC.Ebb{1}.Rstate = 2; PDE.BC.Ebb{1}.delta = 0;
% % % % %  
% % % % %    % BCs: u_x(x=b)
             PDE.BC.Ebb{2}.coeff = [zeros(1);eye(1)];
             PDE.BC.Ebb{2}.Rstate = 2; PDE.BC.Ebb{2}.delta = 1; 
             PDE.BC.Ebb{2}.D = 1;

% % % % 
% % % %      % BCs: inhomogeneous inputs (connect through ODE interconnected signals)
             PDE.BC.Ebv=[1 0 0 0;0 1 0 0];
             PDE.ODE.Dvw=eye(4);

            % Non-polynomial in space forcing needs to be entered through
            % symbolic PDE.PDE.Bpw_nonpol matrix 
            PDE.PDE.Bpw_nonpol=sx*zeros(sum(PDE.n.n_pde),PDE.n.nw);  
            PDE.PDE.Bpw_nonpol(1,3:4)=[sin(pi*sx) visc*pi^2*sin(pi*sx)];

%           Exact solution, initial conditions and inhomogeneous inputs  

           % Initial conditions for the primary states of the PDE
           uinput.ic.PDE(1)=  sin(pi*sx);
           uinput.w(1)= sin(pi*a)*sqrt(st+1);
           uinput.w(2)= pi*cos(pi*b)*sqrt(st+1);
           uinput.w(3)=0.5/sqrt(st+1);
           uinput.w(4)=sqrt(st+1);

          % Exact solution (optional)
          uinput.exact(1) = sin(pi*sx)*sqrt(st+1);
          
%-----------------------------------------------------------------------------
%  Example 10 - Convection-diffusion equation with Dirichlet-Dirichlet boundary
%  conditions
%-----------------------------------------------------------------------------

  case 10
        
% u(x,t)=sin(pi*(x-ct))*exp(-visc*pi^2*t) - solution for u_t+c*ux=visc*uxx
% Function sin(pi*x); 
% This is for Dirichlet-Dirichlet test

            a=0;b=1; PDE.dom=[a b];          
            visc = 0.5; c = 4;   

%             % Batch format
%                PDE.n0=0; PDE.n1=0; PDE.n2=1; PDE.nw=2; 
%                PDE.A0=0; PDE.A1=-c; PDE.A2=visc;
%                PDE.B=[1 0 0 0;0 1 0 0];
%                PDE.Bw=[1 0;0 1];

            % Terms format
% 
                  PDE.n.n_pde=[0,0,1];
                  PDE.n.nw=2;
                  PDE.n.nv=2;

                  PDE.PDE.A{1}.Lstate=2;
                  PDE.PDE.A{1}.Rstate=2;
                  PDE.PDE.A{1}.D=2;
                  PDE.PDE.A{1}.I=0;
                  PDE.PDE.A{1}.coeff=visc;

                  PDE.PDE.A{2}.Lstate=2;
                  PDE.PDE.A{2}.Rstate=2;
                  PDE.PDE.A{2}.D=1;
                  PDE.PDE.A{2}.I=0;
                  PDE.PDE.A{2}.coeff=-c;
% % % % % 
% % % % %        % BCs: u(x=a)       
                 PDE.BC.Ebb{1}.coeff = [eye(1);zeros(1)];
                 PDE.BC.Ebb{1}.Rstate = 2; PDE.BC.Ebb{1}.delta = 0;
% % % % %  
% % % % %        % BCs: u(x=b)
                 PDE.BC.Ebb{2}.coeff = [zeros(1);eye(1)];
                 PDE.BC.Ebb{2}.Rstate = 2; PDE.BC.Ebb{2}.delta = 1; 
% % % % 
% % % % %          % BCs: inhomogeneous inputs (connect through ODE interconnected signals)
                  PDE.BC.Ebv=[1 0;0 1];
                  PDE.ODE.Dvw=[1 0;0 1];

%           Exact solution, initial conditions and inhomogeneous inputs  
           % Initial conditions for the primary states of the PDE
          uinput.ic.PDE(1)= sin(pi*sx);
          uinput.w(1)= sin(pi*(a-c*st))*exp(-visc*pi^2*st);
          uinput.w(2)= sin(pi*(b-c*st))*exp(-visc*pi^2*st);

          % Exact solution (optional)
          uinput.exact(1) = sin(pi*(sx-c*st))*exp(-visc*pi^2*st);


%--------------------------------------------
%  Parabolic equation with all A0, A1, A2 terms and forcing present
%--------------------------------------------
%--------------------------------------------
% Example 11 - Dirichlet-Dirichlet with linear forcing in time
%--------------------------------------------

    case 11

%  u(x,t)=t sin(pi*x) is a solution for u_t=A0_pde*u + A1_pde*u_x+A2_pde*u_xx+f(x,t), where
%  f(x,t)=sin(pi*x) - A0_pde*t*sin(pi*x)-A1_pde*t*pi*cos(pi*x)+t*A2_pde*pi^2*sin(pi*x) 
         
% Diriclet-Dirichlet boundary conditions

            a=1.25;b=2.5;PDE.dom=[a b]; 
        
%           % Batch format
%             PDE.n0=0; PDE.n1=0; PDE.n2=1; 
%             PDE.nw=4;
%             PDE.A0=4; PDE.A1=2; PDE.A2=0.5;
%             PDE.B=[1 0 0 0;0 1 0 0];
%             PDE.Bw=[1 0 0 0;0 1 0 0];
% 
%           % Non-polynomial in space forcing needs to be entered through
%             % symbolic PDE.B21_nonpol matrix
%             PDE.B21_nonpol=sx*zeros(PDE.n0+PDE.n1+PDE.n2,PDE.nw);  
%             PDE.B21_nonpol(1,3:4)=[sin(pi*sx) -PDE.A0*sin(pi*sx)-PDE.A1*pi*cos(pi*sx)+PDE.A2*pi^2*sin(pi*sx)];
%          

            % Terms format
                  PDE.n.n_pde=[0,0,1];
                  PDE.n.nw=4;
                  PDE.n.nv=4;

                  PDE.PDE.A{1}.Lstate=2;
                  PDE.PDE.A{1}.Rstate=2;
                  PDE.PDE.A{1}.D=2;
                  PDE.PDE.A{1}.I=0;
                  PDE.PDE.A{1}.coeff=0.5;

                  PDE.PDE.A{2}.Lstate=2;
                  PDE.PDE.A{2}.Rstate=2;
                  PDE.PDE.A{2}.D=1;
                  PDE.PDE.A{2}.I=0;
                  PDE.PDE.A{2}.coeff=2;

                  PDE.PDE.A{3}.Lstate=2;
                  PDE.PDE.A{3}.Rstate=2;
                  PDE.PDE.A{3}.I=0;
                  PDE.PDE.A{3}.coeff=4;
% % % % % 
% % % % %        % BCs: u(x=a)       
                 PDE.BC.Ebb{1}.coeff = [eye(1);zeros(1)];
                 PDE.BC.Ebb{1}.Rstate = 2; PDE.BC.Ebb{1}.delta = 0;
% % % % %  
% % % % %        % BCs: u(x=b)
                 PDE.BC.Ebb{2}.coeff = [zeros(1);eye(1)];
                 PDE.BC.Ebb{2}.Rstate = 2; PDE.BC.Ebb{2}.delta = 1; 
% % % % 
% % % % %          % BCs: inhomogeneous inputs (connect through ODE interconnected signals)
                 PDE.BC.Ebv=[1 0 0 0;0 1 0 0];
                 PDE.ODE.Dvw=eye(4);

            % Non-polynomial in space forcing needs to be entered through
            % symbolic PDE.PDE.Bpw_nonpol matrix 
            PDE.PDE.Bpw_nonpol=sx*zeros(sum(PDE.n.n_pde),PDE.n.nw);  
            PDE.PDE.Bpw_nonpol(1,3:4)=[sin(pi*sx) -PDE.PDE.A{3}.coeff*sin(pi*sx)-PDE.PDE.A{2}.coeff*pi*cos(pi*sx)+PDE.PDE.A{1}.coeff*pi^2*sin(pi*sx)];


%           Exact solution, initial conditions and inhomogeneous inputs  
           % Initial conditions for the primary states of the PDE
             uinput.ic.PDE(1)=  0;
             uinput.w(1)= st*sin(pi*a);
             uinput.w(2)= st*sin(pi*b);
             uinput.w(3)=1+0*st;
             uinput.w(4)=st;
                                   
            % Exact solution (optional)
           uinput.exact(1) = st*sin(pi*sx);
             
%-------------------------------------------
% Example 12 - Neumann (left) - Dirichlet (right) boundary conditions with
% linear forcing in time
%--------------------------------------------

    case 12

%  u(x,t)=t sin(pi*x) is a solution for u_t=A0_pde*u + A1_pde*u_x+A2_pde*u_xx+f(x,t), where
%  f(x,t)=sin(pi*x) - A0_pde*t*sin(pi*x)-A1_pde*t*pi*cos(pi*x)+t*A2_pde*pi^2*sin(pi*x) 
         
% Neumann (left) - Dirichlet (right) boundary conditions
% %          
            a=1.7;b=1.8;  PDE.dom=[a b]; 

            % Batch format
%              PDE.n0=0; PDE.n1=0; PDE.n2=1; 
%              PDE.nw=4;
%              PDE.A0=4; PDE.A1=2; PDE.A2=0.5;
%              PDE.B=[0 0 1 0;0 1 0 0];
%              PDE.Bw=[1 0 0 0;0 1 0 0];
%  
%  
%              % Non-polynomial in space forcing needs to be entered through
%              % symbolic PDE.B21_nonpol matrix
%              PDE.B21_nonpol=sx*zeros(PDE.n0+PDE.n1+PDE.n2,PDE.nw);  
%              PDE.B21_nonpol(1,3:4)=[sin(pi*sx) -PDE.A0*sin(pi*sx)-PDE.A1*pi*cos(pi*sx)+PDE.A2*pi^2*sin(pi*sx)];
% %             
              % Terms format
                  PDE.n.n_pde=[0,0,1];
                  PDE.n.nw=4;
                  PDE.n.nv=4;

                  PDE.PDE.A{1}.Lstate=2;
                  PDE.PDE.A{1}.Rstate=2;
                  PDE.PDE.A{1}.D=2;
                  PDE.PDE.A{1}.I=0;
                  PDE.PDE.A{1}.coeff=0.5;

                  PDE.PDE.A{2}.Lstate=2;
                  PDE.PDE.A{2}.Rstate=2;
                  PDE.PDE.A{2}.D=1;
                  PDE.PDE.A{2}.I=0;
                  PDE.PDE.A{2}.coeff=2;

                  PDE.PDE.A{3}.Lstate=2;
                  PDE.PDE.A{3}.Rstate=2;
                  PDE.PDE.A{3}.I=0;
                  PDE.PDE.A{3}.coeff=4;
% % % % % 
% % % % %        % BCs: u_x(x=a)       
                 PDE.BC.Ebb{1}.coeff = [eye(1);zeros(1)];
                 PDE.BC.Ebb{1}.Rstate = 2; PDE.BC.Ebb{1}.delta = 0;
                 PDE.BC.Ebb{1}.D = 1;
% % % % %  
% % % % %        % BCs: u(x=b)
                 PDE.BC.Ebb{2}.coeff = [zeros(1);eye(1)];
                 PDE.BC.Ebb{2}.Rstate = 2; PDE.BC.Ebb{2}.delta = 1; 
% % % % 
% % % % %          % BCs: inhomogeneous inputs (connect through ODE interconnected signals)
                 PDE.BC.Ebv=[1 0 0 0;0 1 0 0];
                 PDE.ODE.Dvw=eye(4);

            % Non-polynomial in space forcing needs to be entered through
            % symbolic PDE.PDE.Bpw_nonpol matrix 
            PDE.PDE.Bpw_nonpol=sx*zeros(sum(PDE.n.n_pde),PDE.n.nw);  
            PDE.PDE.Bpw_nonpol(1,3:4)=[sin(pi*sx) -PDE.PDE.A{3}.coeff*sin(pi*sx)-PDE.PDE.A{2}.coeff*pi*cos(pi*sx)+PDE.PDE.A{1}.coeff*pi^2*sin(pi*sx)];


%           Exact solution, initial conditions and inhomogeneous inputs  
           % Initial conditions for the primary states of the PDE
             uinput.ic.PDE(1)= 0;
             uinput.w(1)= st*pi*cos(pi*a);
             uinput.w(2)= st*sin(pi*b);
             uinput.w(3)=1+0*st;
             uinput.w(4)=st;
 
            % Exact solution (optional)
           uinput.exact(1) = st*sin(pi*sx);

%-------------------------------------------
% Example 13 - Neumann (left) -Dirichlet (right) boundary conditions, with
% solution and forcing nonlinear in time
%--------------------------------------------

    case 13
             
%  u(x,t)=sqrt(t+1) sin(pi*x) is a solution for u_t=A0_pde*u + A1_pde*u_x+A2_pde*u_xx+f(x,t), where
%  f(x,t)=1/(2*sqrt(t+1))*sin(pi*x) - A0*sqrt(t+1)*sin(pi*x)-A1*sqrt(t+1)*pi*cos(pi*x)+sqrt(t+1)*A2*pi^2*sin(pi*x) 
         
% Neumann (left) - Dirichlet (right) boundary conditions
            a=1.7;b=1.8; PDE.dom=[a b];   

            % Batch format
%             PDE.n0=0; PDE.n1=0; PDE.n2=1; 
%             PDE.nw=4;
%             PDE.A0=4; PDE.A1=2; PDE.A2=0.5;
%             PDE.B=[0 0 1 0;0 1 0 0];
%             PDE.Bw=[1 0 0 0;0 1 0 0];
% 
%           % Non-polynomial in space forcing needs to be entered through
%             % symbolic PDE.B21_nonpol matrix
%             PDE.B21_nonpol=sx*zeros(PDE.n0+PDE.n1+PDE.n2,PDE.nw);  
%             PDE.B21_nonpol(1,3:4)=[sin(pi*sx) -PDE.A0*sin(pi*sx)-PDE.A1*pi*cos(pi*sx)+PDE.A2*pi^2*sin(pi*sx)];
%            
            % Terms format
                  PDE.n.n_pde=[0,0,1];
                  PDE.n.nw=4;
                  PDE.n.nv=4;

                  PDE.PDE.A{1}.Lstate=2;
                  PDE.PDE.A{1}.Rstate=2;
                  PDE.PDE.A{1}.D=2;
                  PDE.PDE.A{1}.I=0;
                  PDE.PDE.A{1}.coeff=0.5;

                  PDE.PDE.A{2}.Lstate=2;
                  PDE.PDE.A{2}.Rstate=2;
                  PDE.PDE.A{2}.D=1;
                  PDE.PDE.A{2}.I=0;
                  PDE.PDE.A{2}.coeff=2;

                  PDE.PDE.A{3}.Lstate=2;
                  PDE.PDE.A{3}.Rstate=2;
                  PDE.PDE.A{3}.I=0;
                  PDE.PDE.A{3}.coeff=4;
% % % % % 
% % % % %        % BCs: u_x(x=a)       
                 PDE.BC.Ebb{1}.coeff = [eye(1);zeros(1)];
                 PDE.BC.Ebb{1}.Rstate = 2; PDE.BC.Ebb{1}.delta = 0;
                 PDE.BC.Ebb{1}.D = 1;
% % % % %  
% % % % %        % BCs: u(x=b)
                 PDE.BC.Ebb{2}.coeff = [zeros(1);eye(1)];
                 PDE.BC.Ebb{2}.Rstate = 2; PDE.BC.Ebb{2}.delta = 1; 
% % % % 
% % % % %          % BCs: inhomogeneous inputs (connect through ODE interconnected signals)
                 PDE.BC.Ebv=[1 0 0 0;0 1 0 0];
                 PDE.ODE.Dvw=eye(4);

            % Non-polynomial in space forcing needs to be entered through
            % symbolic PDE.PDE.Bpw_nonpol matrix 
            PDE.PDE.Bpw_nonpol=sx*zeros(sum(PDE.n.n_pde),PDE.n.nw);  
            PDE.PDE.Bpw_nonpol(1,3:4)=[sin(pi*sx) -PDE.PDE.A{3}.coeff*sin(pi*sx)-PDE.PDE.A{2}.coeff*pi*cos(pi*sx)+PDE.PDE.A{1}.coeff*pi^2*sin(pi*sx)];
            
%           Exact solution, initial conditions and inhomogeneous inputs  
            % Initial conditions for the primary states of the PDE
             uinput.ic.PDE(1)= sin(pi*sx);

             % Forcing
             uinput.w(1)= (st+1)^0.5*pi*cos(pi*a);
             uinput.w(2)= (st+1)^0.5*sin(pi*b);
             uinput.w(3)=0.5*(st+1)^(-0.5);
             uinput.w(4)=(st+1)^0.5;
             
            % Exact solution (optional)

            uinput.exact(1) = (st+1)^0.5*sin(pi*sx);

%-------------------------------------------------------------
% Euler-Bernoulli beam equation examples u_tt=-c u_xxxx
%-------------------------------------------------------------
%-------------------------------------------------------------
% Example 14 - exact eigenmodes of a cantilever beam vibration
%-------------------------------------------------------------
             case 14
                 
 % NOTE: don't change the domain limits for this test case. Beta values
 % below correspond to the beam of length L=2 only

% % % %  EB-beam equation in second-order formulation
% % % %  NOTE: this formulation provides the highest accuracy (solution is
% % % %  more accurate at lower N compared to the 4th-order formulations in examples 15 and 16)

% % %    Primary states are state(1)=u_t; state(2)=u_xx; 
%        Fundamental states are state(1)_xx=u_txx, 
%        state(2)_xx=u_xxxx.
%        Boundary conditions
%        Dirichlet: state(1): u_t at x=a
%        Dirichlet: state(2): u_xx at x=b
%        Neumann:   state(1)_x: u_tx at x=a
%        Neumann:   state(2)_x: u_xxx at x=b

            a=0;b=2; PDE.dom=[a b]; 
            L=2;
            c=2;
            % Batch format
            PDE.n0=0; PDE.n1=0; PDE.n2=2; PDE.nw=0; 
            PDE.A0=0; PDE.A1=0; PDE.A2=[0 -c;1 0];
            PDE.B=[1 0 0 0 0 0 0 0;0 0 0 1 0 0 0 0; 0 0 0 0 1 0 0 0; 0 0 0 0 0 0 0 1];

%             % Terms format
%             PDE.n.n_pde=[0,0,2];
%             PDE.PDE.A{1}.Lstate=2;
%             PDE.PDE.A{1}.Rstate=2;
%             PDE.PDE.A{1}.I=0;
%             PDE.PDE.A{1}.coeff=[0 -c; 1 0];
%             PDE.PDE.A{1}.D=2;

%             % BCs: 0 = u1(a)
%              PDE.BC.Ebb{1}.coeff = [1,0;0,0;0,0;0,0];
%              PDE.BC.Ebb{1}.Rstate = 2; PDE.BC.Ebb{1}.delta = 0; 

%              % BCs: 0 = u2(b)
%              PDE.BC.Ebb{2}.coeff = [0,0;0,1;0,0;0,0];
%              PDE.BC.Ebb{2}.Rstate = 2; PDE.BC.Ebb{2}.delta = 1;           
            
%              % BCs: 0 = u1_x(a)
%              PDE.BC.Ebb{3}.coeff = [0,0;0,0;1,0;0,0];
%              PDE.BC.Ebb{3}.Rstate = 2; PDE.BC.Ebb{3}.delta = 0; PDE.BC.Ebb{3}.D = 1;
%             
%              % BCs: 0 = u2_x(b)
%              PDE.BC.Ebb{4}.coeff = [0,0;0,0;0,0;0,1];
%              PDE.BC.Ebb{4}.Rstate = 2; PDE.BC.Ebb{4}.delta = 1; PDE.BC.Ebb{4}.D = 1;

%           Exact solution, initial conditions and inhomogeneous inputs  
% %           
% %             
%               beta=0.937552034355981;  % First eigenmode
%               beta=2.347045566487088;  % Second eigenmode
%               beta=3.927378719118806;  % Third eigenmode
                beta=5.497770367437734;  % Fourth eigenmode
                omega=beta^2*sqrt(c);
                coeff=(cos(beta*L)+cosh(beta*L))/(sin(beta*L)+sinh(beta*L));
              
%             exact solution is as follows:
             mode = cosh(beta*sx)-cos(beta*sx)+coeff*(sin(beta*sx)-sinh(beta*sx));
             u = mode*cos(omega*st);
% % %           
% %  %        exact solution for primary states u_t(x,t), u_xx(x,t)
             uinput.exact(1) =  diff(u,st);
             uinput.exact(2) =  diff(u,sx,2);

% % % % % % % Initial conditions are on primary states u_t and u_xx, which are exact solutions
% % % % % % % evaluated at time t=0;

               uinput.ic.PDE(1) =  subs(uinput.exact(1),st,0);
               uinput.ic.PDE(2) =  subs(uinput.exact(2),st,0);


%-------------------------------------------------------------
% Example 15 - exact eigenmodes of a cantilever beam vibration
% Formulation as a fourth-order equation (available only in terms format)  
%-------------------------------------------------------------
               case 15

 % Euler-Bernoulli beam equation examples u_tt=-c u_xxxx
 % Used states are u1=u_xxxx, u2=u_t
 % Dynamics
 % u1_t = u2_xxxx
 % u2_t = -c u_1
                 
 % NOTE: don't change the domain limits for this test case. Beta values
 % below correspond to the beam of length L=2 only

            a=0;b=2; PDE.dom=[a b]; 
            L=2;
            c=2;
  % Terms format
            PDE.n.n_pde=[1,0,0,0,1];

            PDE.PDE.A{1}.Lstate=0;
            PDE.PDE.A{1}.Rstate=4;
            PDE.PDE.A{1}.I=0;
            PDE.PDE.A{1}.coeff=1;
            PDE.PDE.A{1}.D=4;

            PDE.PDE.A{2}.Lstate=4;
            PDE.PDE.A{2}.Rstate=0;
            PDE.PDE.A{2}.I=0;
            PDE.PDE.A{2}.coeff=-c;

%           Boundary conditions: 
% %            u(x=a) = 0
% %            u_x(x=a) = 0 
% %            u_xx(x=b) = 0 
% %            u_xxx(x=b) = 0 

%          Formulated in terms of used states: all four boundary conditions
%          are on state u2=u_t
%          u1(x=a)=0; u1_x(x=1)=0; u1_xx(x=b)=0; u1_xxx(x=b)=0;

%          BCs u1(x=a)=0
           PDE.BC.Ebb{1}.coeff = [1;0;0;0];
           PDE.BC.Ebb{1}.Rstate = 4; PDE.BC.Ebb{1}.delta = 0; 

%          BCs u1_x(x=a)=0
           PDE.BC.Ebb{2}.coeff = [0;1;0;0];
           PDE.BC.Ebb{2}.Rstate = 4; PDE.BC.Ebb{2}.delta = 0; PDE.BC.Ebb{2}.D = 1;

%          BCs u1_xx(x=b)=0
           PDE.BC.Ebb{3}.coeff = [0;0;1;0];
           PDE.BC.Ebb{3}.Rstate = 4; PDE.BC.Ebb{3}.delta = 1; PDE.BC.Ebb{3}.D = 2; 

%          BCs u1_xxx(x=b)=0
           PDE.BC.Ebb{4}.coeff = [0;0;0;1];
           PDE.BC.Ebb{4}.Rstate = 4; PDE.BC.Ebb{4}.delta = 1; PDE.BC.Ebb{4}.D = 3; 

%            Exact solution, initial conditions and inhomogeneous inputs  
% %           
% %             
%               beta=0.937552034355981;  % First eigenmode
%               beta=2.347045566487088;  % Second eigenmode
%               beta=3.927378719118806;  % Third eigenmode
                beta=5.497770367437734;  % Fourth eigenmode
                omega=beta^2*sqrt(c);
                coeff=(cos(beta*L)+cosh(beta*L))/(sin(beta*L)+sinh(beta*L));
              
%             exact solution is as follows:
             mode = cosh(beta*sx)-cos(beta*sx)+coeff*(sin(beta*sx)-sinh(beta*sx));
             u = mode*cos(omega*st);
% % %           

% %  %        exact solution for primary states u_xxxx(x,t), u_t(x,t)
%            Note: states with lower differentiability should be listed first
%            (i.e. n0 states before n1 states before n2 states etc.)
             uinput.exact(1) =  diff(u,sx,4);
             uinput.exact(2) =  diff(u,st);


% % % % % % % Initial conditions are on primary states u_t and u_xxxx, which are exact solutions
% % % % % % % evaluated at time t=0;

               uinput.ic.PDE(1) =  subs(uinput.exact(1),st,0);
               uinput.ic.PDE(2) =  subs(uinput.exact(2),st,0);

%-------------------------------------------------------------
% Example 16 - exact eigenmodes of a cantilever beam vibration
% Another formulation as a fourth-order equation (available only in terms format)  
% This formulation would not require a back reconstruction of the original
% solution out of the derivatives
%-------------------------------------------------------------
               case 16

 % Euler-Bernoulli beam equation examples u_tt=-c u_xxxx
 % Used states are u1=u, u2=u_t
 % Dynamics
 % u1_t = u2
 % u2_t = -c u_1_xxxx
                 
 % NOTE: don't change the domain limits for this test case. Beta values
 % below correspond to the beam of length L=2 only

            a=0;b=2; PDE.dom=[a b]; 
            L=2;
            c=2;
  % Terms format
            PDE.n.n_pde=[0,0,0,0,2];

            PDE.PDE.A{1}.Lstate=4;
            PDE.PDE.A{1}.Rstate=4;
            PDE.PDE.A{1}.I=0;
            PDE.PDE.A{1}.coeff=[0 1;0 0];
            PDE.PDE.A{1}.D=0;

            PDE.PDE.A{2}.Lstate=4;
            PDE.PDE.A{2}.Rstate=4;
            PDE.PDE.A{2}.I=0;
            PDE.PDE.A{2}.coeff=[0 0;-c 0];
            PDE.PDE.A{2}.D=4;

%           Boundary conditions: 
% %            u(x=a) = 0
% %            u_x(x=a) = 0 
% %            u_xx(x=b) = 0 
% %            u_xxx(x=b) = 0 

%          Formulated in terms of used states: 
%          u1(x=a)=0; u1_x(x=1)=0; u1_xx(x=b)=0; u1_xxx(x=b)=0;
%          u2(x=a)=0; u2_x(x=1)=0; u2_xx(x=b)=0; u2_xxx(x=b)=0;

%          BCs u1(x=a)=0
           PDE.BC.Ebb{1}.coeff = [1 0;0 0;0 0;0 0;0 0;0 0;0 0;0 0];
           PDE.BC.Ebb{1}.Rstate = 4; PDE.BC.Ebb{1}.delta = 0; 

%          BCs u2(x=a)=0
           PDE.BC.Ebb{2}.coeff = [0 0;0 1;0 0;0 0;0 0;0 0;0 0;0 0];
           PDE.BC.Ebb{2}.Rstate = 4; PDE.BC.Ebb{2}.delta = 0; 

%          BCs u1_x(x=a)=0
           PDE.BC.Ebb{3}.coeff = [0 0;0 0;1 0;0 0;0 0;0 0;0 0;0 0];
           PDE.BC.Ebb{3}.Rstate = 4; PDE.BC.Ebb{3}.delta = 0; PDE.BC.Ebb{3}.D = 1;

%          BCs u2_x(x=a)=0
           PDE.BC.Ebb{4}.coeff = [0 0;0 0;0 0;0 1;0 0;0 0;0 0;0 0];
           PDE.BC.Ebb{4}.Rstate = 4; PDE.BC.Ebb{4}.delta = 0; PDE.BC.Ebb{4}.D = 1;

%          BCs u1_xx(x=b)=0
           PDE.BC.Ebb{5}.coeff = [0 0;0 0;0 0;0 0;1 0;0 0;0 0;0 0];
           PDE.BC.Ebb{5}.Rstate = 4; PDE.BC.Ebb{5}.delta = 1; PDE.BC.Ebb{5}.D = 2; 

%          BCs u2_xx(x=b)=0
           PDE.BC.Ebb{6}.coeff = [0 0;0 0;0 0;0 0;0 0;0 1;0 0;0 0];
           PDE.BC.Ebb{6}.Rstate = 4; PDE.BC.Ebb{6}.delta = 1; PDE.BC.Ebb{6}.D = 2; 

%          BCs u1_xxx(x=b)=0
           PDE.BC.Ebb{7}.coeff = [0 0;0 0;0 0;0 0;0 0;0 0;1 0;0 0];
           PDE.BC.Ebb{7}.Rstate = 4; PDE.BC.Ebb{7}.delta = 1; PDE.BC.Ebb{7}.D = 3; 

%          BCs u2_xxx(x=b)=0
           PDE.BC.Ebb{8}.coeff = [0 0;0 0;0 0;0 0;0 0;0 0;0 0;0 1];
           PDE.BC.Ebb{8}.Rstate = 4; PDE.BC.Ebb{8}.delta = 1; PDE.BC.Ebb{8}.D = 3; 

%            Exact solution, initial conditions and inhomogeneous inputs  
% %           
% %             
%               beta=0.937552034355981;  % First eigenmode
%               beta=2.347045566487088;  % Second eigenmode
%               beta=3.927378719118806;  % Third eigenmode
                beta=5.497770367437734;  % Fourth eigenmode
                omega=beta^2*sqrt(c);
                coeff=(cos(beta*L)+cosh(beta*L))/(sin(beta*L)+sinh(beta*L));
              
%             exact solution is as follows:
             mode = cosh(beta*sx)-cos(beta*sx)+coeff*(sin(beta*sx)-sinh(beta*sx));
             u = mode*cos(omega*st);
% % %           

% %  %       exact solution for primary states u(x,t), u_t(x,t)
%            Note: states with lower differentiability should be listed first
%            (i.e. n0 states before n1 states before n2 states etc.)
             uinput.exact(1) =  u;
             uinput.exact(2) =  diff(u,st);


% % % % % % % Initial conditions are on primary states u_t and u_xxxx, which are exact solutions
% % % % % % % evaluated at time t=0;

               uinput.ic.PDE(1) =  subs(uinput.exact(1),st,0);
               uinput.ic.PDE(2) =  subs(uinput.exact(2),st,0);



%-------------------------------------------------------------
% Forced Euler-Bernoulli beam equation examples u_tt=-c u_xxxx+f(x,t)
%-------------------------------------------------------------
%-------------------------------------------------------------
% Example 17 - exact solution with forcing, solution is polynomial in time
%-------------------------------------------------------------
    case 17

% % %       u(x,t)=x^3*t^2 solution for u_tt=-c u_xxxx+2x^3;
% % Primary states are state(1)=u_t; state(2)=u_xx; Fundamental states are 
%   state(1)_xx=u_txx, state(2)_xx=u_xxxx.
%             

            a=-1;b=1; PDE.dom=[a b];   
            c=2;

%             % Batch format
%              PDE.n0=0; PDE.n1=0; PDE.n2=2; PDE.nw=5;
%              PDE.A0=0; PDE.A1=0; PDE.A2=[0 -c;1 0];
%              PDE.B=[1 0 0 0 0 0 0 0;0 0 0 1 0 0 0 0; 0 0 0 0 1 0 0 0; 0 0 0 0 0 0 0 1];
%              PDE.Bw=[1 0 0 0 0;0 1 0 0 0; 0 0 1 0 0;0 0 0 1 0];

              % PDE.B21=s*zeros(PDE.n0+PDE.n1+PDE.n2,PDE.nw);  
              % PDE.B21(1,5)=2*s^3;
            
            % Terms format
            PDE.n.n_pde=[0,0,2];
            PDE.n.nw=5;
            PDE.n.nv=5;
            PDE.PDE.A{1}.Lstate=2;
            PDE.PDE.A{1}.Rstate=2;
            PDE.PDE.A{1}.I=0;
            PDE.PDE.A{1}.coeff=[0 -c; 1 0];
            PDE.PDE.A{1}.D=2;

            PDE.BC.Ebv=[1 0 0 0 0;0 1 0 0 0; 0 0 1 0 0;0 0 0 1 0];
            PDE.ODE.Dvw=eye(5);

            PDE.PDE.Bpv=s1*zeros(sum(PDE.n.n_pde),PDE.n.nv);  
            PDE.PDE.Bpv(1,5)=2*s1^3;

            %  Boundary conditions
            % BCs: 0 = u1(a)
             PDE.BC.Ebb{1}.coeff = [1,0;0,0;0,0;0,0];
             PDE.BC.Ebb{1}.Rstate = 2; PDE.BC.Ebb{1}.delta = 0; 

             % BCs: 0 = u2(b)
             PDE.BC.Ebb{2}.coeff = [0,0;0,1;0,0;0,0];
             PDE.BC.Ebb{2}.Rstate = 2; PDE.BC.Ebb{2}.delta = 1;           
            
             % BCs: 0 = u1_x(a)
             PDE.BC.Ebb{3}.coeff = [0,0;0,0;1,0;0,0];
             PDE.BC.Ebb{3}.Rstate = 2; PDE.BC.Ebb{3}.delta = 0; PDE.BC.Ebb{3}.D = 1;
            
             % BCs: 0 = u2_x(b)
             PDE.BC.Ebb{4}.coeff = [0,0;0,0;0,0;0,1];
             PDE.BC.Ebb{4}.Rstate = 2; PDE.BC.Ebb{4}.delta = 1; PDE.BC.Ebb{4}.D = 1;


%           Exact solution, initial conditions and inhomogeneous inputs  
            % exact solution of the PDE
            u = sx.^3*st^2;   
% %           
% %  %        exact solution for primary states u_t(x,t), u_xx(x,t)
             uinput.exact(1) =  diff(u,st);
             uinput.exact(2) =  diff(u,sx,2);
% %           
% % % % % % % Initial conditions on primary states are exact solutions
% % % % % % % for primary states u_t and u_xx evaluated at time t=0;

               uinput.ic.PDE(1)= subs(uinput.exact(1),st,0);
               uinput.ic.PDE(2)= subs(uinput.exact(2),st,0);
               
               % Boundary conditions: need four boundary conditions on primary states 
% %           
               uinput.w(1)= 2*a^3*st;       % Dirichlet: state(1): u_t at x=a
               uinput.w(2)= 6*b*st^2;       % Dirichlet: state(2): u_xx at x=b
               uinput.w(3)= 6*a^2*st;       % Neumann:   state(1)_x: u_tx at x=a
               uinput.w(4)= 6*st^2;         % Neumann:   state(2)_x: u_xxx at x=b

                % If forcing is polynomial in space, it is done through B21
                % matrix (in batch format) or through PDE.PDE.Bpv (in terms
                % format)
                % 5th input is distributed forcing
               uinput.w(5)=1;
% %            


%-------------------------------------------------------------
% Example 18 - exact solution with forcing, solution is polynomial in time.
% Fromulates as 4th order equation
%-------------------------------------------------------------
    case 18

 % Euler-Bernoulli beam equation u_tt=-c u_xxxx+f(x,t)
 % Used states are u1=u, u2=u_t
 % Dynamics
 % u1_t = u2
 % u2_t = -c u_1_xxxx+f(x,t)

            a=-1;b=1; PDE.dom=[a b];   
            c=2;

  % Terms format
            PDE.n.n_pde=[0,0,0,0,2];

            PDE.n.nw=9;
            PDE.n.nv=9;

            PDE.PDE.A{1}.Lstate=4;
            PDE.PDE.A{1}.Rstate=4;
            PDE.PDE.A{1}.I=0;
            PDE.PDE.A{1}.coeff=[0 1;0 0];
            PDE.PDE.A{1}.D=0;

            PDE.PDE.A{2}.Lstate=4;
            PDE.PDE.A{2}.Rstate=4;
            PDE.PDE.A{2}.I=0;
            PDE.PDE.A{2}.coeff=[0 0;-c 0];
            PDE.PDE.A{2}.D=4;

             PDE.BC.Ebv=zeros(PDE.n.nv-1,PDE.n.nv);
             PDE.BC.Ebv(1:PDE.n.nv-1,1:PDE.n.nv-1)=eye(PDE.n.nv-1);
             PDE.ODE.Dvw=eye(PDE.n.nv);
% 
             PDE.PDE.Bpv=s1*zeros(sum(PDE.n.n_pde),PDE.n.nv);  
             PDE.PDE.Bpv(2,PDE.n.nv)=2*s1^3;


%           Boundary conditions: 
% %            u(x=a) 
% %            u_x(x=a) 
% %            u_xx(x=b) 
% %            u_xxx(x=b) 

%          Formulated in terms of used states: 
%          u1(x=a); u1_x(x=1); u1_xx(x=b); u1_xxx(x=b);
%          u2(x=a); u2_x(x=1); u2_xx(x=b); u2_xxx(x=b);

%          BCs u1(x=a)=0
           PDE.BC.Ebb{1}.coeff = [1 0;0 0;0 0;0 0;0 0;0 0;0 0;0 0];
           PDE.BC.Ebb{1}.Rstate = 4; PDE.BC.Ebb{1}.delta = 0; 

%          BCs u2(x=a)=0
           PDE.BC.Ebb{2}.coeff = [0 0;0 1;0 0;0 0;0 0;0 0;0 0;0 0];
           PDE.BC.Ebb{2}.Rstate = 4; PDE.BC.Ebb{2}.delta = 0; 

%          BCs u1_x(x=a)=0
           PDE.BC.Ebb{3}.coeff = [0 0;0 0;1 0;0 0;0 0;0 0;0 0;0 0];
           PDE.BC.Ebb{3}.Rstate = 4; PDE.BC.Ebb{3}.delta = 0; PDE.BC.Ebb{3}.D = 1;

%          BCs u2_x(x=a)=0
           PDE.BC.Ebb{4}.coeff = [0 0;0 0;0 0;0 1;0 0;0 0;0 0;0 0];
           PDE.BC.Ebb{4}.Rstate = 4; PDE.BC.Ebb{4}.delta = 0; PDE.BC.Ebb{4}.D = 1;

%          BCs u1_xx(x=b)=0
           PDE.BC.Ebb{5}.coeff = [0 0;0 0;0 0;0 0;1 0;0 0;0 0;0 0];
           PDE.BC.Ebb{5}.Rstate = 4; PDE.BC.Ebb{5}.delta = 1; PDE.BC.Ebb{5}.D = 2; 

%          BCs u2_xx(x=b)=0
           PDE.BC.Ebb{6}.coeff = [0 0;0 0;0 0;0 0;0 0;0 1;0 0;0 0];
           PDE.BC.Ebb{6}.Rstate = 4; PDE.BC.Ebb{6}.delta = 1; PDE.BC.Ebb{6}.D = 2; 

%          BCs u1_xxx(x=b)=0
           PDE.BC.Ebb{7}.coeff = [0 0;0 0;0 0;0 0;0 0;0 0;1 0;0 0];
           PDE.BC.Ebb{7}.Rstate = 4; PDE.BC.Ebb{7}.delta = 1; PDE.BC.Ebb{7}.D = 3; 

%          BCs u2_xxx(x=b)=0
           PDE.BC.Ebb{8}.coeff = [0 0;0 0;0 0;0 0;0 0;0 0;0 0;0 1];
           PDE.BC.Ebb{8}.Rstate = 4; PDE.BC.Ebb{8}.delta = 1; PDE.BC.Ebb{8}.D = 3; 


%           Exact solution, initial conditions and inhomogeneous inputs  
            % exact solution of the PDE
            u = sx.^3*st^2;   
            u_t=diff(u,st);
% %           
% %  %        exact solution for primary states u_t(x,t), u_xx(x,t)
             uinput.exact(1) =  u;
             uinput.exact(2) =  u_t;
% %           
% % % % % % % Initial conditions on primary states are exact solutions
% % % % % % % for primary states u_t and u_xx evaluated at time t=0;

               uinput.ic.PDE(1)= subs(uinput.exact(1),st,0);
               uinput.ic.PDE(2)= subs(uinput.exact(2),st,0);
               
               % Boundary conditions: need four boundary conditions on primary states 
% %           
               uinput.w(1)=subs(u,sx,a); 
               uinput.w(2)=subs(u_t,sx,a); 
               uinput.w(3)=subs(diff(u,sx),sx,a); 
               uinput.w(4)= subs(diff(u_t,sx),sx,a);  
               uinput.w(5)=subs(diff(u,sx,2),sx,b); 
               uinput.w(6)= subs(diff(u_t,sx,2),sx,b);   
               uinput.w(7)=subs(diff(u,sx,3),sx,b); 
               uinput.w(8)= subs(diff(u_t,sx,3),sx,b);        

                % If forcing is polynomial in space, it is done through B21
                % matrix (in batch format) or through PDE.PDE.Bpv (in terms
                % format)
                % Last input is distributed forcing
               uinput.w(PDE.n.nv)=1;

%--------------------------------------------
% Example 19 - exact solution with forcing, solution is exponentially decaying in time
%--------------------------------------------
    case 19

% % %       u(x,t)=sin x e^(-t) solution for u_tt=-c u_xxxx+(c+1)sino e^(-t);
%   Primary states are state(1)=u_t; state(2)=u_xx; Fundamental states are 
%   state(1)_xx=u_txx, state(2)_xx=u_xxxx.

            a=2;b=3;PDE.dom=[a b];   
            c=2;

            % Batch format
%             PDE.n0=0; PDE.n1=0; PDE.n2=2; 
%             PDE.nw=5;
%             PDE.A0=0; PDE.A1=0; PDE.A2=[0 -c;1 0];
%             PDE.B=[1 0 0 0 0 0 0 0;0 0 0 1 0 0 0 0; 0 0 0 0 1 0 0 0; 0 0 0 0 0 0 0 1];
%             PDE.Bw=[eye(4) zeros(4,1)];
% 
%           % Non-polynomial in space forcing needs to be entered through
%             % symbolic PDE.B21_nonpol matrix
%             PDE.B21_nonpol=sx*zeros(PDE.n0+PDE.n1+PDE.n2,PDE.nw);  
%             PDE.B21_nonpol(1,5)=(c+1)*sin(sx);

            % Terms format

            PDE.n.n_pde=[0,0,2];
            PDE.n.nw=5;
            PDE.n.nv=5;
            PDE.PDE.A{1}.Lstate=2;
            PDE.PDE.A{1}.Rstate=2;
            PDE.PDE.A{1}.I=0;
            PDE.PDE.A{1}.coeff=[0 -c; 1 0];
            PDE.PDE.A{1}.D=2;

            PDE.BC.Ebv=[1 0 0 0 0;0 1 0 0 0; 0 0 1 0 0;0 0 0 1 0];
            PDE.ODE.Dvw=eye(5);

             %  Boundary conditions
            % BCs: 0 = u1(a)
             PDE.BC.Ebb{1}.coeff = [1,0;0,0;0,0;0,0];
             PDE.BC.Ebb{1}.Rstate = 2; PDE.BC.Ebb{1}.delta = 0; 

             % BCs: 0 = u2(b)
             PDE.BC.Ebb{2}.coeff = [0,0;0,1;0,0;0,0];
             PDE.BC.Ebb{2}.Rstate = 2; PDE.BC.Ebb{2}.delta = 1;           
            
             % BCs: 0 = u1_x(a)
             PDE.BC.Ebb{3}.coeff = [0,0;0,0;1,0;0,0];
             PDE.BC.Ebb{3}.Rstate = 2; PDE.BC.Ebb{3}.delta = 0; PDE.BC.Ebb{3}.D = 1;
            
             % BCs: 0 = u2_x(b)
             PDE.BC.Ebb{4}.coeff = [0,0;0,0;0,0;0,1];
             PDE.BC.Ebb{4}.Rstate = 2; PDE.BC.Ebb{4}.delta = 1; PDE.BC.Ebb{4}.D = 1;


            % Non-polynomial in space forcing needs to be entered through
            % symbolic PDE.PDE.Bpw_nonpol matrix
            PDE.PDE.Bpw_nonpol=sx*zeros(sum(PDE.n.n_pde),PDE.n.nw);  
            PDE.PDE.Bpw_nonpol(1,5)=(c+1)*sin(sx);
            
%           Exact solution, initial conditions and inhomogeneous inputs 

% exact solution of the PDE
            
            u = sin(sx)*exp(-st);   
            
% %  %        exact solution for primary states u_t(x,t), u_xx(x,t)
             uinput.exact(1) =  diff(u,st);
             uinput.exact(2) =  diff(u,sx,2);
% %           

% % % % % % % Initial conditions on primary states are exact solutions
% % % % % % % for primary states u_t and u_xx evaluated at time t=0;
               
               uinput.ic.PDE(1)= subs(uinput.exact(1),st,0);
               uinput.ic.PDE(2)= subs(uinput.exact(2),st,0);
               
          
               % Boundary conditions: need four boundary conditions on primary states 
% %           
               uinput.w(1)= -sin(a)*exp(-st);  % Dirichlet: state(1): u_t at x=a
               uinput.w(2)= -sin(b)*exp(-st);  % Dirichlet: state(2): u_xx at x=b
               uinput.w(3)= -cos(a)*exp(-st); % Neumann:   state(1)_x: u_tx at x=a
               uinput.w(4)= -cos(b)*exp(-st);  % Neumann:   state(2)_x: u_xxx at x=b
               uinput.w(5)= exp(-st);


%--------------------------------------------
% Example 20 - exact solution with forcing, solution is exponentially decaying in time
% Formulated as fourth order equation
%--------------------------------------------
    case 20

        % Euler-Bernoulli beam equation u_tt=-c u_xxxx+f(x,t)
 % Used states are u1=u, u2=u_t
 % Dynamics
 % u1_t = u2
 % u2_t = -c u_1_xxxx+f(x,t)

            a=-1;b=1; PDE.dom=[a b];   
            c=2;

  % Terms format
            PDE.n.n_pde=[0,0,0,0,2];

            PDE.n.nw=9;
            PDE.n.nv=9;

            PDE.PDE.A{1}.Lstate=4;
            PDE.PDE.A{1}.Rstate=4;
            PDE.PDE.A{1}.I=0;
            PDE.PDE.A{1}.coeff=[0 1;0 0];
            PDE.PDE.A{1}.D=0;

            PDE.PDE.A{2}.Lstate=4;
            PDE.PDE.A{2}.Rstate=4;
            PDE.PDE.A{2}.I=0;
            PDE.PDE.A{2}.coeff=[0 0;-c 0];
            PDE.PDE.A{2}.D=4;

             PDE.BC.Ebv=zeros(PDE.n.nv-1,PDE.n.nv);
             PDE.BC.Ebv(1:PDE.n.nv-1,1:PDE.n.nv-1)=eye(PDE.n.nv-1);
             PDE.ODE.Dvw=eye(PDE.n.nv);
% 
             PDE.PDE.Bpw_nonpol=sx*zeros(sum(PDE.n.n_pde),PDE.n.nv);  
             PDE.PDE.Bpw_nonpol(2,PDE.n.nv)=(c+1)*sin(sx);


%           Boundary conditions: 
% %            u(x=a) 
% %            u_x(x=a) 
% %            u_xx(x=b) 
% %            u_xxx(x=b) 

%          Formulated in terms of used states: 
%          u1(x=a); u1_x(x=1); u1_xx(x=b); u1_xxx(x=b);
%          u2(x=a); u2_x(x=1); u2_xx(x=b); u2_xxx(x=b);

%          BCs u1(x=a)=0
           PDE.BC.Ebb{1}.coeff = [1 0;0 0;0 0;0 0;0 0;0 0;0 0;0 0];
           PDE.BC.Ebb{1}.Rstate = 4; PDE.BC.Ebb{1}.delta = 0; 

%          BCs u2(x=a)=0
           PDE.BC.Ebb{2}.coeff = [0 0;0 1;0 0;0 0;0 0;0 0;0 0;0 0];
           PDE.BC.Ebb{2}.Rstate = 4; PDE.BC.Ebb{2}.delta = 0; 

%          BCs u1_x(x=a)=0
           PDE.BC.Ebb{3}.coeff = [0 0;0 0;1 0;0 0;0 0;0 0;0 0;0 0];
           PDE.BC.Ebb{3}.Rstate = 4; PDE.BC.Ebb{3}.delta = 0; PDE.BC.Ebb{3}.D = 1;

%          BCs u2_x(x=a)=0
           PDE.BC.Ebb{4}.coeff = [0 0;0 0;0 0;0 1;0 0;0 0;0 0;0 0];
           PDE.BC.Ebb{4}.Rstate = 4; PDE.BC.Ebb{4}.delta = 0; PDE.BC.Ebb{4}.D = 1;

%          BCs u1_xx(x=b)=0
           PDE.BC.Ebb{5}.coeff = [0 0;0 0;0 0;0 0;1 0;0 0;0 0;0 0];
           PDE.BC.Ebb{5}.Rstate = 4; PDE.BC.Ebb{5}.delta = 1; PDE.BC.Ebb{5}.D = 2; 

%          BCs u2_xx(x=b)=0
           PDE.BC.Ebb{6}.coeff = [0 0;0 0;0 0;0 0;0 0;0 1;0 0;0 0];
           PDE.BC.Ebb{6}.Rstate = 4; PDE.BC.Ebb{6}.delta = 1; PDE.BC.Ebb{6}.D = 2; 

%          BCs u1_xxx(x=b)=0
           PDE.BC.Ebb{7}.coeff = [0 0;0 0;0 0;0 0;0 0;0 0;1 0;0 0];
           PDE.BC.Ebb{7}.Rstate = 4; PDE.BC.Ebb{7}.delta = 1; PDE.BC.Ebb{7}.D = 3; 

%          BCs u2_xxx(x=b)=0
           PDE.BC.Ebb{8}.coeff = [0 0;0 0;0 0;0 0;0 0;0 0;0 0;0 1];
           PDE.BC.Ebb{8}.Rstate = 4; PDE.BC.Ebb{8}.delta = 1; PDE.BC.Ebb{8}.D = 3; 

           % exact solution of the PDE
            
            u = sin(sx)*exp(-st);  
            u_t=diff(u,st);
            
% %  %        exact solution for primary states u(x,t), u_t(x,t)
             uinput.exact(1) = u;
             uinput.exact(2) = u_t;
% %           

% % % % % % % Initial conditions on primary states are exact solutions
% % % % % % % for primary states u_t and u_xx evaluated at time t=0;
               
               uinput.ic.PDE(1)= subs(uinput.exact(1),st,0);
               uinput.ic.PDE(2)= subs(uinput.exact(2),st,0);

               % First 8 inputs are inhomogeneous boundary conditions

               uinput.w(1)=subs(u,sx,a); 
               uinput.w(2)=subs(u_t,sx,a); 
               uinput.w(3)=subs(diff(u,sx),sx,a); 
               uinput.w(4)= subs(diff(u_t,sx),sx,a);  
               uinput.w(5)=subs(diff(u,sx,2),sx,b); 
               uinput.w(6)= subs(diff(u_t,sx,2),sx,b);   
               uinput.w(7)=subs(diff(u,sx,3),sx,b); 
               uinput.w(8)= subs(diff(u_t,sx,3),sx,b);   
                % Last input is distributed forcing
               uinput.w(PDE.n.nv)=exp(-st);


              
%---------------------------------------
% Hyperbolic.PDE equation examples
%----------------------------------------
%---------------------------------------
% Transport equation u_t+c u_x=0
%----------------------------------------
%--------------------------------------------
% Example 21 - traveling wave solution
%--------------------------------------------
    case 21
% u(x,t)=sin(x-ct) solution for u_t+c u_x=0;
% % % % % 
            c=4;
            a=-1;b=1;   PDE.dom=[a b]; 

            % Batch format
            PDE.n0=0; PDE.n1=1; PDE.n2=0;  PDE.nw=1;
            PDE.A0=0; PDE.A1=-c; PDE.A2=0;
            if (c>=0)
            xb=a;
            PDE.B=[1 0];
            else
            xb=b;
            PDE.B=[0 1];
            end
            PDE.Bw=1;

            % Terms format
%              PDE.n.n_pde=[0,1];
%              PDE.n.nw=1;
%              PDE.n.nv=1;
%              PDE.PDE.A{1}.Lstate=1;
%              PDE.PDE.A{1}.Rstate=1;
%              PDE.PDE.A{1}.D=1;
%              PDE.PDE.A{1}.I=0;
%              PDE.PDE.A{1}.coeff=-c;
% 
%               % BCs: u(x=a or x=b depending on sign of c)       
%              PDE.BC.Ebb{1}.coeff = 1;
%              PDE.BC.Ebb{1}.Rstate = 1; 
%              if (c>=0)
%              xb=a;
%              PDE.BC.Ebb{1}.delta = 0;
%              else
%              xb=b;
%              PDE.BC.Ebb{1}.delta = 1;
%              end
            
% % % % 
% % % %      % BCs: inhomogeneous inputs (connect through ODE interconnected signals)
%              PDE.BC.Ebv=1;
%              PDE.ODE.Dvw=1;
                

%        Exact solution, initial conditions and inhomogeneous inputs  
         uinput.exact(1) = sin(sx-c*st);
% % % %  Initial conditions are on primary state;
         uinput.ic.PDE(1)=  sin(sx);
         uinput.w(1)= sin(xb-c*st);
% %       
%    
%--------------------------------------------
% Example 22 - Gauss bump propagating solution
%--------------------------------------------
    case 22
% Gauss exponential 1/(sigma*sqrt(2*pi))*exp(-0.5*(x-ct+1)/sigma)^2) solution for u_t+cu_x=0;
% % 
            sigma=0.2;
            mu=0;
            c=4;
            a=-1;b=1;  PDE.dom=[a b];  

           % Batch format
            PDE.n0=0; PDE.n1=1; PDE.n2=0; PDE.nw=1;
            PDE.A0=0; PDE.A1=-c; PDE.A2=0;
            if (c>=0)
            xb=a;
            PDE.B=[1 0];
            else
            xb=b;
            PDE.B=[0 1];
            end
            PDE.Bw=1;

                        % Terms format
%              PDE.n.n_pde=[0,1];
%              PDE.n.nw=1;
%              PDE.n.nv=1;
%              PDE.PDE.A{1}.Lstate=1;
%              PDE.PDE.A{1}.Rstate=1;
%              PDE.PDE.A{1}.I=0;
%              PDE.PDE.A{1}.coeff=-c;
%              PDE.PDE.A{1}.D=1;

              % BCs: u(x=a or x=b depending on sign of c)       
%              PDE.BC.Ebb{1}.coeff = 1;
%              PDE.BC.Ebb{1}.Rstate = 1; 
%              if (c>=0)
%              xb=a;
%              PDE.BC.Ebb{1}.delta = 0;
%              else
%              xb=b;
%              PDE.BC.Ebb{1}.delta = 1;
%              end
            
% % % % 
% % % %      % BCs: inhomogeneous inputs (connect through ODE interconnected signals)
%              PDE.BC.Ebv=1;
%              PDE.ODE.Dvw=1;
                
%        Exact solution, initial conditions and inhomogeneous inputs 
  uinput.exact(1) = 1/(sigma*sqrt(2*pi))*exp(-0.5*((sx-c*st-mu)/sigma).^2);
% % % %  Initial conditions are on primary state;
  uinput.ic.PDE(1)= subs(uinput.exact(1),st,0);

  uinput.w(1)=1/(sigma*sqrt(2*pi))*exp(-0.5*((xb-c*st-mu)/sigma)^2);
% % % %  
% % 

%---------------------------------------
% Wave equation u_tt=c^2 u_xx
%----------------------------------------
%--------------------------------------------
% Example 23 - standing wave solution
% Boundary conditions 
% B=[1 0 0 0;0 0 0 1], u_t(a,t), u_x(b,t)
%--------------------------------------------

    case 23

% u(x,t)=sin(pi*x)*cos(c*pi*t) solution for u_tt=c^2 u_xx;
% u_t is first primary state, u_x is second primary state
% Boundary conditions: u_t(a,t)=-c*pi*sin(pi*a)*sin(c*pi*t),
% u_x(b,t)=pi*cos(pi*b)*cos(c*pi*st);

            a=-1;b=1; PDE.dom=[a b];
            c=2;  

%           Batch format
%               PDE.n0=0; PDE.n1=2; PDE.n2=0; PDE.nw=2;
%               PDE.A0=0; PDE.A1=[0 c^2; 1 0]; PDE.A2=0;
%               PDE.B=[1 0 0 0;0 0 0 1];
%               PDE.Bw=eye(2);
% 
% %           Terms format
             PDE.n.n_pde=[0,2];
             PDE.n.nw=2;
             PDE.n.nv=2;
             PDE.PDE.A{1}.Lstate=1;
             PDE.PDE.A{1}.Rstate=1;
             PDE.PDE.A{1}.I=0;
             PDE.PDE.A{1}.coeff=[0 c^2; 1 0];
             PDE.PDE.A{1}.D=1;

                 % BCs: Dirichlet on u_t - state 1 (x=a)       
                 PDE.BC.Ebb{1}.coeff = [1 0;0 0];
                 PDE.BC.Ebb{1}.Rstate = 1; PDE.BC.Ebb{1}.delta = 0;
% % % % %  
% % % % %        % BCs: Dirichlet on u_x - state 2 (x=b)
                 PDE.BC.Ebb{2}.coeff = [0 0; 0 1];
                 PDE.BC.Ebb{2}.Rstate = 1; PDE.BC.Ebb{2}.delta = 1; 
% % % % 
% % % %          % BCs: inhomogeneous inputs (connect through ODE interconnected signals)
                 PDE.BC.Ebv=[1 0;0 1];
                 PDE.ODE.Dvw=[1 0;0 1];

%        Exact solution, initial conditions and inhomogeneous inputs 
%         u_t and u_x are the primary states 
          uinput.exact(1) = -c*pi*sin(pi*sx)*sin(c*pi*st);
          uinput.exact(2) = pi*cos(pi*sx)*cos(c*pi*st);
          
%       Initial conditions are on the primary states
          uinput.ic.PDE(1)= subs(uinput.exact(1),st,0);
          uinput.ic.PDE(2)= subs(uinput.exact(2),st,0);
          uinput.w(1)= -c*pi*sin(pi*a)*sin(c*pi*st);
          uinput.w(2)= pi*cos(pi*b)*cos(c*pi*st);

%------------------------------------------

% Example 24 - another standing wave solution
% Boundary conditions

% B=[1 0 0 0;0 0 0 1], u_t(a,t), u_x(b,t)
%--------------------------------------------

    case 24

% % u(x,t)=sin(pi*x)*sin(c*pi*t) solution for u_tt=c^2 u_xx;
% % Boundary conditions: u_t(a,t)=c*pi*sin(pi*a)*cos(c*pi*t), u_x(b,t)=pi*cos(pi*b)*sin(c*pi*t)
% % % % 
            a=-1;b=1;  PDE.dom=[a b];  
            c=-2;

            % Batch format
%             PDE.n0=0; PDE.n1=2; PDE.n2=0; PDE.nw=2;
%             PDE.A0=0; PDE.A1=[0 c^2; 1 0]; PDE.A2=0;
%             PDE.B=[1 0 0 0;0 0 0 1];
%             PDE.Bw=eye(2);

            %   Terms format
             PDE.n.n_pde=[0,2];
             PDE.n.nw=2;
             PDE.n.nv=2;
             PDE.PDE.A{1}.Lstate=1;
             PDE.PDE.A{1}.Rstate=1;
             PDE.PDE.A{1}.I=0;
             PDE.PDE.A{1}.coeff=[0 c^2; 1 0];
             PDE.PDE.A{1}.D=1;

                 % BCs: Dirichlet on u_t - state 1 (x=a)       
                 PDE.BC.Ebb{1}.coeff = [1 0;0 0];
                 PDE.BC.Ebb{1}.Rstate = 1; PDE.BC.Ebb{1}.delta = 0;
% % % % %  
% % % % %        % BCs: Dirichlet on u_x - state 2 (x=b)
                 PDE.BC.Ebb{2}.coeff = [0 0; 0 1];
                 PDE.BC.Ebb{2}.Rstate = 1; PDE.BC.Ebb{2}.delta = 1; 
% % % % 
% % % %          % BCs: inhomogeneous inputs (connect through ODE interconnected signals)
                 PDE.BC.Ebv=[1 0;0 1];
                 PDE.ODE.Dvw=[1 0;0 1];

%        Exact solution, initial conditions and inhomogeneous inputs 
%         ut and ux are the primary states 
          uinput.exact(1) = c*pi*sin(pi*sx)*cos(c*pi*st);
          uinput.exact(2) = pi*cos(pi*sx)*sin(c*pi*st);

%         Initial conditions are on the primary states
          uinput.ic.PDE(1)= subs(uinput.exact(1),st,0);
          uinput.ic.PDE(2)= subs(uinput.exact(2),st,0);
          uinput.w(1)= c*pi*sin(pi*a)*cos(c*pi*st);
          uinput.w(2)= pi*cos(pi*b)*sin(c*pi*st);

%-----------------------------------------
% Example 25 - traveling wave solution
% Boundary conditions

% B=[1 0 0 0;0 0 0 1], u_t(a,t), u_x(b,t)
%--------------------------------------------

    case 25

% %       u(x,t)=sin(x-ct) solution for u_tt=c^2 u_xx;
% % % Boundary conditions: u_t(a,t)=-c*cos(a-ct), u_x(b,t)=cos(b-ct)
            a=-1;b=1;   PDE.dom=[a b]; 
            c=2;

            % Batch format
%             PDE.n0=0; PDE.n1=2; PDE.n2=0; PDE.nw=2;
%             PDE.A0=0; PDE.A1=[0 c^2; 1 0]; PDE.A2=0;
%             PDE.B=[1 0 0 0;0 0 0 1];
%             PDE.Bw=eye(2);

             % Terms format   
             PDE.n.n_pde=[0,2];
             PDE.n.nw=2;
             PDE.n.nv=2;
             PDE.PDE.A{1}.Lstate=1;
             PDE.PDE.A{1}.Rstate=1;
             PDE.PDE.A{1}.I=0;
             PDE.PDE.A{1}.coeff=[0 c^2; 1 0];
             PDE.PDE.A{1}.D=1;

                 % BCs: Dirichlet on u_t - state 1 (x=a)       
                 PDE.BC.Ebb{1}.coeff = [1 0;0 0];
                 PDE.BC.Ebb{1}.Rstate = 1; PDE.BC.Ebb{1}.delta = 0;
% % % % %  
% % % % %        % BCs: Dirichlet on u_x - state 2 (x=b)
                 PDE.BC.Ebb{2}.coeff = [0 0; 0 1];
                 PDE.BC.Ebb{2}.Rstate = 1; PDE.BC.Ebb{2}.delta = 1; 
% % % % 
% % % %          % BCs: inhomogeneous inputs (connect through ODE interconnected signals)
                 PDE.BC.Ebv=[1 0;0 1];
                 PDE.ODE.Dvw=[1 0;0 1];

%        Exact solution, initial conditions and inhomogeneous inputs 
%         u_t and u_x are the primary states 
          uinput.exact(1) = -c*cos(sx-c*st);
          uinput.exact(2) = cos(sx-c*st);
%           
% %       Initial conditions are on the primary states
          uinput.ic.PDE(1)= subs(uinput.exact(1),st,0);
          uinput.ic.PDE(2)= subs(uinput.exact(2),st,0);
          uinput.w(1)= -c*cos(a-c*st);
          uinput.w(2)=  cos(b-c*st);

%-----------------------------------------
% Example 26 - traveling wave solution with characteristic boundary
% conditions at the right end
% Use these boundary conditions with c>0 for enhanced stability 
% Boundary conditions

% B=[1 0 0 0;0 0 1 c], u_t(a,t), u_t(b,t)+cu_x(b,t)=0
%--------------------------------------------

    case 26

% %       u(x,t)=sin(x-ct) solution for u_tt=c^2 u_xx;
% % % Boundary conditions: u_t(a,t)=-c*cos(a-ct), u_t(b,t)+c u_x(b,t)=0;
            a=-1;b=1;  PDE.dom=[a b];  
            c=2;

            % Batch format
%             PDE.n0=0; PDE.n1=2; PDE.n2=0; PDE.nw=1;
%             PDE.A0=0; PDE.A1=[0 c^2; 1 0]; PDE.A2=0;
%             PDE.B=[1 0 0 0;0 0 1 c];
%             PDE.Bw=[1;0];

            % Terms format   
             PDE.n.n_pde=[0,2];
             PDE.n.nw=1;
             PDE.n.nv=1;
             PDE.PDE.A{1}.Lstate=1;
             PDE.PDE.A{1}.Rstate=1;
             PDE.PDE.A{1}.I=0;
             PDE.PDE.A{1}.coeff=[0 c^2; 1 0];
             PDE.PDE.A{1}.D=1;

                 % BCs: Dirichlet on u_t - state 1 (x=a)       
                 PDE.BC.Ebb{1}.coeff = [1 0;0 0];
                 PDE.BC.Ebb{1}.Rstate = 1; PDE.BC.Ebb{1}.delta = 0;
% % % % %  
% % % % %        % BCs: u_t+c u_x =0 (x=b)
                 PDE.BC.Ebb{2}.coeff = [0 0; 1 c];
                 PDE.BC.Ebb{2}.Rstate = 1; PDE.BC.Ebb{2}.delta = 1; 
% % % % 
% % % %          % BCs: inhomogeneous inputs (connect through ODE interconnected signals)
                 PDE.BC.Ebv=[1;0];
                 PDE.ODE.Dvw=1;
% % % % % 
%        Exact solution, initial conditions and inhomogeneous inputs 
%         u_t and u_x are the primary states 
          uinput.exact(1) = -c*cos(sx-c*st);
          uinput.exact(2) = cos(sx-c*st);
%           
%        Initial conditions are on the primary states
          uinput.ic.PDE(1)= subs(uinput.exact(1),st,0);
          uinput.ic.PDE(2)= subs(uinput.exact(2),st,0);
          uinput.w(1)= -c*cos(a-c*st);


%-----------------------------------------
% Example 27 - propagating Gauss bump
% Boundary conditions 

% B=[1 0 0 0;0 0 0 1], u_t(a,t), u_x(b,t)
%--------------------------------------------
    case 27

% Propagating Gauss bump 1/(sigma*sqrt(2*pi))*exp(-0.5*(x-ct-mu)/sigma)^2) solution for u_tt=c^2u_xx;
% % % 

            a=-1;b=1;  PDE.dom=[a b];  
            c=4;

        % Batch format
% 
            PDE.n0=0; PDE.n1=2; PDE.n2=0; PDE.nw=2;
            PDE.A0=0; PDE.A1=[0 c^2; 1 0]; PDE.A2=0;
            PDE.B=[1 0 0 0;0 0 0 1];
            PDE.Bw=eye(2);

%            % Terms format   


%              PDE.n.n_pde=[0,2];
%              PDE.n.nw=2;
%              PDE.n.nv=2;
%              PDE.PDE.A{1}.Lstate=1;
%              PDE.PDE.A{1}.Rstate=1;
%              PDE.PDE.A{1}.I=0;
%              PDE.PDE.A{1}.coeff=[0 c^2; 1 0];
%              PDE.PDE.A{1}.D=1;
% 
%                  % BCs: Dirichlet on u_t - state 1 (x=a)       
%                  PDE.BC.Ebb{1}.coeff = [1 0;0 0];
%                  PDE.BC.Ebb{1}.Rstate = 1; PDE.BC.Ebb{1}.delta = 0;
% % % % % %  
% % % % % %        % BCs: Dirichlet on u_x - state 2 (x=b)
%                  PDE.BC.Ebb{2}.coeff = [0 0; 0 1];
%                  PDE.BC.Ebb{2}.Rstate = 1; PDE.BC.Ebb{2}.delta = 1; 
% % % % % 
% % % % %          % BCs: inhomogeneous inputs (connect through ODE interconnected signals)
%                  PDE.BC.Ebv=[1 0;0 1];
%                  PDE.ODE.Dvw=[1 0;0 1];


            % Case with added regulated outputs

            % Batch format
%                PDE.n0=0; PDE.n1=2; PDE.n2=0; PDE.nw=2;
%                PDE.A0=0; PDE.A1=[0 c^2; 1 0]; PDE.A2=0;
%                PDE.B=[1 0 0 0;0 0 0 1];
%                PDE.Bw=eye(2);
%                PDE.nz=2;
%                PDE.Ca1=[1 2;3 4];
%                PDE.D11=[5 6; 7 8];
%                PDE.nb=1;
%                PDE.Ca2=[1 2];
%                PDE.D21=[-1 3];

%            % Terms format   
%              PDE.n.n_pde=[0,2];
%              PDE.n.nw=2;
%              PDE.n.nv=2;
%              PDE.n.nz=2;
%              PDE.n.nb=1;
%              PDE.n.nr=7;
%              PDE.PDE.A{1}.Lstate=1;
%              PDE.PDE.A{1}.Rstate=1;
%              PDE.PDE.A{1}.I=0;
%              PDE.PDE.A{1}.coeff=[0 c^2; 1 0];
%              PDE.PDE.A{1}.D=1;
% 
%              PDE.PDE.Crp{1}.Rstate=1;
%              PDE.PDE.Crp{1}.coeff=[0 0;0 0;0 0;0 0;1 2;3 4;1 2];
% 
%              PDE.ODE.Dzw=[5 6; 7 8];
%              PDE.ODE.Dyw=[-1 3];
%              PDE.ODE.Dyr=[0 0 0 0 0 0 1];
%              PDE.ODE.Dzr=[0 0 0 0 1 0 0;0 0 0 0 0 1 0];
% 
%                  % BCs: Dirichlet on u_t - state 1 (x=a)       
%                  PDE.BC.Ebb{1}.coeff = [1 0;0 0];
%                  PDE.BC.Ebb{1}.Rstate = 1; PDE.BC.Ebb{1}.delta = 0;
% % % % % %  
% % % % % %        % BCs: Dirichlet on u_x - state 2 (x=b)
%                  PDE.BC.Ebb{2}.coeff = [0 0; 0 1];
%                  PDE.BC.Ebb{2}.Rstate = 1; PDE.BC.Ebb{2}.delta = 1; 
% % % % % 
% % % % %          % BCs: inhomogeneous inputs (connect through ODE interconnected signals)
%                  PDE.BC.Ebv=[1 0;0 1];
%                  PDE.ODE.Dvw=[1 0;0 1];

%        Exact solution, initial conditions and inhomogeneous inputs 
            sigma=0.2;
            mu=0;
% % 
    u = 1/(sigma*sqrt(2*pi))*exp(-0.5*((sx-c*st-mu)/sigma).^2);

%         u_t and u_x are the primary states 
    uinput.exact(1) =c*(sx-c*st-mu)/sigma^2.*u;   
    uinput.exact(2) =-1*(sx-c*st-mu)/sigma^2.*u;
  
%         Initial conditions are on the primary states
          uinput.ic.PDE(1)= subs(uinput.exact(1),st,0);
          uinput.ic.PDE(2)= subs(uinput.exact(2),st,0);
    
    uinput.w(1)=subs(uinput.exact(1),sx,a);
    uinput.w(2)=subs(uinput.exact(2),sx,b);


    
  %-----------------------------------------
% Example 28 - propagating Gauss bump with characteristic boundary
% conditions at the right end
% Use these boundary conditions with c>0 for enhanced stability 
% Boundary conditions 
% B=[1 0 0 0;0 0 1 c], u_t(a,t), u_t(b,t)+cu_x(b,t)=0

%--------------------------------------------
    case 28

% Propagating Gauss bump 1/(sigma*sqrt(2*pi))*exp(-0.5*(x-ct-mu)/sigma)^2) solution for u_tt=c^2u_xx;
% % % 

            a=-1;b=1;  PDE.dom=[a b];  
            c=4;

%            Batch format
            PDE.n0=0; PDE.n1=2; PDE.n2=0; PDE.nw=1;
            PDE.A0=0; PDE.A1=[0 c^2; 1 0]; PDE.A2=0;
            PDE.B=[1 0 0 0;0 0 1 c];
            PDE.Bw=[1;0];

            % Terms format
%             PDE.n.n_pde=[0,2];
%              PDE.n.nw=1;
%              PDE.n.nv=1;
%              PDE.PDE.A{1}.Lstate=1;
%              PDE.PDE.A{1}.Rstate=1;
%              PDE.PDE.A{1}.I=0;
%              PDE.PDE.A{1}.coeff=[0 c^2; 1 0];
%              PDE.PDE.A{1}.D=1;
% 
%                  % BCs: Dirichlet on u_t - state 1 (x=a)       
%                  PDE.BC.Ebb{1}.coeff = [1 0;0 0];
%                  PDE.BC.Ebb{1}.Rstate = 1; PDE.BC.Ebb{1}.delta = 0;
% % % % % %  
% % % % % %        % BCs: u_t+c u_x =0 (x=b)
%                  PDE.BC.Ebb{2}.coeff = [0 0; 1 c];
%                  PDE.BC.Ebb{2}.Rstate = 1; PDE.BC.Ebb{2}.delta = 1; 
% % % % % 
% % % % %          % BCs: inhomogeneous inputs (connect through ODE interconnected signals)
%                  PDE.BC.Ebv=[1;0];
%                  PDE.ODE.Dvw=1;

%        Exact solution, initial conditions and inhomogeneous inputs 
            sigma=0.2;
            mu=0;
% % 
    u = 1/(sigma*sqrt(2*pi))*exp(-0.5*((sx-c*st-mu)/sigma).^2);

%         u_t and u_x are the primary states 
    uinput.exact(1) =c*(sx-c*st-mu)/sigma^2.*u;   
    uinput.exact(2) =-1*(sx-c*st-mu)/sigma^2.*u;

%       Initial conditions are on the primary states
          uinput.ic.PDE(1)= subs(uinput.exact(1),st,0);
          uinput.ic.PDE(2)= subs(uinput.exact(2),st,0);
    
    uinput.w(1)=subs(uinput.exact(1),sx,a); 
  
%-----------------------------------------
%  Example 29 - splitting Gauss bump
%  Boundary conditions
%  B=[1 0 0 0;0 0 0 1], u_t(a,t), u_x(b,t)
%--------------------------------------------
    case 29

% Splitting Gauss bump 1/2(sigma*sqrt(2*pi))*exp(-0.5*(x-ct-mu)/sigma)^2)+1/2(sigma*sqrt(2*pi))*exp(-0.5*(x+ct-mu)/sigma)^2) solution for u_tt=c^2u_xx;

            a=-0.6;b=0.6;  PDE.dom=[a b]; 
            c=-2;

            % Batch format
%             PDE.n0=0; PDE.n1=2; PDE.n2=0; PDE.nw=2;
%             PDE.A0=0; PDE.A1=[0 c^2; 1 0]; PDE.A2=0;
%             PDE.B=[1 0 0 0;0 0 0 1];
%             PDE.Bw=eye(2);

            % Terms format
             PDE.n.n_pde=[0,2];
             PDE.n.nw=2;
             PDE.n.nv=2;
             PDE.PDE.A{1}.Lstate=1;
             PDE.PDE.A{1}.Rstate=1;
             PDE.PDE.A{1}.I=0;
             PDE.PDE.A{1}.coeff=[0 c^2; 1 0];
             PDE.PDE.A{1}.D=1;

                 % BCs: Dirichlet on u_t - state 1 (x=a)       
                 PDE.BC.Ebb{1}.coeff = [1 0;0 0];
                 PDE.BC.Ebb{1}.Rstate = 1; PDE.BC.Ebb{1}.delta = 0;
% % % % %  
% % % % %         BCs: Dirichlet on u_x - state 2 (x=b)
                 PDE.BC.Ebb{2}.coeff = [0 0; 0 1];
                 PDE.BC.Ebb{2}.Rstate = 1; PDE.BC.Ebb{2}.delta = 1; 
% % % % 
% % % %          % BCs: inhomogeneous inputs (connect through ODE interconnected signals)
                 PDE.BC.Ebv=[1 0;0 1];
                 PDE.ODE.Dvw=[1 0;0 1];


%        Exact solution, initial conditions and inhomogeneous inputs 
            sigma=0.2;
            mu=0;
% % 
    uminus =  0.5/(sigma*sqrt(2*pi))*exp(-0.5*((sx-c*st-mu)/sigma).^2);
    uplus=  0.5/(sigma*sqrt(2*pi))*exp(-0.5*((sx+c*st-mu)/sigma).^2);
    u = uminus+uplus;
    
%         u_t and u_x are the primary states 
    uinput.exact(1) = c*(sx-c*st-mu)/sigma^2.*uminus-c*(sx+c*st-mu)/sigma^2.*uplus;   
    uinput.exact(2) = -(sx-c*st-mu)/sigma^2.*uminus-(sx+c*st-mu)/sigma^2.*uplus;
    
%       Initial conditions are on the primary states
          uinput.ic.PDE(1)= subs(uinput.exact(1),st,0);
          uinput.ic.PDE(2)= subs(uinput.exact(2),st,0);
    uinput.w(1)=subs(uinput.exact(1),sx,a);
    uinput.w(2)=subs(uinput.exact(2),sx,b);    
    

%-----------------------------------------
%  Example 30 - This example is added to document a solution of PDEs in a system form. 
%  NOTE: PDEs in this example are uncoupled and are simply solved together to test a matrix implementation.  
%  This example combines u_t=-lambda u equation, transport equation and heat equation
%  into one system. 
%--------------------------------------------
    case 30
        
        % State 1 u(x,t)=u(x,0)exp(-lambda t) solution for u_t=-lambda u;
        
        % State 2 u(x,t)=sin(x-ct) solution for u_t+c u_x=0;
        
        % State 3 u(x,t)=sin(5*pi/4*x+alpha)*exp(-visc*pi^2*t) - solution for u_t=visc*uxx
        
% For state 3: solution function is sin(5*pi/4*x+alpha)
% Dirichlet-Dirichlet boundary conditions

% Boundary conditions: u(a,t)=sin(5*pi/4*a+alpha)*exp(-visc*(5*pi/4)^2*t), u(b,t)=sin(5*pi/4*b+alpha)*exp(-visc*(5*pi/4)^2*t)
% % 
% This would give a zero boundary condition on the left, and time-dependent
% % % boundary condition on the right

            c=2; visc=0.5; lambda=1;
            a=-1;b=1;   PDE.dom=[a b]; 

            % Batch format
%             PDE.n0=1; PDE.n1=1; PDE.n2=1; PDE.nw=3;
%             PDE.A0=[-lambda 0 0;0 0 0;0 0 0]; PDE.A1=[0 0;-c 0;0 0]; PDE.A2=[0;0;visc];
%             PDE.Bw=eye(3);
%            
%   % State 1: n0 state. No boundary conditions required on this state
%             uinput.ic.PDE(1)=cos(sx);
%             uinput.exact(1)=uinput.ic.PDE(1)*exp(-lambda*st);
%               
%   % State 2: n1 state. One boundary condition required on this state
%             if (c>=0)
%             xb=a;
%             PDE.B(1,:)=[1 0 0 0 0 0];
%             else
%             xb=b;
%             PDE.B(1,:)=[0 1 0 0 0 0];
%             end
%             uinput.exact(2) = sin(sx-c*st);
%             uinput.ic.PDE(2)=  sin(sx);
%             uinput.w(1)= sin(xb-c*st);
% % %       
% 
%  % State 3: n2 state. Two boundary conditions required on this state
%            PDE.B(2:3,:)=[0 0 1 0 0 0;0 0 0 1 0 0];
%            alpha=pi/8;
%            uinput.exact(3) =  sin(5*pi/4*sx+alpha)*exp(-visc*(5*pi/4)^2*st);
%            uinput.ic.PDE(3)= subs(uinput.exact(3),st,0);
%            uinput.w(2)=sin(5*pi/4*a+alpha)*exp(-visc*(5*pi/4)^2*st);
%            uinput.w(3)=sin(5*pi/4*b+alpha)*exp(-visc*(5*pi/4)^2*st);

            % Terms format
             PDE.n.n_pde=[1,1,1];
             PDE.n.nw=3;
             PDE.n.nv=3;

             PDE.PDE.A{1}.Lstate=0;
             PDE.PDE.A{1}.Rstate=0;
             PDE.PDE.A{1}.I=0;
             PDE.PDE.A{1}.coeff=-lambda;
             PDE.PDE.A{1}.D=0;

             PDE.PDE.A{2}.Lstate=1;
             PDE.PDE.A{2}.Rstate=1;
             PDE.PDE.A{2}.I=0;
             PDE.PDE.A{2}.coeff=-c;
             PDE.PDE.A{2}.D=1;

             PDE.PDE.A{3}.Lstate=2;
             PDE.PDE.A{3}.Rstate=2;
             PDE.PDE.A{3}.I=0;
             PDE.PDE.A{3}.coeff=visc;
             PDE.PDE.A{3}.D=2;

             % Boundary conditions

%        Exact solution, initial conditions and inhomogeneous inputs 
           
  % State 1: n0 state. No boundary conditions required on this state
            uinput.ic.PDE(1)=cos(sx);
            uinput.exact(1)=uinput.ic.PDE(1)*exp(-lambda*st);
              
  % State 2: n1 state. One boundary condition required on this state     
            PDE.BC.Ebb{1}.coeff = [1;0;0];
            PDE.BC.Ebb{1}.Rstate = 1; 
            if (c>=0)
            xb=a;
            PDE.BC.Ebb{1}.delta = 0;
            else
            xb=b;
            PDE.BC.Ebb{1}.delta = 1;
            end
         uinput.exact(2) = sin(sx-c*st);
         uinput.ic.PDE(2)=  sin(sx);
         uinput.w(1)= sin(xb-c*st);
% %       

 % State 3: n2 state. Two boundary conditions required on this state
            
              % BCs on state n2: u(x=a)       
             PDE.BC.Ebb{2}.coeff = [0; 1; 0];
             PDE.BC.Ebb{2}.Rstate = 2; PDE.BC.Ebb{2}.delta = 0;
  
              % BCs on state n2:  u(x=b)
             PDE.BC.Ebb{3}.coeff = [0;0;1];
             PDE.BC.Ebb{3}.Rstate = 2; PDE.BC.Ebb{3}.delta = 1; 

           alpha=pi/8;
           uinput.exact(3) =  sin(5*pi/4*sx+alpha)*exp(-visc*(5*pi/4)^2*st);
           uinput.ic.PDE(3)= subs(uinput.exact(3),st,0);
           uinput.w(2)=sin(5*pi/4*a+alpha)*exp(-visc*(5*pi/4)^2*st);
           uinput.w(3)=sin(5*pi/4*b+alpha)*exp(-visc*(5*pi/4)^2*st);

           PDE.BC.Ebv=eye(PDE.n.nw);
           PDE.ODE.Dvw=eye(PDE.n.nw);
         
%-------------------------------------------------------------------------------------
%  Example 31 - This example combines u_t=-lambda u equation, 
%  transport equation and heat equation
%  into one system and has forcing
%-------------------------------------------------------------------------------------
    case 31
        
        % State 1 u(x,t)=cos(x)cos(t) solution for u_t=-lambda u+lambda cos(x)cos(t)-cos(x)sin(t);
        
        % State 2 u(x,t)=x^10*t^5 solution for u_t+c u_x=5*x^10*t^4+10*c*x^9*t^5;
        
        % State 3 u(x,t)=sin(pi*x)*sqrt(t+1) - solution for
        % u_t=visc*uxx+f(x,t), where f(x,t)=0.5*(t+t0)^{-1/2}*sin(pi*x)+visc*pi^2*sin(pi*x)*sqrt(t+t0)  
        
            c=10; visc=[1 3]; lambda=1;
            a=2;b=3;   PDE.dom=[a b]; 
            t0=1;

            % Batch format
            PDE.n0=1; PDE.n1=1; PDE.n2=2; 
            PDE.nw=11;
            PDE.no=0; PDE.nu=0;
            PDE.A0=[-lambda 0 0 0;0 0 0 0;0 0 0 0;0 0 0 0]; 
            PDE.A1=[0 0 0; -c 0 0;0 0 0;0 0 0]; 
            PDE.A2=[0 0; 0 0; visc(1) 0;0 visc(2)];
            PDE.Bw=[eye(5) zeros(5,6)];

            % Polynomial in space forcing enters through
            % pvar PDE.B21 matrix

            PDE.B21=s1*zeros(PDE.n0+PDE.n1+PDE.n2,PDE.nw);
            PDE.B21(2,6)=5*s1^10;
            PDE.B21(2,7)=10*c*s1^9;

             % Non-polynomial in space forcing enters through
            % symbolic PDE.B21_nonpol matrix
            PDE.B21_nonpol=sx*zeros(PDE.n0+PDE.n1+PDE.n2,PDE.nw); 
            PDE.B21_nonpol(3:4,8)=sin(pi*sx);
            PDE.B21_nonpol(3,9) =  visc(1)*pi^2*sin(pi*sx);
            PDE.B21_nonpol(4,9) =  visc(2)*pi^2*sin(pi*sx);
            PDE.B21_nonpol(1,10)=lambda*cos(sx);
            PDE.B21_nonpol(1,11)=-cos(sx);

          
  % State 1: n0 state. No boundary conditions required on this state
            uinput.ic.PDE(1)=cos(sx);
            uinput.exact(1)=uinput.ic.PDE(1)*cos(st);
            
  % State 2: n1 state. One boundary condition required on this state  
            
            xb=a;
            PDE.B(1,:)=[1 0 0 0 0 0 0 0 0 0];
           
 
         uinput.exact(2) = sx^10*st^5;
         uinput.ic.PDE(2)=  0;
         
         uinput.w(1)= xb^10*st^5;
% %       

 % State 3: n2 state. Two boundary conditions required on this state
             PDE.B(2:5,:)=[0 0 1 0 0 0 0 0 0 0;0 0 0 0 1 0 0 0 0 0;0 0 0 1 0 0 0 0 0 0;0 0 0 0 0 1 0 0 0 0];
             
           uinput.exact(3) =  sin(pi*sx)*sqrt(st+t0);
           uinput.ic.PDE(3)=  subs(uinput.exact(3),st,0);;
           uinput.exact(4) =  sin(pi*sx)*sqrt(st+t0);
           uinput.ic.PDE(4)= subs(uinput.exact(4),st,0);;
           
           uinput.w(2)=sin(pi*a)*sqrt(st+t0);
           uinput.w(3)=sin(pi*b)*sqrt(st+t0);
           
           uinput.w(4)=sin(pi*a)*sqrt(st+t0);
           uinput.w(5)=sin(pi*b)*sqrt(st+t0);


            uinput.w(6)=st^4;
            uinput.w(7)=st^5;
            uinput.w(8)=0.5/sqrt(st+t0);
            uinput.w(9)=sqrt(st+t0);
            uinput.w(10)=cos(st);
            uinput.w(11)=sin(st);
            
%----------------------------------------------------------------------------
% Examples of PDE/ODE coupling
%----------------------------------------------------------------------------

%----------------------------------------------------------------------------
%Example 32 - 2 PDE states + 2 ODE states. Dynamics uncoupled.
% ODE here is not very interesting. Seems like nothing happens to ODE -
% check
%----------------------------------------------------------------------------

    case 32
            
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

% Forcing is done through B21 matrix
            
            a=-1;b=1; PDE.dom=[a b]; 
            PDE.n0=0; PDE.n1=0; PDE.n2=2; PDE.nw=5; PDE.no=2; PDE.nu=0;
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
               

               % If forcing is polynomial in space, it is done through B21 matrix
               % 5th input is distributed forcing
               uinput.w(5)=1;
               
               
               PDE.B21=s1*zeros(PDE.n0+PDE.n1+PDE.n2,PDE.nw);  
               PDE.B21(1,5)=2*s1^3;
               
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
%Example 33 - Heat Equation coupled with unstable ODE at the boundary 
%----------------------------------------

    case 33

 PDE.n0=0;PDE.n1=0;PDE.n2=1; PDE.no =1;
 a=0; b=1; PDE.dom=[a b]; 
%
 PDE.A0=0; PDE.A1=0; PDE.A2=1; PDE.A = 1; PDE.E =0;
%
%%% % x(a)=0, x(b)=0   - Unstable always (this line should remain commented)
 PDE.B = [1 0 0 0; 0 1 0 0]; PDE.Bx = [0;1];
 
 uinput.ic.ODE=2;
 uinput.ic.PDE=0;
 
 
%----------------------------------------
% Example 34:  x_t = -x(t-tau) %stable if tau<=pi/2 
%----------------------------------------
 

  case 34
  PDE.no = 1; PDE.n0=0; PDE.n1 =1; PDE.n2 =0; 
  PDE.A = 0; PDE.E = 0;
  tau = pi/2 -5e-1; 
  PDE.A0 = 0; PDE.A1 = 1/tau; PDE.A2 = zeros(PDE.n2);
  PDE.E0 = -[1 0]; PDE.Ea = 0; PDE.Eb = 0;
  PDE.B=[0 1];
  PDE.Bx = [1];
  a = -1; b =0; PDE.dom = [-1 0];
  uinput.ic.ODE=[1];
  uinput.ic.PDE=[0];
  
  
%----------------------------------------
% Example 35 - Example 1 from Amritam's paper: stable
%----------------------------------------
    case 35
  % % 
 a=0; b= 1; PDE.dom = [a b];
 PDE.no = 1; PDE.n0=0; PDE.n1=0; PDE.n2 =1;
 PDE.A = -3; 
 PDE.A2 = 1; 
 PDE.B = [0, 0, 1, 0; 0, 1, 0, 0]; 
 PDE.Bx = [0; 0];
 uinput.ic.ODE=[1];
 uinput.ic.PDE=[0];
 
 
%----------------------------------------
% Example 36 - Example 2 from Amritam's paper: PIESIM gets unstable
%----------------------------------------
 
    case 36

 lambda = pi^2-0.01; %(not stable here)
 a=0; b= 1; PDE.dom = [a b];
 PDE.no = 4; PDE.n0=0; PDE.n1=0; PDE.n2 =2;
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
 PDE.nz=2;
 PDE.Ca1=[1 2;3 4];
 PDE.nb=1;
 PDE.Ca2=[1 2];
 
 uinput.ic.ODE=[1 1 1 1];
 uinput.ic.PDE=[0 0];
 
  
%----------------------------------------
% Example 37
%----------------------------------------
 
    case 37
 
 p = 1; q = 0; qc = 3;

A = [0, -1/4, -1/5, 1/5, 1/6;
1/2, 1, -4, 9/2, 7/2;
-9/4, -1/2, -14, 23, 16;
-1/5, -1/2, -11/4, 1/10, 5/4;
-4/3, -4/3, -9, 9, 5/2];

B = [-7/2, -3/2, -1/10, 1/2, 1]';
C = [1/10, -1/3, -4, 7/8, 7/8];


 PDE.nw = 0;   PDE.nb = 0;   PDE.nz = 0;   PDE.no = 5;   PDE.nu = 0;
 PDE.n0 = 0;   PDE.n1 = 0;   PDE.n2 = 1;
 PDE.dom = [0,1];

 PDE.A2 = p;   PDE.A0 = qc-q;
 PDE.B22 = 0;    
 
 PDE.A = A; PDE.E0 = [B, zeros(5,3)];
 
 
 PDE.B = [0 0 1 0;
            0 1 0 0];
 PDE.Bx = [zeros(1,5);C];
 
 uinput.ic.ODE=[1 1 1 1 1];
 uinput.ic.PDE=[0];

%----------------------------------------
% Example 38 - example 2 from  Lhachemi & Prieur 2021 -
% Stable for both Dirichlet-Dirichlet and Dirichlet-Neumann
%----------------------------------------
 
    case 38

p = 1; q = 0; qc = 3;

A = [-1/4, -1/6, 2, 1, 1/12;
-3/2, -3/2, 5, 5, 1/6;
3/2, -4, -15/2, -5, -1/3;
-13/2, 22, 22, -14, -1/2;
1/7, -1/2, -1/2, 1/5, -5/2];

B = [-5/4, 2/3, 1/6, -1/6, 0]';
C = [-2/5, -5/4, 3/2, 1/3, 1/40];
    

 PDE.nw = 0;   PDE.nb = 0;   PDE.nz = 0;   PDE.no = 5;   PDE.nu = 0;
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

%----------------------------------------
% Cylindrical Coordinates: Variable Substitution Examples
%----------------------------------------

% Examples 39–40 use the substitution r = s^2 to remove the 1/r singularity
% in axisymmetric PDEs. After substitution, all terms become regular in s.

%----------------------------------------
% Example 39 - Axisymmetric Diffusion
% Neumann(left) and Dirichlet(right)
%----------------------------------------

% Original PDE in r: u_t = alpha[(1/r)*u_r + u_rr] 
% Original BCs:
%     r = 0:  u_r(0,t) = 0          (Neumann)
%     r = 1:  u(1,t)   = 0          (Dirichlet)

% Variable substitution:  r = s^2,  s ∈ [0,1] 
% Transformed PDE in s: u_t = 4*alpha*(u_s + s*u_ss) 
% Transformed BCs:
%     s = 1:  u(1,t) = 0
%     s = 0:  u_s(0,t) = -(j01^2 / 4) * exp(-j01^2 * t)
%             (this is the transformed Neumann condition at r=0)
% Exact solution: u(s,t) = J0(j01 * sqrt(s)) * exp(-j01^2 * t)
% where J0 is the Bessel function of the first kind (order 0) and j01 is
% its first zero.    

case 39
   
    a = 0.0001; b = 1;         
    PDE.dom = [a b];
    alpha = 1;  
    j01 = 2.4048;

    x1 = pde_var(s1, [a, b]);    
    x = x1;
    w = pde_var('in');
    
    Dyn = diff(x, t) == 4*alpha*(diff(x, s1) + s1 * diff(x, s1, 2));

    BCs = [subs(diff(x1, s1), s1, a) == w;   
               subs(x1, s1, b) == 0];                                    

    PDE = initialize([Dyn; BCs]);
    
%   Exact solution, initial conditions, and inhomogeneous inputs
    uinput.exact(1) = besselj(0, j01 * sqrt(sx)) * exp(-j01^2* st);
    uinput.ic.PDE = besselj(0, j01 * sqrt(sx));
    uinput.w(1) = -j01^2/4 * exp(-j01^2* st);

%----------------------------------------
% Example 40 - Axisymmetric Diffusion Reaction
% Neumann(left) and Dirichlet(right)
%----------------------------------------

% Original PDE in r: u_t = alpha[lambda*u + (1/r)*u_r + u_rr] 
% Original BCs:
%     r = 0:  u_r(0,t) = 0          (Neumann)
%     r = 1:  u(1,t)   = 0          (Dirichlet)

% Variable substitution:  r = s^2,  s ∈ [0,1] 
% Transformed PDE in s: u_t = lambda*u + 4*alpha*(u_s + s*u_ss) 
% Transformed BCs:
%     s = 1:  u(1,t) = 0
%     s = 0:  u_s(0,t) = -(j01^2 / 4) * exp(-(j01^2-lambda)*t)
%             (this is the transformed Neumann condition at r=0)
% Exact solution: u(s,t) = J0(j01 * sqrt(s)) * exp(-(j01^2-lambda)*t)
% where J0 is the Bessel function of the first kind (order 0) and j01 is
% its first zero.    

case 40

    a = 0.0001; b = 1;        
    PDE.dom = [a b];
    alpha = 1;
    lambda = 1;
    j01 = 2.4048;

    x1 = pde_var(s1, [a, b]);    
    x = x1;
    w = pde_var('in');
    
    Dyn = diff(x, t) == lambda * x + 4 * alpha * (diff(x, s1) + s1 * diff(x, s1, 2));

    BCs = [subs(diff(x1, s1), s1, a) == w;   
               subs(x1, s1, b) == 0];                                    

    PDE = initialize([Dyn; BCs]);
    
%   Exact solution, initial conditions, and inhomogeneous inputs
    uinput.exact(1) = besselj(0, j01 * sqrt(sx)) * exp(-(j01^2 - lambda) * st);
    uinput.ic.PDE = besselj(0, j01 * sqrt(sx));
    uinput.w(1) = -j01^2/4 * exp(-(j01^2 - lambda) * st);

%----------------------------------------
% Cylindrical Coordinates: Weighted Formulation Examples
%----------------------------------------

% Examples 41–43 illustrate the "weighted formulation" approach. 
% In cylindrical coordinates, many PDEs contain singular terms
% such as (1/r)*u_r or (1/r^2)*u. The weighted formulation avoids
% singularities by multiplying the entire PDE by the highest
% power of r that clears all singular terms.

%----------------------------------------
% Example 41 - Axisymmetric Diffusion
% Neumann(left) and Dirichlet(right)
%----------------------------------------

% Original PDE: u_t = alpha*( (1/r)u_r + u_rr)
% Weighted PDE: r*u_t = alpha*(u_r + r*u_rr) 
% BCs:
%     r = 0:  u_r(0,t) = 0          (Neumann)
%     r = 1:  u(1,t)   = 0          (Dirichlet)
% Exact solution: u(r,t) = J0(j01*r)*exp(-alpha*j01^2*t)
% Initial condition: u(r, 0) = J0(j01*r)

case 41

    a = 0.0001; b = 1;
    PDE.dom = [a b];
    alpha = 4;
    j01 = 2.4048;

    x1 = pde_var(s1, [a, b]);    
    x = x1;

    Dyn = diff(x, t) == alpha*(diff(x, s1) + s1 * diff(x, s1, 2));
    
    BCs = [subs(diff(x1, s1), s1, a) == 0;   
               subs(x1, s1, b) == 0];                                    

    PDE = initialize([Dyn; BCs]); 

    % Enable weighted formulation
    uinput.weighted = true;

    % Specify the polynomial weight. Since the PDE has the
    % highest singularity of order 1, we use weight = 1.
    uinput.weight = 1;

    uinput.exact = besselj(0, j01*sx)*exp(-alpha*j01^2* st);
    uinput.ic.PDE = besselj(0, j01*sx);

%----------------------------------------
% Example 42 - Axisymmetric Diffusion Reaction
% Neumann(left) and Dirichlet(right)
%----------------------------------------

% Original PDE: u_t = alpha*(lambda*u + (1/r)u_r + u_rr)
% PDE in r: r*u_t = alpha*(r*lambda*u + u_r + r*u_rr) 
% BCs:
%     r = 0:  u_r(0,t) = 0          (Neumann)
%     r = 1:  u(1,t)   = 0          (Dirichlet)
% Exact solution: u(r,t) = J0(j01*r)*exp(-alpha*(j01^2-lambda)*t) 
% Initial condition: u(r, 0) = J0(j01*r)

case 42

    a = 0.0001; b = 1;
    PDE.dom = [a b];
    alpha = 1;
    lambda = 1;
    j01 = 2.4048;

    x1 = pde_var(s1, [a, b]);    
    x = x1;

    Dyn = diff(x, t) == alpha*(diff(x, s1) + s1*diff(x, s1, 2) + s1*lambda*x);
    
    BCs = [subs(diff(x1, s1), s1, a) == 0;   
               subs(x1, s1, b) == 0];  

    PDE = initialize([Dyn; BCs]); 

    % Enable weighted formulation
    uinput.weighted = true;

    % Specify the polynomial weight. Since the PDE has the
    % highest singularity of order 1, we use weight = 1.
    uinput.weight = 1;
    
    uinput.exact = besselj(0, j01*sx)*exp(-alpha*(j01^2 - lambda)*st);
    uinput.ic.PDE = besselj(0, j01*sx);

%----------------------------------------
% Example 43 - Poloidal Magnetic Flux Gradient(Tokamak)
% This example is inspired by the model(PDE) presented in A.Gahlawat et al., 
% "Bootstrap Current Optimization in Tokamaks Using Sum-Of-Squares Polynomials"
%----------------------------------------

% Case I - Ramp-up Plasma Current I(t) = I0*(1 - exp(-t/tau))
% Original PDE: u_t = D(u_rr + (1/r)u_r - (1/r^2)u) + f(r,t)
% where f(r,t) = (I0/tau)*r*exp(-t/tau);
% Weighted PDE: r^2*u_t = D(r^2*u_rr + r*u_r - u) + r^2*f(r,t)
% BCs:
%     r = 0:  u(0,t) = 0          
%     r = 1:  u(1,t)   = I(t)
% IC: u(r,0) = r*I0
% Exact Solution: u(r, t) = r*I0*(1 - exp(-t/tau))

% IMPORTANT: In this case, since the highest singularity in the PDE is of the 
% order 1/r^2, we multiply the PDE throughout by r^2. To ensure the time 
% derivative is multiplied appropriately, we will need to modify the M 
% operator. To do this, uncomment the lines 140-157 and ensure
% that M.R.R0 = s1^2 * eye(m)

case 43

    a = 0; b = 1;
    PDE.dom = [a b];
    D = 0.04;
    I0 = 0.37;    % flat-top plasma current amplitude 
    tau = 15;     % ramp-up time constant 
    
    x1 = pde_var(s1, [a, b]);
    x = x1;
    w1 = pde_var('in', 1, s1, [a, b]);
    w2 = pde_var('in');
    
    Dyn = diff(x, t) == D*(s1^2*diff(x, s1, 2) + s1*diff(x, s1) - x) + s1^2*w1;
    
    BCs = [ subs(x1, s1, a) == 0;
            subs(x1, s1, b) == w2 ];
    
    PDE = initialize([Dyn; BCs]);

    % Enable weighted formulation
    uinput.weighted = true;

    % Specify the polynomial weight. Since the PDE has the
    % highest singularity of order 2, we use weight = 2.
    uinput.weight = 2;
    
    uinput.exact = sx*I0*(1 - exp(-st/tau));
    uinput.ic.PDE = sx*I0;
    uinput.w(1) = (I0/tau)*sx*exp(-st/tau); % Forcing term
    uinput.w(2) = I0 * (1 - exp(-st/tau));  % Boundary input
    
end % cases
% %            
