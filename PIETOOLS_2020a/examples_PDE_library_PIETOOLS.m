%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% examples_PDE_library_PIETOOLS.m     PIETOOLS 2020a
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This library file contains the definition of some common ODE-PDE systems
% drawn from the literature. To use one of the examples, simply uncomment
% the variable definitions for that example. Please make sure no other
% example is uncommented, or you will have problems.
%
% The examples are grouped into various problem types. These types are as
% follows:
% 1. Stability Analysis Tests
% 2. Hinf-gain Analysis
% 3. Optimal Control Problems
% 4. Optimal Estimator Design Problems
%
% Within each type of examples, we have grouped the examples by the type of
% system they represent. Specifically, we have:
%
% 1. Hyperbolic Transport, Balance Laws, Conservation Equations
% 2. Diffusion and Heat-Equation Type Systems
% 3. Beam Type Equations
% 4. Wave Equations
%
% The examples are typically called using a line in the main PIETOOLS_PDE.m
% file. Simply save your changes to the library file, uncomment the line
% examples_PDE_library_PIETOOLS.m in PIETOOLS_PDE.m and run PIETOOLS_PDE.m.
% Of course, you can also run the library file directly.
%
% Most examples in each problem type include the option to solve the LPI
% according to their problem type. Hence, it is recommended that you also
% uncomment the toggle as well. e.g. stability=1, Hinf_gain=1, etc.
%
% When relevant, we also include citations for each example, indicating the
% sources. The bibtex for each citation is included at the end of the
% library file.
%
% If you wish to include a new example in our library, please send it to us
% and we will include it in the next release. Please also include the
% citation information, if available.
%
% NOTE: At present, PIETOOLS does not support inputs at the boundary for
% solving the Hinf optimal control problem. Support for this option will be
% included in a future release.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% STABILITY TEST Examples%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% % % % %------------------------------------------------------------------
% % % % % 1. Hyperbolic Transport, Balance Laws, Conservation Equations
% % % % %------------------------------------------------------------------
% % Transport equation u_t=u_{s}
% stability=1;
% stability_dual=1;
% 
% n0=0; n1=1; n2 =0;
% a=0; b= 1;
% 
% A0 = 0; A1 = -1; 
% 
% on = eye(n1); ze = zeros(n1);
% B = [on ze];
%
% % % %--------------------------------------------------------------------
%%%% Examples from Lamare [1]
%%%% There are several examples from [1] included here. Uncomment the
%%%% example you desire, then ALSO uncomment the problem definition at the
%%%% end
% % Example 5.1
% Gm1=[.2];
% Gm2=[-.3];
% Gm3=[.6];
% Gm4=[.1];
% Lm=[-3 0;0 1];
% Fm=[.2 -.3; .6 .1];
%
% % Example 5.2
% Gm1=[.1];
% Gm2=[-.8];
% Gm3=[.6];
% Gm4=[-.4];
% Lm=[-1 0;0 1];
% Fm=[-.3 .1; .1 -.3];
%
% % Example 5.3
% Gm1=[-.2202];
% Gm2=[1.3955];
% Gm3=[-.0596];
% Gm4=[.2090];
% Lm=[-1 0;0 2];
% Fm=[-.1 .1; .5 -.8];
%
% % Example 5.4
% Gm1=[.5];
% Gm2=[-.4];
% Gm3=[.2];
% Gm4=[.8];
% Lm=[-2 0;0 1];
% Fm=[-.6565 -.3743; -.113 -.6485];
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Uncomment for Examples from [2]
% stability=1
%%% stability_dual=1
% n0=0;n1=2;n2=0;
% a=0; b= 1;
% A0 = Fm; 
% A1= -Lm;
% A2=0*Fm;
% ny1=size(Gm1,1);ny2=size(Gm3,1);
% on0 = eye(ny1);on1 = eye(ny2); zer12 = zeros(ny1,ny2);
% B = [-Gm1 zer12 on0 -Gm2;
%      -Gm3 on1 zer12' -Gm4];
%
% % % %--------------------------------------------------------------------
% % Example from Diagne [2]
%
% stability=1
%%% stability_dual=1
%
% n0=0;n1=2;n2=0;
% a=0; b= 1;
%
% A0 = Mm; 
% A1= -Lm;
% A2=0*Fm;
%
% ny1=size(K00,1);ny2=size(K10,1);
% on0 = eye(ny1);on1 = eye(ny2); zer12 = zeros(ny1,ny2);
% B = [ on0 -K01 -K00 zer12;
%     zer12' -K11 -K10  on1;];
% 
% % % %--------------------------------------------------------------------
% % Example from Saba [3]
% stability=1
%%% stability_dual=1
%
% n0=0;n1=2;n2=0;
% a=0; b= 1;
%
% r1=.8;r2=1.1;sig1=2.3;sig2=-3.5;pb=.5;qb=-.7; % Stable with stripped settings
% r1=.5;r2=1.1;sig1=1;qb=1.2; sig2=-.1;pb=0; % Stable with stripped settings
% r1=.5;r2=1.1;sig1=1;qb=1.2; 
% sig2=1;pb=-.7; 
% sig2=.663;pb=0; % max sig2=.663 
% sig2=2.048;pb=-.4; % max sig2=1.049 
% A0 = [0 sig1; sig2 0]; 
% A1= [-1/r1 0; 0 1/r2];
%
%%%% BC: x(a)=qbx(b), pbx_s(a)=x_s(b) (this line should remain commented)
% B = [1 -qb 0 0;
%      0 0 -pb 1];
%
% % % % %------------------------------------------------------------------
% % % % % 2. Diffusion and Heat-Equation Type Systems
% % % % %------------------------------------------------------------------
% % Scalable Diffusion Equation on [0,1] adapted from Ahmadi 2015 [5]
% stability=1
%
% ne=1; lam=9.86;  
% n0=0;n1=0;n2=ne;
% a=0;b=1;
%
% A0=lam*eye(ne);
% A1=0*eye(ne);
% A2=1*eye(ne);
%
% % x(a)=0, x(b)=0   - Stable for lam<pi^2=9.8696 (this line should remain commented)
% B=[eye(ne) zeros(ne) zeros(ne) zeros(ne);
%     zeros(ne) eye(ne) zeros(ne) zeros(ne)];
%
% % % %--------------------------------------------------------------------
% % Scalable Diffusion Equation Problem 1 from Gahlawat_2017 [4]
%
% stability=1
% stability_dual=1
%
% ne=1; lam=2.466
% n0=0;n1=0;n2=ne;
% a=0;b=1;
%
% A0=lam*eye(ne);
% A1=0*eye(ne);
% A2=1*eye(ne);
%
% % x(a)=0, x_s(b)=0   - Unstable for \lam>2.467 (this line should remain commented)
% B=[eye(ne) zeros(ne) zeros(ne) zeros(ne);
%     zeros(ne)  zeros(ne) zeros(ne) eye(ne)];
%
% % % %--------------------------------------------------------------------
% % Diffusion Equation Problem 2 from Gahlawat_2017 [4]
% stability=1
% stability_dual=1
%
% ne=1; lam=4.66
% n0=0;n1=0;n2=ne;
% a=0;b=1;
%
% A2=s^3-s^2+2;
% A1=3*s^2-2*s;
% A0=-.5*s^3+1.3*s^2-1.5*s+.7+lam;
%
% % x(a)=0, x_s(b)=0   - Unstable for \lam>4.66 (this line should remain commented)
% B=[eye(ne) zeros(ne) zeros(ne) zeros(ne);
%     zeros(ne)  zeros(ne) zeros(ne) eye(ne)];
%
% % % %--------------------------------------------------------------------
% % Heat Equation coupled with unstable ODE at the boundary 
% stability=1
% stability_dual=1
%
% n0=0;n1=0;n2=1; no =1;
% a=0;b=1;
%
% A0=0; A1=0; A2=1; A = 1; E =0;
%
%%% % x(a)=0, x(b)=0   - Unstable always (this line should remain commented)
% B = [1 0 0 0; 0 1 0 0]; Bxo = [0;1];
%
% % % %--------------------------------------------------------------------
% % Example Diffusion-reaction from Ahmadi [6]
% stability=1
% stability_dual=0
%
% ne=2; R=2.93 % [2.93 2.94]
% n0 =0; n1=0; n2 =ne;
% a=0;b=1;
%
% A0=[1 1.5;
%     5 .2 ];
% A1=0*eye(ne);
% A2=inv(R)*eye(ne);
%
%%% % x(a)=0, x_s(b)=0   - Stable for R<2.7 (this line should remain commented)
% B=[eye(ne) zeros(ne) zeros(ne) zeros(ne);
%     zeros(ne) eye(ne) zeros(ne) zeros(ne) ];
%
% % % %--------------------------------------------------------------------
% % Diffusion-Reaction example from [5]
% stability=1
% stability_dual=0
%
% n0=0;n1=0;n2=3;
% a=0; b= 1;
% R = (21+9); 
%
% A0 = [0 0 0; s 0 0; s^2 -s^3 0]; A1 = zeros(np); 
% A2= (1/R)*eye(np);
%
%%% % x(a)=0, x_s(b)=0   - Stable for R<21 (this line should remain commented)
% B = [eye(n2) zeros(n2) zeros(n2) zeros(n2);
%      zeros(n2) eye(n2) zeros(n2) zeros(n2)];
%
% % % % %------------------------------------------------------------------
% % % % % 3. Beam Type Equations
% % % % %------------------------------------------------------------------
% % E-B beam equation [8]
% stability=1
% stability_dual=1
%
% c=.1%.01;%c=EI/mu
% n0=0;n1=0;n2=2;
% a=0;b=1;
%
% A0=[0 0; 0 0];
% A1=[0 0; 0 0];
% A2=[0 -c; 1 0];
%
% % x1(a)=0 x2(b)=0 x_1s(a)=0 x_2s(b)=0 (this line should remain commented)
% B=[ 1 0 0 0 0 0 0 0;
%     0 0 0 1 0 0 0 0;
%     0 0 0 0 1 0 0 0;
%     0 0 0 0 0 0 0 1];
%
% % % %--------------------------------------------------------------------
% % Timoschenko beam equation (hyperbolic) [8]
% % The states are u1=wt, u2=wx-phi, u3=phit, u4=phix
% stability=1
%
% k=1; aa=1;II=1; g=1; E=1;%5-30kPa
% r=1; ne=2;c=.01;%c=EI/mu
%
% n0=0;n1=4;n2=0; 
% a=0;b=1;
%
% A0=[0 0 0 0; 
%     0 0 -k*aa*g 0;
%     0 1/r/II 0 0;
%     0 0 0 0];
% A1=[0 1/r/aa 0 0; 
%     k*aa*g 0 0 0;
%     0 0 0 1/r/II;
%     0 0 E*II 0];
% A2=zeros(4,0);
%
% % x1(b)=0 x2(b)=0 x_1s(a)=0 x_2s(b)=0 (this line should remain commented)
% B=[ 1 0 0 0 0 0 0 0;
%     0 0 1 0 0 0 0 0;
%     0 0 0 0 0 0 0 1;
%     0 0 0 0 0 1 0 0];
%
%
% % % %--------------------------------------------------------------------
% % Timoschenko Beam equation Example (hyperbolic/diffusive) - unstable [8]
% % The states are u1=wt, u2=wx, u3=phit,u4=phi
% stability=1
%
% n0=0;n1=3;n2=1;
% a=0; b=1;
%
% A0=[0 0 0 0;
%     0 0 0 0;
%     0 1 0 -1;
%     0 0 1 0];
% A1=[0 1 0 -1;
%     1 0 0 0;
%     0 0 0 0;
%     0 0 0 0];
% A2=[0;0;1;0];
%
% %  u4(0)=0, u4x(L)=0, u3(0)=0, u1(0)=0, u2(L)-u4(L)=0 (this line should remain commented)
% B=[1 0 0 0 0 0 0 0 0 0;
%    0 0 1 0 0 0 0 0 0 0;
%    0 0 0 0 0 0 1 0 0 0;
%    0 0 0 0 0 0 0 0 0 1;
%    0 0 0 0 1 0 0 -1 0 0];
%
% % % % %------------------------------------------------------------------
% % % % % 4. Wave Equations
% % % % %------------------------------------------------------------------
% % % %--------------------------------------------------------------------
% % Boundary-Damped Wave equation (Hyperbolic) [8]
% stability=1
%
% n0=0;n1=2;n2=0;
%  a=0;b=1; k=1;
%
%  A0=[0 0; 0 0];
%  A1=[0 1; 1 0];
%  A2=zeros(2,0);
%
% % x(a)=0 kx_s(a)=-x_s(b)   (this line should remain commented)
%  B=[0 1 0 0;
%     0 0 k 1];
%
%
% % % %--------------------------------------------------------------------
% % Test 7.5c - Datko Boundary-Damped Wave equation [9]
% stability=1
%
% n0=0;n1=3;n2=0;
%  a=0;b=1; k=1;ad=1;
%
%  A0=[-2*ad -ad^2; 1 0];
%  A1=[0 0; 0 0];
%  A2=[1; 0];
%
% %  u_x(L)=-k u_t(L) (this line should remain commented)
%  B=[0 0 1 0 0 0;
%      1 0 0 0 0 0;
%      0 k 0 0 0 1];
%
% % % %--------------------------------------------------------------------
% % Test 7.5d - Datko Boundary-Damped Wave equation (Hyperbolic) [9]
% stability=1
%
% n0=0;n1=3;n2=0;
%  a=0;b=1; k=1;ad=1;
%
%  A0=[-2*ad -ad^2 0; 1 0 0;0 0 0];
%  A1=[0 0 1; 0 0 0;1 0 0];
%  A2=zeros(3,0);
%
% % x(0)=0 x(a)=K*x(L) x_s(a)=Kx_s(L)  (this line should remain commented)
%  B=[1 0 0 0 0 0;
%      0 1 0 0 0 0;
%      0 0 0 k 0 1];
%
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Hinf gain %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% % % % %------------------------------------------------------------------
% % % % % 1. Hyperbolic/Transport/Balance Type Systems
% % % % %------------------------------------------------------------------
% % Example pure transport equation 1D u_t + u_x=0
% Hinf_gain=1
% Hinf_gain_dual=1
%
% n0=0; n1=1; n2=0; nw=1; nz=1;
% a=0;b=1;
%
% A0 = zeros(np); A1=-ones(np,n1+n2); A2=zeros(np,n2);
% B21 = eye(nw); B21=1;
% C10 = zeros(nz,2*n1+4*n2); Ca1 = [1]; Cb1 =zeros(nz,n1+n2);
% D11 = zeros(nz,nw);
%
% on = eye(n1); zer = zeros(n1);
% B = [on zer];
% % gamma = 0.5 (this line should remain commented)
%
% % % %--------------------------------------------------------------------
% % Example tip damped wave equation u_tt= u_xx
% Hinf_gain=1
% Hinf_gain_dual=1
%
% n0=0; n1=2; n2=0; nw=1; nz=1;
% a=0; b=1;
%
% k =0.5; 
% A1 = [0 1; 1 0]; A0 = zeros(np); A2 = zeros(np,n2);
% B21 = [1; 0];
% C10 = zeros(nz,2*n1+4*n2); Ca1 = [1 0]; Cb1 =zeros(nz,n1+n2); D11 = zeros(nz,nw);
%
% on = eye(n2); zer = zeros(n2);
% B = [0 1 0 0; 0 0 k 1];
% % gamma = 2 for k = 0.5 (this line should remain commented)
%
% % % % %------------------------------------------------------------------
% % % % % 2. Diffusive/Heat Equation Type Systems
% % % % %------------------------------------------------------------------
% % Heat equation with distributed input u_t=u_ss
% Hinf_gain=1
%
% n0=0; n1=0; n2=1; nw=1; nz=1;
% a=0;b=1;
%
% A0=0; A1=0; A2=1; B21= s;
% C10 =[0 0 0 0]; Ca1 = 1; Cb1 =0; D11 = 0;  
%
% on = eye(n2); zer = zeros(n2);
% B = [on zer zer zer;
%      zer zer zer on];
% % gamma= 0.3333  (this line should remain commented)
%
% % % %--------------------------------------------------------------------
% %  Parabolic PDE example from [12]
% Hinf_gain=1
%%% Hinf_gain_dual=1
%
% n0=0; n1=0; n2=1; nw=1; nz=1;
% lamb = 4.6; 
% a=0; b= 1;
%
% A0 = -0.5*s^3+1.3*s^2-1.5*s+0.7+lamb;
% A1 = 3*s^2-2*s; A2 = s^3-s^2+2;
% B21 = s;
% C10 =[0,1,0,0]; Ca1 = 0; Cb1 =zeros(nz,n1+n2); D11 = zeros(nz,nw);
%
% on = eye(n2); zer = zeros(n2);
% B = [on zer zer zer;
%      zer zer zer on];
% % gamma =15.147 for lamb = 4.6 (this line should remain commented)
%
% % % %--------------------------------------------------------------------
% % Diffusion-Reaction equation from [12]
Hinf_gain=1;

n0=0; n1=0; n2=1; nw=1; nz=1;
lamb = (1-1e-2)*pi^2; 
a=0; b= 1;

A0=lamb; A1=0; A2=1;
B21 = 1;
C10 =[0 0 0 0]; Ca1 = 1; Cb =0; D11 = zeros(nz,nw);

B = [1 0 0 0;
     0 1 0 0];
% % gamma = 8.1069 for lamb = (1-1e-2)*pi^2 (this line should remain commented)

% % % %--------------------------------------------------------------------
% % Diffusion-reation pde from [6]
% Hinf_gain=1
%
% R=2.6; 
% n0=0; n1=0; n2=2; nw=1; nz=1;
% a=0; b= 1;
%
% A0 = [1 1.5;
%       5 0.2]; 
% A1= zeros(n1+n2);
% A2=(1/R)*[1 0; 0 1];
% B21 = [s;s];
% C10 =zeros(nz,2*n1+4*n2); Ca1 = [1 0]; Cb1 =zeros(nz,n1+n2); D11 = zeros(nz,nw);
%
% on = eye(n2); zer = zeros(n2);
% B = [on zer zer zer;
%      zer on zer zer];
% %  gamma = 0.8102 for R=2.6 (this line should remain commented)
% 
% % % %--------------------------------------------------------------------
% % Diffusion-reaction pde from [12]
% Hinf_gain=1
%
% n0=0; n1=0; n2=3; nw=1; nz=3;
% R = (21-1e-3); 
% a=0; b= 1;
%
% A0 = [0 0 0; s 0 0; s^2 -s^3 0]; A1 = zeros(n2); 
% A2= (1/R)*eye(n2);
% B21 = s*[1; 1; 1];
% C10 =zeros(nz,2*n1+4*n2); Ca1 = eye(n2); Cb1 =zeros(nz,n1+n2); D11 = zeros(nz,nw);
%
% on = eye(n2); zer = zeros(n2);
% B = [on zer zer zer;
%      zer on zer zer];
% % gamma =  4.23, for R = 21-1e-3 (this line should remain commented)
%
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Hinf optimal observer %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % % % %------------------------------------------------------------------
% % % % % 1. Hyperbolic Transport, Balance Laws, Conservation Equations
% % % % %------------------------------------------------------------------
% % % Example  - Answer: 1.0012
%
% nw = 1; ny = 1; nz = 1; no = 0;
% n0=0; n1 =1; n2 =0; np = n0+n1+n2; 
% 
% Hinf_estimator=1
% A0 = zeros(np); A1 = 1; A2 = zeros(n2);
% C10 = [0 0]; Ca1 = 1; Cb1 = zeros(nz,n1+n2);
% C20 = [0 0]; Ca2 = 1; Cb2 = zeros(ny,n1+n2);
% B21 = 1; D11 = 1; 
%   
% B=[0 1];
% a = 0; b =1;
% % answer 1.0012
%
% % % % %------------------------------------------------------------------
% % % % % 2. Diffusive/Heat Equation Type Systems
% % % % %------------------------------------------------------------------
% % % Example  - answer 1.0045 
% nw = 1; ny = 1; nz = 1; no = 0;
% n0=0; n1 =0; n2 =1; np = n0+n1+n2; 
% 
% Hinf_estimator=1
% A0 = zeros(np); A1 = 0; A2 = 1;
% C10 = [0 0 0 0]; Ca1 = 1; Cb1 = zeros(nz,n1+n2);
% C20 = [0 0 0 1]; Ca2 = zeros(ny,np); Cb2 = zeros(ny,n1+n2);
% B21 = 1; D11 = 1; 
% 
% B=[0 1 0 0; 1 0 0 0];
% a = 0; b =1;
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Hinf optimal controller %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % % % %------------------------------------------------------------------
% % % % % 2. Diffusive/Heat Equation Type Systems
% % % % %------------------------------------------------------------------
% Example 1, stabilizing controller, we set z=0, w=0 
% Hinf_control=1
% lamb = 10;
% nw = 0; ny = 0; nz = 0; no = 0; nu = 1;
% n0=0; n1 =0; n2 =1; np = n0+n1+n2; 
% A0 = lamb; A1 = zeros(np,n1+n2); A2 = 1;
% Ca1 = [0]; Cb1 = zeros(nz,n1+n2); %D12 = 1;
% B21 = 0; B22 = 1;
% on = eye(n2); zer = zeros(n2);
% 
% B=[on zer zer zer;
%    zer on zer zer];
% a = 0; b =1;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% References %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % [1] - 
% @article{lamare2016optimisation,
%   title={An optimisation approach for stability analysis and controller synthesis of linear hyperbolic systems},
%   author={Lamare, Pierre-Olivier and Girard, Antoine and Prieur, Christophe},
%   journal={ESAIM: Control, Optimisation and Calculus of Variations},
%   volume={22},
%   number={4},
%   pages={1236--1263},
%   year={2016},
%   publisher={EDP Sciences}
% }
% % [2] - 
% @article{diagne2012lyapunov,
%   title={Lyapunov exponential stability of 1-D linear hyperbolic systems of balance laws},
%   author={Diagne, Ababacar and Bastin, Georges and Coron, Jean-Michel},
%   journal={Automatica},
%   volume={48},
%   number={1},
%   pages={109--114},
%   year={2012},
%   publisher={Elsevier}
% }
% % [3] - 
% @article{saba2019stability,
%   title={Stability Analysis for a Class of Linear 2x2 Hyperbolic PDEs Using a Backstepping Transform},
%   author={Saba, David Bou and Argomedo, Federico Bribiesca and Auriol, Jean and Di Loreto, Michael and Di Meglio, Florent},
%   journal={IEEE Transactions on Automatic Control},
%   year={2019},
%   publisher={IEEE}
% }
% % [4] - 
% @article{gahlawat2016convex,
%   title={A convex sum-of-squares approach to analysis, state feedback and output feedback control of parabolic PDEs},
%   author={Gahlawat, Aditya and Peet, Matthew M},
%   journal={IEEE Transactions on Automatic Control},
%   volume={62},
%   number={4},
%   pages={1636--1651},
%   year={2016},
%   publisher={IEEE}
% }
% 
% % [5] - 
% @article{valmorbida2015stability,
%   title={Stability analysis for a class of partial differential equations via semidefinite programming},
%   author={Valmorbida, Giorgio and Ahmadi, Mohamadreza and Papachristodoulou, Antonis},
%   journal={IEEE Transactions on Automatic Control},
%   volume={61},
%   number={6},
%   pages={1649--1654},
%   year={2015},
%   publisher={IEEE}
% }
% 
% % [6] - 
% @inproceedings{valmorbida2014semi,
%   title={Semi-definite programming and functional inequalities for distributed parameter systems},
%   author={Valmorbida, Giorgio and Ahmadi, Mohamadreza and Papachristodoulou, Antonis},
%   booktitle={53rd IEEE conference on decision and control},
%   pages={4304--4309},
%   year={2014},
%   organization={IEEE}
% }
% % [8] - 
% @article{peet2019discussion,
%   title={Discussion paper: A new mathematical framework for representation and analysis of coupled pdes},
%   author={Peet, Matthew M and Shivakumar, Sachin and Das, Amritam and Weiland, Seip},
%   journal={IFAC-PapersOnLine},
%   volume={52},
%   number={2},
%   pages={132--137},
%   year={2019},
%   publisher={Elsevier}
% }
% 
% % [9] - 
% @article{datko1986example,
%   title={An example on the effect of time delays in boundary feedback stabilization of wave equations},
%   author={Datko, Richard and Lagnese, John and Polis, MP},
%   journal={SIAM Journal on Control and Optimization},
%   volume={24},
%   number={1},
%   pages={152--156},
%   year={1986},
%   publisher={SIAM}
% }
% % [12] -
% @inproceedings{shivakumar2019computing,
%   title={Computing input-ouput properties of coupled linear pde systems},
%   author={Shivakumar, Sachin and Peet, Matthew M},
%   booktitle={2019 American Control Conference (ACC)},
%   pages={606--613},
%   year={2019},
%   organization={IEEE}
% }
