% date: 8/29/2019
close all; clear all; path(pathdef); clc;
addpath(genpath('.')) % makes sure any local version  of SOSTOOLS_MOD in the current folder is at the head of the path
stability=0; stability_dual=0; Hinf_gain=0; Hinf_gain_dual=0; Hinf_control=0; Hinf_estimator=0;Hinf_control_boundary=0;control_boundary=0;
pvar s theta;

% Solver file for coupled ODE-PDE systems
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This program determines stability, hinf-norm and designs hinf-optimal 
% observer of a linear coupled ODE-PDE which is defined in the format given below.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% System Definition:
% \dot [x0(t) ]  =A x0(t)+ E0 [x_2(a) + int(Ea(s) [x_1(t,s)] ds,a,b) + int(Eb(s) [x_2s(t,s)] ds,a,b)  + [B11   ]w(t) + [B12   ]u(t)
%                              x_2(b)             [x_2(t,s)]                     [x_3s(t,s)]            
%                              x_3(a)             [x_3(t,s)]
%                              x_3(b)
%                              x_3s(a)
%                              x_3s(b)]  


% \dot  [x(t,s)] = E(s) x0(t)+ A0(s) [x_1(t,s)] + A1(s) [x_2s(t,s)] + A2(s) [x_3ss(t,s)]+ [B21(s)]w(t)+ [B22(s)]u(t)
%                                    [x_2(t,s)]         [x_3s(t,s)]
%                                    [x_3(t,s)] 
% z(t) = C1 x0(t) + C10[x_2(a) + int(Ca1(s)[x_1(t,s)] ds,a,b) + int(Cb1(s)[x_2s(t,s)] ds,a,b) + [D11]w(t) + [D12]u(t)
%                       x_2(b)             [x_2(t,s)]                     [x_3s(t,s)]            
%                       x_3(a)             [x_3(t,s)]
%                       x_3(b)
%                       x_3s(a)
%                       x_3s(b)] 
%                        
% y(t) = C2 x0(t) + C20[x_2(a) + int(Ca2(s)[x_1(t,s)] ds,a,b) + int(Cb2(s)[x_2s(t,s)] ds,a,b) + [D21]w(t) + [D22]u(t)
%                       x_2(b)             [x_2(t,s)]                     [x_3s(t,s)]            
%                       x_3(a)             [x_3(t,s)]
%                       x_3(b)
%                       x_3s(a)
%                       x_3s(b)] 
%                                 
            
% % Boundary conditions are of the form 
% % B[x_2(a)       
% %   x_2(b)
% %   x_3(a)
% %   x_3(b)
% %   x_3s(a)
% %   x_3s(b)]= Bx0 x0(t)+ Bw w(t)+ Bu u(t)


% User Inputs:
% no     -  number of ODE states  
% n1     -  number of undifferentiated PDE states
% n2     -  number of single differentiated PDE states
% n3     -  number of double differentiated PDE states
% nu     -  number of inputs
% nw     -  number of disturbances
% ny     -  number of observed outputs
% nz     -  number of regulated outputs
% A0(s)  -  matrix function of s of dimension n1+n2+n3 x n1+n2+n3 
% A1(s)  -  matrix function of s of dimension n1+n2+n3 x n2+n3
% A2(s)  -  matrix function of s of dimension n1+n2+n3 x n3
% B      -  matrix of dimension n2+2*n3 x n1+n2+n3 (must be full row rank)
% Bx0     -  matrix of dimension n2+2*n3 x no (must be full row rank)
% Bw     -  matrix of dimension n2+2*n3 x nw (must be full row rank)
% Bu     -  matrix of dimension n2+2*n3 x nu (must be full row rank)
% B11    -  matrix of dimension no x nw 
% B12    -  matrix of dimension no x nu 
% B21(s) -  polynomial of dimension np x nw 
% B22(s) -  polynomial of dimension np x nu
% C1     -  matrix of dimension nz x no
% C10    -  matrix of dimension nz x 2n2+4n3 
% Ca1(s) -  polynomial of dimension nz x np
% Cb1(s) -  polynomial of dimension nz x n2+n3 
% D11    -  matrix of dimension nz x nw
% D12    -  matrix of dimension nz x nu
% C2     -  matrix of dimension ny x no
% C20    -  matrix of dimension ny x 2n2+4n3 
% Ca2(s) -  polynomial of dimension ny x np
% Cb2(s) -  polynomial of dimension ny x n2+n3 
% D21    -  matrix of dimension ny x nw 
% D22    -  matrix of dimension ny x nu 
% a,b    -  \R x \R - interval of the domain - s \in [a,b]
% eppos  -  for strict positivity
% epsneg (optional) - epsneg >0 for strictness of the negativity of derivative 




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% internal variables:
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


eppos = 1e-2;
eppos2 = 1e-2;
epneg = 1e-3; 

override1=1;
override2=1;
options1.sep = 0; %this is to select separable case, R1=R2 if sep=1
options1.exclude=[0 0 0 0]; 
options12.exclude = [0 0 0 0]; 
options12.psatz=1;

n_order1 = 1; 
n_order2 = 1; % This is supposed to be an accuracy/non-balancing degree
n_order3 = 1;
n_order4 = max([n_order2,n_order3]);

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Things you should probably not change
opts.psatz = 1;
opts.pure = 0;
options2.exclude= [0 1 0 0];
options3.exclude = [0 1 0 0]; 
options3.psatz=1;

Dup=6;


dd1 = {n_order1, [n_order2, n_order3, n_order4],[n_order2, n_order3, n_order4]}; 
dd12 = {n_order1, [n_order2, n_order3, n_order4],[n_order2, n_order3, n_order4]}; 
dd2 = {n_order1+1, [n_order2+Dup, n_order3+Dup, n_order4+Dup], [n_order2+Dup, n_order3+Dup, n_order4+Dup]};
dd3 = {n_order1, [n_order2+Dup-1, n_order3+Dup-1, n_order4+Dup-1], [n_order2+Dup-1, n_order3+Dup-1, n_order4+Dup-1]};
ddZ=[2*n_order1 2*n_order2 2*n_order3];


a=0;b=1;
X=[a,b];
nx1=0;nx2=2;

varlist = [s; theta];
prog = sosprogram(varlist);


disp('Parameterizing Positive Lyapunov Operator');
%nx1=10
%nx2=10

[prog, P1op] = sos_posopvar(prog, [nx1 ,nx2],X,s,theta,dd1,options1);

if override1~=1
    [prog, P2op] = sos_posopvar(prog, [nx1 ,nx2],X,s,theta,dd12,options12);
    Pop=P1op+P2op;
else
    Pop=P1op;
end
%toc

%sospos_RL2RL_MP_7_24_19
% enforce strict positivity on the operator
Pop.P = Pop.P+eppos*eye(nx1);
Pop.R.R0 = Pop.R.R0+eppos2*eye(nx2);  

%Assemble the big operator
% Pheq = [ -gamma*I  D'       B*PH
%          D         -gamma*I C
%          H*PB      C*       A*P*H+H'*PA]
% these Eop Aop are defined for the transport equation with periodic bc.
% x_t = x_s and x(a)=x(b)
opvar Eop Eop2 Aop Aop2;
% Aop2 is the PI form for A1(s)
% Eop2 refers to int(x(theta), theta, a, s) with a Q2 term for the x(a)
% boundary

Eop.R.R0 = eye(nx2);
Eop.R.R1 = -(1/(b-a))*eye(nx2);Eop.R.R2 = -1/(b-a)*eye(nx2);
Eop.I = X; 
Eop2.I = X; Eop2.R.R1 = eye(nx2);
Eop = Eop2*Eop; Eop.Q2 = eye(nx2);


Aop.R.R0=eye(nx2); Aop.R.R1=-1/(b-a)*eye(nx2); Aop.R.R2=-1/(b-a)*eye(nx2);
Aop.I =X;
Aop2.I = X; Aop2.R.R0 = [0 1;1 0]; Aop2.P = [];
Aop = Aop2*Aop;
disp('Constructing the derivative inequality');

Dop = [Eop'*(Pop*Aop)+Aop'*(Pop*Eop)]; 
    


%getting multiplier and kernel for derivative lyapunov function
disp('Parameterizing the negative operators');


% [prog, Deop] = sos_posopvar(prog, [nx1 ,nx2],X,s,theta,dd2,options2);
% [prog, De2op] = sos_posopvar(prog, [nx1 ,nx2],X,s,theta,dd3,options3);


% derivative negativity
% constraints
disp('Setting up the equality constraints');
tic
% prog = sosopeq(prog,Deop+De2op+Dop); %Dop=-Deop-De2op
prog = sosopineq(prog,-Dop,opts);
toc
% choosing a different solver if needed
% option.solver = 'mosek'; 

%solving the sos program
prog = sossolve(prog); 