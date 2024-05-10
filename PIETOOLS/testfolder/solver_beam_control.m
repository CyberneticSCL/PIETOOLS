% date: 8/29/2019
close all; clear all; path(pathdef); clc;
addpath(genpath('.')) % makes sure any local version  of SOSTOOLS_MOD in the current folder is at the head of the path
stability=0; stability_dual=0; Hinf_gain=0; Hinf_gain_dual=0; Hinf_control=1; Hinf_estimator=0;Hinf_control_boundary=0;control_boundary=0;
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
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% NOTE: All cases were tested and found to work using the following optional 
% and degree choices(which are relatively high)
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%







%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%USER INPUT START%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Hinf gain %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% % Example 1, 0.0833  
% % output y(t) = int(x(s),a,b) ds
% 
% Hinf_gain=1
% % Hinf_gain_dual=1
n1=0;n2=4;n3=0;
nw =1; nz=1; nu =1; 
a=0;b=1;

A0=[0 0 0 0; 0 0 -1 0; 0 1 0 0; 0 0 0 0]; 
A1=[0 1 0 0; 1 0 0 0; 0 0 0 1; 0 0 1 0]; 
B21 = [0 0 1 0]';
B22=[0 0 1 0]';
% Bd1 = [0; 1];
C10 =[0 0 0 0 0 0 0 0]; Ca1 = [1,0,0,0]; Cb1 =zeros(1,4); D11 = 0;  % maybe u_s(L) or u_s(0) are not bounded by u_s is L2? - sachin
on = eye(n3); zer = zeros(n3);
B = [1 0 0 0 0 0 0 0;
     0 0 1 0 0 0 0 0;
     0 0 0 0 0 0 0 1;
     0 0 0 0 0 1 0 0];

setup_PIETOOLS_pde;
executive_PIETOOLS_Hinf_control;

% Conclusion:
if norm(prog.solinfo.info.feasratio-1)<=.1 && ~prog.solinfo.info.numerr
    disp('The System of equations was successfully solved.')
elseif norm(prog.solinfo.info.feasratio-1)<=.1 && prog.solinfo.info.numerr
    disp('The System  of equations was successfully solved. However, Double-check the precision.')
elseif prog.solinfo.info.pinf || prog.solinfo.info.dinf || norm(prog.solinfo.info.feasratio+1)<=.1
    disp('The System is  of equations not solved.')
else
    disp('Unable to definitively determine feasibility. Numerical errors dominating or at the limit of stability.')
end