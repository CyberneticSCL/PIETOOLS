%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% convert_PIETOOLS_PDE.m     PIETOOLS 2020a
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This a setup routine for converting PDEs into PIEs
% Converts the PDE formulation to PIE formulation
% This script uses a method of constructing a PIE
% representation of a PDE as defined in solver_PIETOOLS_PDE and
% error-checked in initialize_PIETOOLS_PDE.
% 
% A Partial Integral Equation is defined by 12 PI operators as
%
% TB1op \dot{w}(t)+TB2op \dot{u}(t)+Top \dot{x}(t)= Aop x(t) + B1op u(t)+ + B2op w(t)
%                                             z(t)= C1op x(t) + D11op u(t)+ + D12op w(t)
%                                             y(t)= C2op x(t) + D21op u(t)+ + D22op w(t)
%
% The formulae were defined in the paper
% ''A generalized LMI formulation for input-output analysis of linear systems of ODEs coupled with PDEs''
% reference: https://arxiv.org/abs/1904.10091

% convert_PIETOOLS_PDE file performs the following two tasks.
% 1) It verifies the dimension compatibility of input parameters of ODE-PDE
% and sets any missing parameters to zero.
% 2) It converts the input ODE-PDE representation to a PIE
% representation. 
% This script takes a user-defined PDE system in the format outlined in the
% header of solver_PIETOOLS_PDE and converts it to a PIE by defining the 11
% PI operators {BT1op,BT2op,Top,Aop,B1op,B2op,C1op,D11op,D12op,C2op,D21op,D22op} 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% The following script performs error checking operations and sets
% undefined operators to empty opbjects.
initialize_PIETOOLS_PDE

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Define auxiliary variables to convert the ODE-PDE to a PIE
X=[a b];
T = [eye(n1) zeros(n1,n2) zeros(n1,n2);
     eye(n1) zeros(n1,n2) zeros(n1,n2);
     zeros(n2,n1) eye(n2) zeros(n2);
     zeros(n2,n1) eye(n2) (b-a)*eye(n2);
     zeros(n2,n1) zeros(n2) eye(n2);
     zeros(n2,n1) zeros(n2) eye(n2)];
Q = [zeros(n1,n0) zeros(n1) zeros(n1,n2);
     zeros(n1,n0) eye(n1) zeros(n1,n2);
     zeros(n2,n0) zeros(n2,n1) zeros(n2);
     zeros(n2,n0) zeros(n2,n1) (b-theta)*eye(n2);
     zeros(n2,n0) zeros(n2,n1) zeros(n2);
     zeros(n2,n0) zeros(n2,n1) eye(n2)];
K = [zeros(n0,n1) zeros(n0,n2) zeros(n0,n2);
     eye(n1) zeros(n1,n2) zeros(n1,n2);
     zeros(n2,n1) eye(n2) (s-a)*eye(n2)];
L0 = [eye(n0) zeros(n0,n1) zeros(n0,n2);
     zeros(n1,n0) zeros(n1) zeros(n1,n2);
     zeros(n2,n0) zeros(n2,n1) zeros(n2)];
L1 = [zeros(n0) zeros(n0,n1) zeros(n0,n2);
     zeros(n1,n0) eye(n1) zeros(n1,n2);
     zeros(n2,n0) zeros(n2,n1) (s-theta)*eye(n2)];
V = [zeros(n1,n1) zeros(n1,n2) zeros(n1,n2);
     zeros(n2,n1) zeros(n2,n2) eye(n2)];
F0 = [zeros(n1,n0) eye(n1) zeros(n1,n2);
      zeros(n2,n0) zeros(n2,n1) zeros(n2)];
F1 = [zeros(n1,n0) zeros(n1) zeros(n1,n2);
      zeros(n2,n0) zeros(n2,n1) eye(n2)];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Test for well-posedness of Boundary conditions

if rcond(B*T)<1e-15
    error('Defined boundary conditions are rank deficient or have prohibited boundary conditions. Your system is likely ill-posed.');
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Converts ODE-PDE to PIE and defines PI operators
disp('Converting ODE-PDE to PIE');

BTinv=inv(B*T); % store repeatedly used variable

opvar H0 H1 Hbf TB1 TB2;  %initialize temporary opvars

% H0 maps [xo x0 x1s x2ss] to [xo x0 x1 x2] 
H0.dim = [nx nx; np np]; 
H0.var1 = s; H0.var2 = theta; 
H0.I = [a b]; 

H0.P = eye(nx); 
H0.Q1 = zeros(nx,np); 
H0.Q2 = K*BTinv*Bx;
H0.R.R0 = L0; 
H0.R.R1 = L1-K*BTinv*B*Q; 
H0.R.R2 = -K*BTinv*B*Q;

% TB1 maps w to [xo x0 x1 x2], TB2 maps u to [xo x0 x1 x2]
TB1.dim = [nx,nw;np,0]; 
TB1.var1 = s; TB1.var2 = theta; 
TB1.I = [a,b];
TB2.dim = [nx,nu;np,0]; 
TB2.var1 = s; TB2.var2 = theta; 
TB2.I = [a,b];

TB1.P = zeros(nx,nw); 
TB1.Q2 = K*BTinv*Bw;

TB2.P = zeros(nx,nu); 
TB2.Q2 = K*BTinv*Bu;

% H1 maps [xo x0 x1s x2ss] to [x1s x2s]
H1.dim = [0 nx; n1+n2 np]; 
H1.var1 = s; H1.var2 = theta; 
H1.I = [a b];

H1.Q2 = V*BTinv*Bx;
H1.R.R0 = F0; 
H1.R.R1 = F1 - V*BTinv*B*Q; 
H1.R.R2 = - V*BTinv*B*Q;

% Hbf maps [xo x0 x1s x2ss] to [x1(0) x1(L) x2(0) x2(L) x2s(0) x2s(L)]
Hbf.dim = [2*n1+4*n2 nx; 0 np]; 
Hbf.var1 = s; Hbf.var2 = theta; 
Hbf.I = [a b];

Hbf.P = T*BTinv*Bx; 
Hbf.Q1 = var_swap(-T*BTinv*B*Q+Q, s, theta); 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Convert parameters related to A's and E's to PI operator
opvar EAop EBop;

% EAop is integral operator on [x0 x1 x2] maps to dot(xo)
EAop.dim = [nx nx;0 np]; 
EAop.var1 = s; EAop.var2 = theta; 
EAop.I = [a b]; 

EAop.Q1 = Ea; 

%Ebop is integral operator on [x1s x2s] maps to dot(xo)
EBop.dim = [nx 0;0 n1+n2]; 
EBop.var1 = s; EBop.var2 = theta; 
EBop.I = [a b]; 

EBop.Q1 = Eb;

% Af1 maps [xo x0 x1s x2ss] to dot(xo)
Af1 = E0*Hbf+EAop*H0+EBop*H1;
Af1.P = Af1.P+A; 
Af1.Q1 = Af1.Q1+[zeros(nx,n0) zeros(nx,n1) Ec];


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
opvar A2p Ip; 
% A2p, temp var, maps x_2ss to dot(x0; x1s; x2ss)
A2p.dim = [0 nx; np np]; 
A2p.var1 = s; A2p.var2 = theta; 
A2p.I = [a b]; 

A2p.R.R0 = [zeros(np,n0) zeros(np,n1) A2]; 

%Ip is temp var, maps [xo x0 x1 x2] to [x0 x1 x2]
Ip.dim = [0 nx; np np]; 
Ip.var1 = s; Ip.var2 = theta; 
Ip.I = [a b]; 

Ip.R.R0 = eye(np); 

% Af2 maps [xo x0 x1s x2ss] to dot(x0; x1s; x2ss)
Af2 = A0*Ip*H0+A1*H1+A2p;
Af2.Q2 = Af2.Q2+E;

% Af maps [xo x0 x1s x2ss] to dot(xo; x0; x1s; x2ss)
% Af= [Af1; = [A     E0*Hbf+EA*H0+EB*H1;
%      Af2]    E     A0*Ip*H0+A1*H1+[zeros(np,n0) zeros(np,n1) A2]];
Af = [Af1; Af2]; 
clear Af1 Af2 EAop EBop A2p;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Convert all C's to PI operators
opvar CA1 CB1 CA2 CB2;
% CA1 and CA2 integral operators on [x0 x1 x2] maps to z/y resp
CA1.dim = [nz nx;0 np]; 
CA1.var1 = s; CA1.var2 = theta; 
CA1.I = [a b]; 

CA2.dim = [ny nx;0 np]; 
CA2.var1 = s; CA2.var2 = theta;
CA2.I = [a b];

CA1.Q1 = Ca1; 
CA2.Q1 = Ca2; 

% CB1 and CB2 integral operators on [x1s x2s] maps to z/y resp
CB1.dim = [nz 0;0 n1+n2]; 
CB1.var1 = s; CB1.var2 = theta; 
CB1.I = [a b]; 

CB2.dim = [ny 0;0 n1+n2]; 
CB2.var1 = s; CB2.var2 = theta; 
CB2.I = [a b];

CB1.Q1 = Cb1; 
CB2.Q1 = Cb2; 

%Assemble final C operator, Cf1 maps to [x; x0; x1s; x2ss] to z
% Cfi = [Ci Ci0*Hbf+CAi*H0+CBi*H1]
%       [0           0           ]
Cf1 = C10*Hbf + CA1*H0 + CB1*H1; 
Cf1.P = Cf1.P + C1; 
Cf1.Q1 = Cf1.Q1+[zeros(nz,n0) zeros(nz,n1) Cc1];

%Assemble final C operator, Cf2 maps to [x; x0; x1s; x2ss] to y
Cf2 = C20*Hbf + CA2*H0 + CB2*H1; 
Cf2.P = Cf2.P + C2; 
Cf2.Q1 = Cf2.Q1+[zeros(ny,n0) zeros(ny,n1) Cc2];

clear CA1 CB1 CA2 CB2;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Convert all B's and D's to PI operators
opvar Bf1 Bf2 Df11 Df12 Df21 Df22;

% Assemble B operator, Bfi's maps w/u to dot([xo;x0;x1;x2])
% Bfi = [Bf1.P;
%        Bf2.Q2]
Bf1.dim = [nx nw; np 0]; 
Bf1.var1 = s;Bf1.var2 = theta; 
Bf1.I = [a b];

Bf2.dim = [nx nu; np 0]; 
Bf2.var1 = s;Bf2.var2 = theta; 
Bf2.I = [a b];

Bf1.P = B11 +E0*T*BTinv*Bw+int(Ea*K*BTinv*Bw,s,a,b)+int(Eb*V*BTinv*Bw,s,a,b);
Bf1.Q2 = B21+A0*K*BTinv*Bw+A1*V*BTinv*Bw;

Bf2.P = B12+E0*T*BTinv*Bu+int(Ea*K*BTinv*Bu,s,a,b)+int(Eb*V*BTinv*Bu,s,a,b);
Bf2.Q2 = B22+A0*K*BTinv*Bu+A1*V*BTinv*Bu;

%Assemble D operators, Dfij's maps w/u to z/y
Df11.dim = [nz nw; 0 0]; 
Df11.var1 =s; Df11.var2 = theta;
Df11.I = [a b]; 

Df12.dim = [nz nu; 0 0]; 
Df12.var1 =s; Df12.var2 = theta;
Df12.I = [a b]; 

Df21.dim = [ny nw; 0 0]; 
Df21.var1 =s; Df21.var2 = theta;
Df21.I = [a b]; 

Df22.dim = [ny nu; 0 0]; 
Df22.var1 =s; Df11.var2 = theta;
Df22.I = [a b]; 

Df11.P = D11+C10*T*BTinv*Bw+int(Ca1*K*BTinv*Bw,s,a,b)+int(Cb1*V*BTinv*Bw,s,a,b); 
Df12.P = D12+C10*T*BTinv*Bu+int(Ca1*K*BTinv*Bu,s,a,b)+int(Cb1*V*BTinv*Bu,s,a,b); 
Df21.P = D21+C20*T*BTinv*Bw+int(Ca2*K*BTinv*Bw,s,a,b)+int(Cb2*V*BTinv*Bw,s,a,b); 
Df22.P = D22+C20*T*BTinv*Bu+int(Ca2*K*BTinv*Bu,s,a,b)+int(Cb2*V*BTinv*Bu,s,a,b); 


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Reassign using standard variable names
Top = H0;
TB1op = TB1; 
TB2op = TB2;
Aop = Af;
C1op = Cf1;
C2op = Cf2;
B1op = Bf1;
B2op = Bf2;
D12op = Df12;
D11op = Df11;
D21op = Df21;
D22op = Df22;
nx1 = nx; nx2 = np;

%remove temporary opvars
clear TB1 TB2 Af Cf1 Cf2 Bf1 Bf2 Df11 Df22 Df12 Df21; 