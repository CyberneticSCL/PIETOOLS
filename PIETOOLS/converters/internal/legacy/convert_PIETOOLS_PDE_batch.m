%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% convert_PIETOOLS_PDE_batch.m     PIETOOLS 2022
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function PIE_out=convert_PIETOOLS_PDE_batch(PDE)
% convert_PIETOOLS_PDE_batch is an alternative version of the PDE converter 
% file which uses the batch input formula and performs the following two tasks.
% 1) It verifies the dimension compatibility of input parameters of ODE-PDE
% and sets any missing parameters to zero.
% 2) It converts the input ODE-PDE representation to a PIE
% representation. 
%
% For reference, the batch input format is listed at the end of this file.
%
% A Partial Integral Equation is defined by 11 PI operators as
%
% BT1op \dot{w}(t)+BT2op \dot{u}(t)+Top \dot{x}(t)= Aop x(t) + B1op u(t)+ + B2op w(t)
%                                             z(t)= C1op x(t) + D11op u(t)+ + D12op w(t)
%                                             y(t)= C2op x(t) + D21op u(t)+ + D22op w(t)
%
% This script takes a user-defined PDE system in the format outlined in the
% header of solver_PIETOOLS_PDE and converts it to a PIE by defining the 11
% PI operators {BT1op,BT2op,Top,Aop,B1op,B2op,C1op,D11op,D12op,C2op,D21op,D22op} 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Copyright (C)2022  M. Peet, S. Shivakumar, D. Jagt
%
% This program is free software; you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation; either version 2 of the License, or
% (at your option) any later version.
%
% This program is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
%
% You should have received a copy of the GNU General Public License
% along with this program; if not, write to the Free Software
% Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% If you modify this code, document all changes carefully and include date
% authorship, and a brief description of modifications
%
% Initial coding MMP  - 5_29_2019
%    MMP - updated to new data structure and made the script a function
% DJ, 09/29/2021: correction at the end to fill in empty parameters with
% appropriate zeros just in case...
% DJ, 12/07/2024: Use new default vars s1 and s1_dum;
% DJ, 12/16/2024: Output PIE as 'pie_struct' object;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% The following script performs error checking operations and sets
% undefined operators to empty opbjects.
PDE=initialize_PIETOOLS_PDE_batch(PDE);

n0=PDE.n0;

pvar s1 s1_dum                                                              % DJ, 12/07/2024

X=PDE.dom;
a=X(1); b=X(2);
nx=PDE.nx; % nx     -  number of ODE states
nw=PDE.nw; % nw     -  number of disturbances
nu=PDE.nu; % nu     -  number of controlled inputs
nz=PDE.nz; % nz     -  number of regulated outputs
ny=PDE.ny; % ny     -  number of observed outputs
%
% NOTE: Any matrix with a 0-dimension should be ommitted
%
% PDE Terms (required)
n0=PDE.n0;  % n0     -  number of undifferentiated PDE states (Required)
n1=PDE.n1;  % n1     -  number of single differentiated PDE states (Required)
n2=PDE.n2;  % n2     -  number of double differentiated PDE states (Required)
X=PDE.dom;  % dom    -  \R x \R - interval of the domain of the spatial variable - s \in [dom(1),dom(2)] (Required)
B=PDE.B;    % B      -  matrix of dimension n1+2*n2 x 2*(n1+2*n2) (Required, must have row rank n1+2n2 and contain no prohibited boundary conditions)
A0=PDE.A0;  % A0(s)  -  matrix function of s of dimension n0+n1+n2 x n0+n1+n2
A1=PDE.A1;  % A1(s)  -  matrix function of s of dimension n0+n1+n2 x n1+n2
A2=PDE.A2;  % A2(s)  -  matrix function of s of dimension n0+n1+n2 x n2
%%%%%% Input coupling (optional)
Bw=PDE.Bw;  % Bw     -  matrix of dimension n1+2*n2 x nw (must be full row rank)       - effect of disturbance on Boundary Conditions
Bu=PDE.Bu;  % Bu     -  matrix of dimension n1+2*n2 x nu (must be full row rank)       - effect of the input on the boundary conditions
B21=PDE.B21;% B21(s) -  polynomial matrix in pvar s of dimension (n0+n1+n2) x nw       - Distributed effect of disturbance in the domain of the PDE
B22=PDE.B22;% B22(s) -  polynomial matrix in pvar s of dimension (n0+n1+n2) x nu       - Distributed effect of input in the domain of the PDE
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ODE Dynamics (all optional)
nx=PDE.nx;  % nx     -  number of ODE states
A=PDE.A;    % A      - matrix of dimension nx x nx                                     - ODE dynamics
B11=PDE.B11;% B11    -  matrix of dimension nx x nw                                    - effect of disturbance on ODE state
B12=PDE.B12;% B12    -  matrix of dimension nx x nu                                    - effect of input on ODE state
%%%%%%%% Coupling of ODE to PDE
Bx=PDE.Bx;  % Bx    -  matrix of dimension n1+2*n2 x nx                                - Effect of ODE state on boundary conditions
E=PDE.E;    % E(s)   - pvar matrix in variable s of dimension n0+n1+n2 x nx            - Distributed effect of ODE in the domain of the PDE
%%%%%%%% Coupling of PDE to ODE
E0=PDE.E0;  % E0     -  matrix of dimension nx x 2*n1+4*n2                             - Effect of boundary values of distributed states on ODE states
Ea=PDE.Ea;  % Ea     -  polynomial matrix in pvar s of dimension nx x (n0+n1+n2)       - kernel of integration for effect of distributed states on ODE states
Eb=PDE.Eb;  % Eb     -  polynomial matrix in pvar s of dimension nx x n1+n2            - kernel of integration for effect of 1st-order derivatives of distributed states on ODE states
Ec=PDE.Ec;  % Ec     -  polynomial matrix in pvar s of dimension nx x n2               - kernel of integration for effect of 2nd-order derivatives of distributed states on ODE states
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% disturbance-related inputs (all optional)
% D21 (see observed outputs), 
% D11 (see regulated output)
% Bw,
% B21 (See PDE dynamics); 
% B11 (See ODE dynamics)
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% control-related inputs (all optional)
% D22 (see observed outputs); 
% D12 (see regulated output)
% Bu,
% B22 (See PDE dynamics); 
% B12 (See ODE dynamics)
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Regulated outputs (all optional)
% nz     -  number of regulated outputs
C1=PDE.C1;% C1     -  matrix of dimension nz x nx                                    - Effect of ODE states on regulated output
D11=PDE.D11;% D11    -  matrix of dimension nz x nw                                    - effect of disturbance on regulated output (avoid for H2-optimal control problems)
D12=PDE.D12;% D12    -  matrix of dimension nz x nu                                    - effect of control input on regulated output (Recommended for realistic controllers)
C10=PDE.C10;% C10    -  matrix of dimension nz x 2*n1+4*n2                             - Effect of boundary values of distributed states on regulated output
Ca1=PDE.Ca1;% Ca1(s) -  polynomial matrix in pvar s of dimension nz x (n0+n1+n2)       - kernel of integration for effect of distributed states on regulated outputs
Cb1=PDE.Cb1;% Cb1(s) -  polynomial matrix in pvar s of dimension nz x n1+n2            - kernel of integration for effect of 1st-order derivatives of distributed states on regulated outputs
Cc1=PDE.Cc1;% Cc1(s) -  polynomial matrix in pvar s of dimension nz x n2               - kernel of integration for effect of 2nd-order derivatives of distributed states on regulated outputs
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Observed outputs (all optional)
% ny     -  number of observed outputs
C2=PDE.C2;% C2     -  matrix of dimension ny x nx                                    - Effect of ODE states on observed output
D21=PDE.D21;% D21    -  matrix of dimension ny x nw                                    - effect of disturbance on observed output (e.g. sensor noise)
D22=PDE.D22;% D22    -  matrix of dimension ny x nu                                    - effect of control input on observed output (rare)
C20=PDE.C20;% C20    -  matrix of dimension ny x 2*n1+4*n2                             - Effect of boundary values of distributed states on observed outputs
Ca2=PDE.Ca2;% Ca2(s) -  polynomial matrix in pvar s of dimension ny x (n0+n1+n2)       - kernel of integration for effect of distributed states on observed outputs
Cb2=PDE.Cb2;% Cb2(s) -  polynomial matrix in pvar s of dimension ny x n1+n2            - kernel of integration for effect of 1st-order derivatives of distributed states on observed outputs
Cc2=PDE.Cc2;% Cc2(s)

% Converts ODE-PDE to PIE and defines PI operators
fprintf('\n --- Converting ODE-PDE to PIE --- \n');


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Define auxiliary variables to convert the ODE-PDE to a PIE
z00=zeros(n0);
z01=zeros(n0,n1);z10=z01';
z11=zeros(n1);
z02=zeros(n0,n2);z20=z02';
z12=zeros(n1,n2);z21=z12';
z22=zeros(n2);
I2=eye(n2);I1=eye(n1);I0=eye(n0);
zwu=zeros(nz,nu);zuw=zwu';
zwx=zeros(nz,nx);zxw=zwx';
Iw=eye(nw);

np=n0+n1+n2;
nrL1=2*n1+4*n2;
ncL1=np;
nrL2=n0+n1+n2+n1+2*n2;
ncL2=np;
 
%%%%%%%%%%%%%%%%%%%% Defining Primal Dynamics %%%%%%%%%%%%%%%%%%%%%%%
% Some terms are not yet compatible:
if isfield(PDE,'Exx')
Exx=PDE.Exx;
else
    Exx = zeros(np,nrL1);
end
Rxx1=zeros(np,nrL2);
Rxx2=zeros(np,nrL2);

Pb=[D11 D12 C1 C10;
    D21 D22 C2 C20;
    B11 B12 A  E0];
Q1b=[Ca1 Cb1 Cc1;
    Ca2 Cb2 Cc2;
    Ea Eb Ec];
Q2b=[B21 B22 E Exx];
 R0b=[A0 A1 A2];
 R1b=Rxx1; R2b=Rxx2;

%%% At this point, we can construct the primal dynamics opvar
opvar Apop;
Apop.dim = [nz+ny+nx,nw+nu+nx+nrL1;np,nrL2]; Apop.var1 = s1; Apop.var2 = s1_dum; Apop.I = X;
Apop.P=Pb;
Apop.Q1=Q1b;
Apop.Q2=Q2b;
Apop.R.R0=R0b;
Apop.R.R1=R1b;
Apop.R.R2=R2b;


%%%%%%%%%%%%%%%%%%%% Converting Lambda1 X  %%%%%%%%%%%%%%%%%%%%%%%
% Another Term not currently supported
Bxx=zeros(n1+2*n2,np);

%%% The following auxiliary matrices are used in this section 
T = [I1 z12 z12;
     I1 z12 z12;
     z21 I2 z22;
     z21 I2 (b-a)*I2;
     z21 z22 I2;
     z21 z22 I2];
 
 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Test for well-posedness of Boundary conditions

if rcond(B*T)<1e-15
    error('Defined boundary conditions are rank deficient or have prohibited boundary conditions. Your system is likely ill-posed.');
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Q = [z10 z11 z12;
     z10 I1 z12;
     z20 z21 z22;
     z20 z21 (b-s1)*I2;
     z20 z21 z22;
     z20 z21 I2];
 
BTinv=inv(B*T);
Btemp=BTinv*[Bw Bu Bx];
Btemp2=BTinv*(Bxx-B*Q);

K = [z01 z02 z02;
     I1 z12 z12;
     z21 I2 (s1-a)*I2];

Q2f=K*Btemp;
R0f=[I0 z01 z02;
     z10 z11 z12;
     z20 z21 z22]; 
R2f=K*var_swap(Btemp2,s1,s1_dum);
R1f=R2f+[z00 z01 z02;
     z10 I1 z12;
     z20 z21 (s1-s1_dum)*I2];
 
 
 
Ib1=[eye(np);   %nRL2(np+n1+2*n2) x np
    zeros(n1+2*n2,np)];
Ib2=[zeros(np);   %nRL2(np+n1+2*n2) x np
    z10 I1 z12;
    z20 z21 I2;
    z20 z21 z22];
Ib3=[zeros(np);   %nRL2(np+n1+2*n2) x np
    z10 z11 z12;
    z20 z21 z22;
    z20 z21 I2];

Phf=[eye(nw+nu+nx);
    T*Btemp];
Q1hf=[zeros(nw+nu+nx,np);
    Q+T*Btemp2] ;
Q2hf=Ib1*Q2f+diff(Ib2*Q2f+diff(Ib3*Q2f,s1),s1);
R0hf=Ib1*[I0 z01 z02;z10 z11 z12; z20 z21 z22]+Ib2*[z00 z01 z02;z10 I1 z12; z20 z21 z22]+Ib3*[z00 z01 z02;z10 z11 z12; z20 z21 I2];
R1hf=Ib1*R1f+diff(Ib2*R1f+diff(Ib3*R1f,s1),s1);
R2hf=Ib1*R2f+diff(Ib2*R2f+diff(Ib3*R2f,s1),s1);

%%% We now construct the Phfop
opvar Phfop;
Phfop.dim = [nw+nu+nx+nrL1, nw+nu+nx;nrL2, np]; Phfop.var1 = s1; Phfop.var2 = s1_dum; Phfop.I = X;
Phfop.P=Phf;
Phfop.Q1=Q1hf;
Phfop.Q2=Q2hf;
Phfop.R.R0=R0hf;
Phfop.R.R1=R1hf;
Phfop.R.R2=R2hf;

%%% We now construct the Ptop
Ptop=Apop*Phfop;

% %%% We now construct the Tbigop
% opvar Tbigop
% Tbigop.dim = [nz+ny+nx, nz+ny+nw+nu+nx;np, np]; Tbigop.var1 = s; Tbigop.var2 = theta; Tbigop.I = X;
% Tbigop.P=[eye(nz+ny) zeros(nz+ny,nw+nu+nx);
%     zeros(nx,nz+ny+nw+nu) eye(nx)];
% Tbigop.Q1=zeros(nz+ny+nx,np);
% Tbigop.Q2=[zeros(np,nz+ny) Q2f];
% Tbigop.R.R0=R0f;
% Tbigop.R.R1=R1f;
% Tbigop.R.R2=R2f;
%%% We now construct the Tbigop (smaller construction)

opvar Tbigop;
Tbigop.dim = [nx, nw+nu+nx;np, np]; Tbigop.var1 = s1; Tbigop.var2 = s1_dum; Tbigop.I = X;
Tbigop.P=[ zeros(nx,nw+nu) eye(nx)];
Tbigop.Q1=zeros(nx,np);
Tbigop.Q2=Q2f;
Tbigop.R.R0=R0f;
Tbigop.R.R1=R1f;
Tbigop.R.R2=R2f;


%%% Now we have to partition Ptop to get the desired pieces
D11op=op_slice(Ptop,1:nz,1:nw);
D21op=op_slice(Ptop,(nz+1):(nz+ny),1:nw);
B1op=op_slice(Ptop,(nz+ny+1):(nz+ny+nx+np),1:nw);
D12op=op_slice(Ptop,1:nz,(nw+1):(nw+nu));
D22op=op_slice(Ptop,(nz+1):(nz+ny),(nw+1):(nw+nu));
B2op=op_slice(Ptop,(nz+ny+1):(nz+ny+nx+np),(nw+1):(nw+nu));
C1op=op_slice(Ptop,1:nz,(nw+nu+1):(nw+nu+nx+np));
C2op=op_slice(Ptop,(nz+1):(nz+ny),(nw+nu+1):(nw+nu+nx+np));
Aop=op_slice(Ptop,(nz+ny+1):(nz+ny+nx+np),(nw+nu+1):(nw+nu+nx+np));
TB1op=op_slice(Tbigop,1:(nx+np),1:nw);
TB2op=op_slice(Tbigop,1:(nx+np),(nw+1):(nw+nu));
Top=op_slice(Tbigop,1:(nx+np),(nw+nu+1):(nw+nu+nx+np));

%nx1 = nx; nx2 = np;

% DJ 09/29/2021: fill in empty parameters with appropriate zeros
opnames = {'Top','TB1op','TB2op','Aop','B1op','B2op','C1op','D11op','D12op','C2op','D21op','D22op'};
for j=1:length(opnames)
    op = opnames{j};
    eval([op,'.dim=',op,'.dim;'])
end

%remove temporary opvars
%clear Apop Phfop Tbigop; 
pie_struct PIE_out;                                                         % DJ, 12/16/2024
PIE_out.dom = X;   PIE_out.vars = [Top.var1,Top.var2];
PIE_out.T = Top;   PIE_out.Tw = TB1op;  PIE_out.Tu = TB2op;
PIE_out.A = Aop;   PIE_out.B1 = B1op;   PIE_out.B2 = B2op;
PIE_out.C1 = C1op; PIE_out.D11 = D11op; PIE_out.D12 = D12op;
PIE_out.C2 = C2op; PIE_out.D21 = D21op; PIE_out.D22 = D22op;
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% System Definition:
% \dot [xo(t) ]  =A xo(t)+ E0 [x_2(a) + int(Ea(s) [x_1(t,s)] ds,a,b) + int(Eb(s) [x_2s(t,s)] ds,a,b) + int(Ec(s)[x_2ss(t,s)] ds,a,b) + [B11] w(t) + [B12] u(t)
%                              x_2(b)             [x_2(t,s)]                     [x_3s(t,s)]
%                              x_3(a)             [x_3(t,s)]
%                              x_3(b)
%                              x_3s(a)
%                              x_3s(b)]
%
% \dot [x(t,s)] = E(s) xo(t)+ A0(s) [x_1(t,s)] + A1(s) [x_2s(t,s)] + A2(s) [x_3ss(t,s)]+ [B21(s)]w(t)+ [B22(s)]u(t)
%                                   [x_2(t,s)]         [x_3s(t,s)]
%                                   [x_3(t,s)]
%
% z(t) = C1 xo(t) + C10[x_2(a) + int(Ca1(s)[x_1(t,s)] ds,a,b) + int(Cb1(s)[x_2s(t,s)] ds,a,b)+ int(Cc1(s)[x_2ss(t,s)] ds,a,b) + [D11]w(t) + [D12]u(t)
%                       x_2(b)             [x_2(t,s)]                     [x_3s(t,s)]
%                       x_3(a)             [x_3(t,s)]
%                       x_3(b)
%                       x_3s(a)
%                       x_3s(b)]
%
% y(t) = C2 xo(t) + C20[x_2(a) + int(Ca2(s)[x_1(t,s)] ds,a,b) + int(Cb2(s)[x_2s(t,s)] ds,a,b)+ int(Cc2(s)[x_2ss(t,s)] ds,a,b) + [D21]w(t) + [D22]u(t)
%                       x_2(b)             [x_2(t,s)]                     [x_3s(t,s)]
%                       x_3(a)             [x_3(t,s)]
%                       x_3(b)
%                       x_3s(a)
%                       x_3s(b)]
%
%
% % Boundary conditions are of the form
% % B[x_2(a)
% %   x_2(b)
% %   x_3(a)
% %   x_3(b)
% %   x_3s(a)
% %   x_3s(b)]= Bx xo(t)+ Bw w(t)+ Bu u(t)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% User-Defined Inputs for specifying the dynamics:
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% NOTE: All dimensions (nx,nw,nu,nz,ny) are assumed to be 0 by default.
% Specifically, (all optional)
% nx     -  number of ODE states
% nw     -  number of disturbances
% nu     -  number of controlled inputs
% nz     -  number of regulated outputs
% ny     -  number of observed outputs
%
% NOTE: Any matrix with a 0-dimension should be ommitted
%
% PDE Terms (required)
% n0     -  number of undifferentiated PDE states (Required)
% n1     -  number of single differentiated PDE states (Required)
% n2     -  number of double differentiated PDE states (Required)
% a,b    -  \R x \R - interval of the domain of the spatial variable - s \in [a,b] (Required)
% B      -  matrix of dimension n1+2*n2 x 2(n1+2n2) (Required, must have row rank n1+2n2 and contain no prohibited boundary conditions)
% A0(s)  -  matrix function of s of dimension n0+n1+n2 x n0+n1+n2
% A1(s)  -  matrix function of s of dimension n0+n1+n2 x n1+n2
% A2(s)  -  matrix function of s of dimension n0+n1+n2 x n2
%%%%%% Input coupling (optional)
% Bw     -  matrix of dimension n1+2*n2 x nw (must be full row rank)       - effect of disturbance on Boundary Conditions
% Bu     -  matrix of dimension n1+2*n2 x nu (must be full row rank)       - effect of the input on the boundary conditions
% B21(s) -  polynomial matrix in pvar s of dimension (n0+n1+n2) x nw       - Distributed effect of disturbance in the domain of the PDE
% B22(s) -  polynomial matrix in pvar s of dimension (n0+n1+n2) x nu       - Distributed effect of input in the domain of the PDE
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ODE Dynamics (all optional)
% no     -  number of ODE states
% A      - matrix of dimension nx x nx                                     - ODE dynamics
% B11    -  matrix of dimension nx x nw                                    - effect of disturbance on ODE state
% B12    -  matrix of dimension nx x nu                                    - effect of input on ODE state
%%%%%%%% Coupling of ODE to PDE
% Bx    -  matrix of dimension n1+2*n2 x nx                                - Effect of ODE state on boundary conditions
% E(s)   - pvar matrix in variable s of dimension n0+n1+n2 x nx            - Distributed effect of ODE in the domain of the PDE
%%%%%%%% Coupling of PDE to ODE
% E0     -  matrix of dimension nx x 2*n1+4*n2                             - Effect of boundary values of distributed states on ODE states
% Ea     -  polynomial matrix in pvar s of dimension nx x (n0+n1+n2)       - kernel of integration for effect of distributed states on ODE states
% Eb     -  polynomial matrix in pvar s of dimension nx x n1+n2            - kernel of integration for effect of 1st-order derivatives of distributed states on ODE states
% Ec     -  polynomial matrix in pvar s of dimension nx x n2               - kernel of integration for effect of 2nd-order derivatives of distributed states on ODE states
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% disturbance-related inputs (all optional)
% D21 (see observed outputs), D11 (see regulated output)
% Bw,B21 (See PDE dynamics); B11 (See ODE dynamics)
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% control-related inputs (all optional)
% D22 (see observed outputs); D12 (see regulated output)
% Bu,B22 (See PDE dynamics); B12 (See ODE dynamics)
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Regulated outputs (all optional)
% nz     -  number of regulated outputs
% C1     -  matrix of dimension nz x nx                                    - Effect of ODE states on regulated output
% D11    -  matrix of dimension nz x nw                                    - effect of disturbance on regulated output (avoid for H2-optimal control problems)
% D12    -  matrix of dimension nz x nu                                    - effect of control input on regulated output (Recommended for realistic controllers)
% C10    -  matrix of dimension nz x 2*n1+4*n2                             - Effect of boundary values of distributed states on regulated output
% Ca1(s) -  polynomial matrix in pvar s of dimension nz x (n0+n1+n2)       - kernel of integration for effect of distributed states on regulated outputs
% Cb1(s) -  polynomial matrix in pvar s of dimension nz x n1+n2            - kernel of integration for effect of 1st-order derivatives of distributed states on regulated outputs
% Cc1(s) -  polynomial matrix in pvar s of dimension nz x n2               - kernel of integration for effect of 2nd-order derivatives of distributed states on regulated outputs
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Observed outputs (all optional)
% ny     -  number of observed outputs
% C2     -  matrix of dimension ny x nx                                    - Effect of ODE states on observed output
% D21    -  matrix of dimension ny x nw                                    - effect of disturbance on observed output (e.g. sensor noise)
% D22    -  matrix of dimension ny x nu                                    - effect of control input on observed output (rare)
% C20    -  matrix of dimension ny x 2*n1+4*n2                             - Effect of boundary values of distributed states on observed outputs
% Ca2(s) -  polynomial matrix in pvar s of dimension ny x (n0+n1+n2)       - kernel of integration for effect of distributed states on observed outputs
% Cb2(s) -  polynomial matrix in pvar s of dimension ny x n1+n2            - kernel of integration for effect of 1st-order derivatives of distributed states on observed outputs
% Cc2(s) -  polynomial matrix in pvar s of dimension ny x n2               - kernel of integration for effect of 2nd-order derivatives of distributed states on observed outputs
%