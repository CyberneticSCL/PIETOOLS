%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% convert_PIETOOLS_PDE_MMP.m     PIETOOLS 2020a
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% convert_PIETOOLS_PDE_MMP is an alternative version of the PDE converter 
% file which performs the following two tasks.
% 1) It verifies the dimension compatibility of input parameters of ODE-PDE
% and sets any missing parameters to zero.
% 2) It converts the input ODE-PDE representation to a PIE
% representation. 
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
% The following script performs error checking operations and sets
% undefined operators to empty opbjects.
initialize_PIETOOLS_PDE


% Converts ODE-PDE to PIE and defines PI operators
disp('Converting ODE-PDE to PIE');


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Define auxiliary variables to convert the ODE-PDE to a PIE
X=[a b];
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
Exx=zeros(np,nrL1);
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
Apop.dim = [nz+ny+nx,nw+nu+nx+nrL1;np,nrL2]; Apop.var1 = s; Apop.var2 = theta; Apop.I = X;
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
     z20 z21 (b-s)*I2;
     z20 z21 z22;
     z20 z21 I2];
 
BTinv=inv(B*T);
Btemp=BTinv*[Bw Bu Bx];
Btemp2=BTinv*Bxx-B*Q;


K = [z01 z02 z02;
     I1 z12 z12;
     z21 I2 (s-a)*I2];

Q2f=K*Btemp;
R0f=[I0 z01 z02;
     z10 z11 z12;
     z20 z21 z22]; 
R2f=K*var_swap(Btemp2,s,theta);
R1f=R2f+[z00 z01 z02;
     z10 I1 z12;
     z20 z21 (s-theta)*I2];
 
 
 
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
Q2hf=Ib1*Q2f+diff(Ib2*Q2f+diff(Ib3*Q2f,s),s);
R0hf=Ib1*[I0 z01 z02;z10 z11 z12; z20 z21 z22]+Ib2*[z00 z01 z02;z10 I1 z12; z20 z21 z22]+Ib3*[z00 z01 z02;z10 z11 z12; z20 z21 I2];
R1hf=Ib1*R1f+diff(Ib2*R1f+diff(Ib3*R1f,s),s);
R2hf=Ib1*R2f+diff(Ib2*R2f+diff(Ib3*R2f,s),s);

%%% We now construct the Phfop
opvar Phfop;
Phfop.dim = [nw+nu+nx+nrL1, nw+nu+nx;nrL2, np]; Phfop.var1 = s; Phfop.var2 = theta; Phfop.I = X;
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
Tbigop.dim = [nx, nw+nu+nx;np, np]; Tbigop.var1 = s; Tbigop.var2 = theta; Tbigop.I = X;
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

% D11op_MMP=D11op;
% D21op_MMP=D21op;
% B1op_MMP=B1op;
% D12op_MMP=D12op;
% D22op_MMP=D22op;
% B2op_MMP=B2op;
% C1op_MMP=C1op;
% C2op_MMP=C2op;
% Aop_MMP=Aop;
% TB1op_MMP=TB1op;
% TB2op_MMP=TB2op;
% Top_MMP=Top;
% clear D11op D21op B1op D12op D22op B2op C1op C2op Aop TB1op TB2op Top;


nx1 = nx; nx2 = np;



%remove temporary opvars
clear Apop Phfop Tbigop; 