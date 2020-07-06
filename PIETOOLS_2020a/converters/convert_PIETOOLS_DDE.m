%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% convert_PIETOOLS_DDE.m     PIETOOLS 2020a
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Converts the DDE formulation to PIE formulation
% This script uses an inefficient method of constructing a PIE
% representation of a DDE as defined in solver_PIETOOLS_DDE and
% error-checked in initialize_PIETOOLS_DDE.
% 
% A Partial Integral Equation is defined by 12 PI operators as
%
% TB1op \dot{w}(t)+TB2op \dot{u}(t)+Top \dot{x}(t)= Aop x(t) + B1op u(t)+ B2op w(t)
%                                             z(t)= C1op x(t) + D11op u(t)+ D12op w(t)
%                                             y(t)= C2op x(t) + D21op u(t)+ D22op w(t)
%
% The formulae were defined in (12) and (15)-(16) in the paper
% ''Representation of Networks and Systems with Delay: DDEs, DDFs, ODE-PDEs and PIEs''
% reference: https://arxiv.org/abs/1910.03881
%
% The use of this script is not generally recommended as it exploits no
% structure of the system. Specifically, all states and inputs are delayed,
% forming a rather large-dimensional delay channel
%
% However, if there are no inputs or these are delayed, the representation
% will not be too bad.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Convert DDE to a 4-PI representation
% The following dimensions should have been defined in the
% initialize_PIETOOLS_DDE script
% nx=size(A0,1);   % number of states
% nw=size(B1,2);   % number of disturbances
% nu=size(B2,2);   % number of inputs
% nz=size(C1,1);  % number of regulated outputs
% ny=size(C2,1);  % number of sensed outputs

%For each delay in tau, the converter will create a delayed channel
%consisting of 
%
% ri(t)=[x(t);
%        w(t);
%        u(t)];
% The outputs of these channels are then picked up by
% Ai, Adi, C1i, C1di, C2i, C2di      for x(t) 
% B1i, B1di D11i, D11di, D21i, D21di for w(t)
% B2i, B2di D12i, D12di, D22i, D22di for u(t)
%


pvar theta
% First we define the Single-Pipe PIE Representation of the DDE
Itau=[];
Tb0=[];
Tb1=[];
Tb2=[];
for i=1:nK
    Itau=blkdiag(Itau,eye(nx+nw+nu)/tau(i)); % double check this
    Tb0=[Tb0;[eye(nx);zeros(nw,nx);zeros(nu,nx)]];
    Tb1=[Tb1;[zeros(nx,nw);eye(nw);zeros(nu,nw)]];
    Tb2=[Tb2;[zeros(nx,nu);zeros(nw,nu);eye(nu)]];
end
%oneK=repmat(eye(nx),nK,1);
Tba=zeros((nw+nu+nx)*nK);
Tbb=-eye((nw+nu+nx)*nK);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% We now define the matrices in Eqns 15-16 in ''Representation of Networks and Systems with Delay: DDEs, DDFs, ODE-PDEs and PIEs'' 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Ab0=A0;
Ab=[];
Bb1=B1;
Bb2=B2;
Cb10=C1;
Cb20=C2;
Cb11=[];
Cb21=[];
Db11=D11;
Db12=D12;
Db21=D21;
Db22=D22;
for i = 1:nK
  XAi=[Ai{i} B1i{i} B2i{i}]+tau(i)*int(subs([Adi{i} B1di{i} B2di{i}],s,tau(i)*s),s,-1,s);
  XC1i=[C1i{i} D11i{i} D12i{i}]+tau(i)*int(subs([C1di{i} D11di{i} D12di{i}],s,tau(i)*s),s,-1,s);
  XC2i=[C2i{i} D21i{i} D22i{i}]+tau(i)*int(subs([C2di{i} D21di{i} D22di{i}],s,tau(i)*s),s,-1,s);
  Ab0 = Ab0 + Ai{i} + int(tau(i)*subs(Adi{i},s,tau(i)*s),s,-1,0);
  Ab=[Ab -XAi];
  Bb1 = Bb1 + B1i{i} + int(tau(i)*subs(B1di{i},s,tau(i)*s),s,-1,0);
  Bb2 = Bb2 + B2i{i} + int(tau(i)*subs(B2di{i},s,tau(i)*s),s,-1,0);
  Db11 = Db11 + D11i{i} + int(tau(i)*subs(D11di{i},s,tau(i)*s),s,-1,0);
  Db12 = Db12 + D12i{i} + int(tau(i)*subs(D12di{i},s,tau(i)*s),s,-1,0);
  Db21 = Db21 + D21i{i} + int(tau(i)*subs(D21di{i},s,tau(i)*s),s,-1,0);
  Db22 = Db22 + D22i{i} + int(tau(i)*subs(D22di{i},s,tau(i)*s),s,-1,0);
  Cb11=[Cb11 -XC1i];
  Cb21=[Cb21 -XC2i];
  Cb10 = Cb10 + C1i{i} + int(tau(i)*subs(C1di{i},s,tau(i)*s),s,-1,0);
  Cb20 = Cb20 + C2i{i} + int(tau(i)*subs(C2di{i},s,tau(i)*s),s,-1,0);
  clear XAi XC1i XC2i
end
nxb=(nx+nw+nu)*nK; % number of infinite-dimensional states

% Initialize the system 4-PI operators
opvar Top Aop C1op C2op B1op B2op D11op D12op D21op D22op TB1op TB2op;
Top.I = [-1 0]; Top.dim = [nx nx; nxb nxb]; Top.var1 = s; Top.var2 = theta;
Aop.I = [-1 0]; Aop.dim = [nx nx; nxb nxb]; Aop.var1 = s; Aop.var2 = theta;
C1op.I = [-1 0]; C1op.dim = [nz nx; 0 nxb]; C1op.var1 = s; C1op.var2 = theta;
C2op.I = [-1 0]; C2op.dim = [ny nx; 0 nxb]; C2op.var1 = s; C2op.var2 = theta;
B1op.I = [-1 0]; B1op.dim=[nx nw; nxb 0]; B1op.var1 = s; B1op.var2 = theta;
B2op.I = [-1 0]; B2op.dim = [nx nu; nxb 0]; B2op.var1 = s; B2op.var2 = theta;
D11op.I = [-1 0]; D11op.dim = [nz nw; 0 0]; D11op.var1 = s; D11op.var2 = theta;
D12op.I = [-1 0]; D12op.dim = [nz nw; 0 0]; D12op.var1 = s; D12op.var2 = theta;
D21op.I = [-1 0]; D21op.dim = [nz nw; 0 0]; D21op.var1 = s; D21op.var2 = theta;
D22op.I = [-1 0]; D22op.dim = [nz nw; 0 0]; D22op.var1 = s; D22op.var2 = theta;
TB1op.I = [-1 0]; TB1op.dim = [nx nw; nxb 0]; TB1op.var1 = s; TB1op.var2 = theta;
TB2op.I = [-1 0]; TB2op.dim = [nx nu; nxb 0]; TB2op.var1 = s; TB2op.var2 = theta;

% We now define the 4-PI system operators
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Start with Top
Top.P=eye(nx);
Top.Q1=zeros(nx,nxb);
Top.Q2=Tb0;%oneK;
Top.R.R0=zeros(nxb);
Top.R.R1=Tba;%zeros(nxb);
Top.R.R2=Tbb;%-eye(nxb);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Next TB1op
TB1op.P=zeros(nx,nw);
TB1op.Q1=[];
TB1op.Q2=Tb1;
TB1op.R.R0=[];
TB1op.R.R1=[];
TB1op.R.R2=[];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Next TB2op
TB2op.P=zeros(nx,nu);
TB2op.Q1=[];
TB2op.Q2=Tb2;
TB2op.R.R0=[];
TB2op.R.R1=[];
TB2op.R.R2=[];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Next Aop
Aop.P=Ab0;
Aop.Q1=Ab;
Aop.Q2=zeros(nxb,nx);
Aop.R.R0=Itau;
Aop.R.R1=zeros(nxb);
Aop.R.R2=zeros(nxb);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Next C1op
C1op.P=Cb10;
C1op.Q1=Cb11;
C1op.Q2=[];
C1op.R.R0=[];
C1op.R.R1=[];
C1op.R.R2=[];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Next C2op
C2op.P=Cb20;
C2op.Q1=Cb21;
C2op.Q2=[];
C2op.R.R0=[];
C2op.R.R1=[];
C2op.R.R2=[];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Next B1op
B1op.P=Bb1;
B1op.Q1=[];
B1op.Q2=zeros(nxb,nw);
B1op.R.R0=[];
B1op.R.R1=[];
B1op.R.R2=[];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Next B2op
B2op.P=Bb2;
B2op.Q1=[];
B2op.Q2=zeros(nxb,nu);
B2op.R.R0=[];
B2op.R.R1=[];
B2op.R.R2=[];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Next D11op
D11op.P=Db11;
D11op.Q1=[];
D11op.Q2=[];
D11op.R.R0=[];
D11op.R.R1=[];
D11op.R.R2=[];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Next D12op
D12op.P=Db12;
D12op.Q1=[];
D12op.Q2=[];
D12op.R.R0=[];
D12op.R.R1=[];
D12op.R.R2=[];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Next D21op
D21op.P=Db21;
D21op.Q1=[];
D21op.Q2=[];
D21op.R.R0=[];
D21op.R.R1=[];
D21op.R.R2=[];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Next D22op
D22op.P=Db22;
D22op.Q1=[];
D22op.Q2=[];
D22op.R.R0=[];
D22op.R.R1=[];
D22op.R.R2=[];


