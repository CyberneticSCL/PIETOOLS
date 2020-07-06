%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% convert_PIETOOLS_DDF.m     PIETOOLS 2020a
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This a setup routine for converting DDFs into PIEs
% Converts the DDF formulation to PIE formulation
% This script uses a method of constructing a PIE
% representation of a DDF as defined in solver_PIETOOLS_DDF and
% error-checked in initialize_PIETOOLS_DDF.
% 
% A Partial Integral Equation is defined by 12 PI operators as
%
% TB1op \dot{w}(t)+TB2op \dot{u}(t)+Top \dot{x}(t)= Aop x(t) + B1op u(t)+ + B2op w(t)
%                                             z(t)= C1op x(t) + D11op u(t)+ + D12op w(t)
%                                             y(t)= C2op x(t) + D21op u(t)+ + D22op w(t)
%
% The formulae were defined in (12) and (15)-(16) in the paper
% ''Representation of Networks and Systems with Delay: DDEs, DDFs, ODE-PDEs and PIEs''
% reference: https://arxiv.org/abs/1910.03881



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Convert DDF to 4-PI representation
% nx=size(A0,1);   % number of states
% nw=size(B1,2);   % number of disturbances
% nu=size(B2,2);   % number of inputs
% nv=size(Cv{1},1);   % number of inputs
% nr{i}=size(Cr{1},1);   % number of inputs
% nz=size(C10,1);  % number of regulated outputs
% ny=size(C20,1);  % number of sensed outputs

% First we define the Single-Pipe PIE Representation of the TDS


%Chv{1}=Cv{1}+double(tau(1)*int(subs(Cvd{1},s,tau(1)*s),s,-1,0));

sumCvx=zeros(nv,nx);
sumDI=zeros(nv,nv);
sumDvw=zeros(nv,nw);
sumDvu=zeros(nv,nu);


for i = 1:nK
    Chv{i}=Cv{i}+double(tau(i)*int(subs(Cvd{i},s,tau(i)*s),s,-1,0));
    sumCvx=sumCvx+Chv{i}*Cr{i};
    sumDI=sumDI+Chv{i}*Drv{i};
    sumDvw=sumDvw+Chv{i}*Br1{i};
    sumDvu=sumDvu+Chv{i}*Br2{i};
end
DI=inv(eye(nv)-sumDI);
Cvx=DI*sumCvx;
Dvw=DI*sumDvw;
Dvu=DI*sumDvu;
for i = 1:nK
    CI{i}=-DI*(Cv{i}+tau(i)*int(subs(Cvd{i},s,theta*tau(i)),theta,-1,s));
end

%%%%%%%%%%%%%%%%%%%%%%%

Tb0=[];
Tb1=[];
Tb2=[];
for i = 1:nK
    Tb0=[Tb0; Cr{i}+Drv{i}*Cvx];
    Tb1=[Tb1; Br1{i}+Drv{i}*Dvw];
    Tb2=[Tb2; Br2{i}+Drv{i}*Dvu];
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Itau=[];
for i=1:nK
    Itau=blkdiag(Itau,eye(nr{i})/tau(i));
end
nxb=size(Itau,2); % number of infinite-dimensional states

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%
Tbai=[];
for i = 1:nK
    Tbai=[Tbai CI{i}];
end
Tbas=[];
for i = 1:nK
    Tbas=[Tbas;Drv{i}*Tbai];
end
Tba=subs(Tbas,s,theta);
Tbb=-eye(nxb)+Tba;

%%%%%%%%%%%%%%%%%%%%%%%%

Ab0=A0+Bv*Cvx;
Bb1=B1+Bv*Dvw;
Bb2=B2+Bv*Dvu;
Cb10=C1+D1v*Cvx;
Cb20=C2+D2v*Cvx;
Db11=D11+D1v*Dvw;
Db12=D12+D1v*Dvu;
Db21=D21+D2v*Dvw;
Db22=D22+D2v*Dvu;

% temp=[];
% for i=1:nK
% temp=[temp CI{i}];
% end
Ab=Bv*Tbai;
Cb11=D1v*Tbai;
Cb21=D2v*Tbai;

% Initialize the system 4-PI operators
opvar Top Aop C1op C2op B1op B2op D11op D12op D21op D22op TB1op TB2op;
Top.I = [-1 0]; Top.dim = [nx nx; nxb nxb]; Top.var1 = s; Top.var2 = theta;
Aop.I = [-1 0]; Aop.dim = [nx nx; nxb nxb]; Aop.var1 = s; Aop.var2 = theta;
C1op.I = [-1 0]; C1op.dim = [nz nx; 0 nxb]; C1op.var1 = s; C1op.var2 = theta;
C2op.I = [-1 0]; C2op.dim = [ny nx; 0 nxb]; C2op.var1 = s; C2op.var2 = theta;
B1op.I = [-1 0]; B1op.dim=[nx nw; nxb 0]; B1op.var1 = s; B1op.var2 = theta;
B2op.I = [-1 0]; B2op.dim = [nx nu; nxb 0]; B2op.var1 = s; B2op.var2 = theta;
D11op.I = [-1 0]; D11op.dim = [nz nw; 0 0]; D11op.var1 = s; D11op.var2 = theta;
D12op.I = [-1 0]; D12op.dim = [nz nu; 0 0]; D12op.var1 = s; D12op.var2 = theta;
D21op.I = [-1 0]; D21op.dim = [ny nw; 0 0]; D21op.var1 = s; D21op.var2 = theta;
D22op.I = [-1 0]; D22op.dim = [ny nw; 0 0]; D22op.var1 = s; D22op.var2 = theta;
TB1op.I = [-1 0]; TB1op.dim = [nx nw; nxb 0]; TB1op.var1 = s; TB1op.var2 = theta;
TB2op.I = [-1 0]; TB2op.dim = [nx nw; nxb 0]; TB2op.var1 = s; TB2op.var2 = theta;


% We now define the 12-opvar system representation
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Start with Aop
Aop.P=Ab0;
Aop.Q1=Ab;
Aop.Q2=zeros(nxb,nx);
Aop.R.R0=Itau;
Aop.R.R1=zeros(nxb);
Aop.R.R2=zeros(nxb);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Next Top
Top.P=eye(nx);
Top.Q1=zeros(nx,nxb);
Top.Q2=Tb0;
Top.R.R0=zeros(nxb);
Top.R.R1=Tba;
Top.R.R2=Tbb;

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

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
nx1=nx;nx2=nxb;
a=-1;b=0;
X = [a b];