%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% convert_PIETOOLS_DDF.m     PIETOOLS 2022
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function PIE=convert_PIETOOLS_DDF(DDF,out_type)
% This a setup routine for converting DDFs into PIEs
% Input argument "type" is optional, and should always be set to ''pie''.
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
% Initial coding MP - 10_01_2020
%  MP - 5_30_2021; converted script to function format
% SS - 9/28, added a dimension correction step to replace 0x0 empty
%  DJ - 09/01/2022: Output result as "pie_struct" object.
% matrices with empty matrices of correct size

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Convert DDF to 4-PI representation
if nargin==2 && ~strcmpi(out_type,'pie')
    error('Second argument must be ''pie'', or must be omitted. Only conversion from DDF to PIE is supported.')
end
% The following dimensions assume the DDE has been initialized, so lets
% make sure.
DDF=initialize_PIETOOLS_DDF(DDF);
 nK=length(DDF.tau);   % number of states
 nx=size(DDF.A0,1);   % number of states
 nw=size(DDF.B1,2);   % number of disturbances
 nu=size(DDF.B2,2);   % number of inputs
 nv=size(DDF.Cv{1},1);   % number of inputs
 for i=1:nK
    nr{i}=size(DDF.Cr{i},1);   % dimension of the original delayed channel, i
 end
 nz=size(DDF.C1,1);  % number of regulated outputs
 ny=size(DDF.C2,1);  % number of sensed outputs

 pvar s theta
% First we define the Single-Pipe PIE Representation of the TDS


%Chv{1}=Cv{1}+double(tau(1)*int(subs(Cvd{1},s,tau(1)*s),s,-1,0));

sumCvx=zeros(nv,nx);
sumDI=zeros(nv,nv);
sumDvw=zeros(nv,nw);
sumDvu=zeros(nv,nu);


for i = 1:nK
    Chv{i}=DDF.Cv{i}+double(DDF.tau(i)*int(subs(DDF.Cvd{i},s,DDF.tau(i)*s),s,-1,0));
    sumCvx=sumCvx+Chv{i}*DDF.Cr{i};
    sumDI=sumDI+Chv{i}*DDF.Drv{i};
    sumDvw=sumDvw+Chv{i}*DDF.Br1{i};
    sumDvu=sumDvu+Chv{i}*DDF.Br2{i};
end
DI=inv(eye(nv)-sumDI);
Cvx=DI*sumCvx;
Dvw=DI*sumDvw;
Dvu=DI*sumDvu;
for i = 1:nK
    CI{i}=-DI*(DDF.Cv{i}+DDF.tau(i)*int(subs(DDF.Cvd{i},s,theta*DDF.tau(i)),theta,-1,s));
end

%%%%%%%%%%%%%%%%%%%%%%%

Tb0=[];
Tb1=[];
Tb2=[];
for i = 1:nK
    Tb0=[Tb0; DDF.Cr{i}+DDF.Drv{i}*Cvx];
    Tb1=[Tb1; DDF.Br1{i}+DDF.Drv{i}*Dvw];
    Tb2=[Tb2; DDF.Br2{i}+DDF.Drv{i}*Dvu];
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Itau=[];
for i=1:nK
    Itau=blkdiag(Itau,eye(nr{i})/DDF.tau(i));
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
    Tbas=[Tbas;DDF.Drv{i}*Tbai];
end
Tba=subs(Tbas,s,theta);
Tbb=-eye(nxb)+Tba;

%%%%%%%%%%%%%%%%%%%%%%%%

Ab0=DDF.A0+DDF.Bv*Cvx;
Bb1=DDF.B1+DDF.Bv*Dvw;
Bb2=DDF.B2+DDF.Bv*Dvu;
Cb10=DDF.C1+DDF.D1v*Cvx;
Cb20=DDF.C2+DDF.D2v*Cvx;
Db11=DDF.D11+DDF.D1v*Dvw;
Db12=DDF.D12+DDF.D1v*Dvu;
Db21=DDF.D21+DDF.D2v*Dvw;
Db22=DDF.D22+DDF.D2v*Dvu;

% temp=[];
% for i=1:nK
% temp=[temp CI{i}];
% end
Ab=DDF.Bv*Tbai;
Cb11=DDF.D1v*Tbai;
Cb21=DDF.D2v*Tbai;

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


% correction of empty dimensions
Top.dim = [nx nx; nxb nxb]; 
Aop.dim = [nx nx; nxb nxb]; 
C1op.dim = [nz nx; 0 nxb]; 
C2op.dim = [ny nx; 0 nxb]; 
B1op.dim=[nx nw; nxb 0]; 
B2op.dim = [nx nu; nxb 0]; 
D11op.dim = [nz nw; 0 0]; 
D12op.dim = [nz nu; 0 0]; 
D21op.dim = [ny nw; 0 0]; 
D22op.dim = [ny nw; 0 0]; 
TB1op.dim = [nx nw; nxb 0]; 
TB2op.dim = [nx nw; nxb 0]; 


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
nx1=nx;nx2=nxb;
a=-1;b=0;
X = [a b];

PIE = pie_struct();
PIE.vars = [Top.var1,Top.var2];
PIE.dom=X;
PIE.T = Top; PIE.Tw = TB1op; PIE.Tu = TB2op;
PIE.A = Aop; PIE.B1 = B1op; PIE.B2 = B2op;
PIE.C1 = C1op; PIE.D11 = D11op; PIE.D12 = D12op;
PIE.C2 = C2op; PIE.D21 = D21op; PIE.D22 = D22op;


end
