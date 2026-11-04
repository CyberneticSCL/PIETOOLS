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
% matrices with empty matrices of correct size
%  DJ - 09/01/2022: Output result as "pie_struct" object.
% DJ, 12/07/2024: Use new default vars (s1,s2) and (s1_dum,s2_dum);
% DJ, 10/29/2026: Set the spatial variables and domain in PIE.space, instead
%                   of the dependent fields PIE.vars and PIE.dom;
% MMP, 11/04/2026: Build the twelve operators as 'copvar' containers of
%                   'sopvar' blocks (PI4_COPVAR_SOP) instead of 'opvar'
%                   objects that INITIALIZE converted with OP2COPVAR_SOP.
%                   PI4_COPVAR_SOP builds the blocks OPVAR2COPVAR builds
%                   from the same fields, so the PIE is unchanged; the DDE,
%                   NDS and DDF conversion no longer constructs an 'opvar',
%                   which is being removed from PIETOOLS. PIE.space.vars
%                   from s1, as T no longer has the field var1.

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

 pvar s1 s1_dum                                                             % DJ, 12/07/2024
% First we define the Single-Pipe PIE Representation of the TDS


%Chv{1}=Cv{1}+double(tau(1)*int(subs(Cvd{1},s,tau(1)*s),s,-1,0));

sumCvx=zeros(nv,nx);
sumDI=zeros(nv,nv);
sumDvw=zeros(nv,nw);
sumDvu=zeros(nv,nu);


for i = 1:nK
    Chv{i}=DDF.Cv{i}+double(DDF.tau(i)*int(subs(DDF.Cvd{i},s1,DDF.tau(i)*s1),s1,-1,0));
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
    CI{i}=-DI*(DDF.Cv{i}+DDF.tau(i)*int(subs(DDF.Cvd{i},s1,s1_dum*DDF.tau(i)),s1_dum,-1,s1));
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
Tba=subs(Tbas,s1,s1_dum);
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

% % % BEGIN change MMP, 11/04/2026: the operators as containers.            % MMP, 11/04/2026
% Each operator [P, Q1; Q2, R]: R^p0 x L2^p1[s1] -> R^q0 x L2^q1[s1], R the % MMP, 11/04/2026
% multiplier R0 plus the integrals of R1 on [-1,s1] and R2 on [s1,0], is    % MMP, 11/04/2026
% built by PI4_COPVAR_SOP from the components of the former opvar, block    % MMP, 11/04/2026
% for block as OP2COPVAR_SOP converted that opvar in INITIALIZE: the PIE    % MMP, 11/04/2026
% is unchanged, and no opvar is constructed.                                % MMP, 11/04/2026
% Arguments: P, Q1, Q2, {R0,R1,R2}, dim = [q0 p0; q1 p1] as the opvar dim.  % MMP, 11/04/2026
% Former code, kept as comments after the new: the opvar declarations and   % MMP, 11/04/2026
% fields (MP 10_01_2020; variables DJ 12/07/2024), dims (SS 9/28).          % MMP, 11/04/2026
I1 = [-1 0];                       % domain of s1, the former opvar field I % MMP, 11/04/2026
Top   = pi4_copvar_sop(eye(nx),zeros(nx,nxb),Tb0,{zeros(nxb),Tba,Tbb},...
                       [nx nx; nxb nxb],s1,s1_dum,I1);                      % MMP, 11/04/2026
Aop   = pi4_copvar_sop(Ab0,Ab,zeros(nxb,nx),{Itau,zeros(nxb),zeros(nxb)},...
                       [nx nx; nxb nxb],s1,s1_dum,I1);                      % MMP, 11/04/2026
B1op  = pi4_copvar_sop(Bb1,[],zeros(nxb,nw),{[],[],[]},...
                       [nx nw; nxb 0],s1,s1_dum,I1);                        % MMP, 11/04/2026
B2op  = pi4_copvar_sop(Bb2,[],zeros(nxb,nu),{[],[],[]},...
                       [nx nu; nxb 0],s1,s1_dum,I1);                        % MMP, 11/04/2026
TB1op = pi4_copvar_sop(zeros(nx,nw),[],Tb1,{[],[],[]},...
                       [nx nw; nxb 0],s1,s1_dum,I1);                        % MMP, 11/04/2026
TB2op = pi4_copvar_sop(zeros(nx,nu),[],Tb2,{[],[],[]},...
                       [nx nu; nxb 0],s1,s1_dum,I1);                        % MMP, 11/04/2026
C1op  = pi4_copvar_sop(Cb10,Cb11,[],{[],[],[]},...
                       [nz nx; 0 nxb],s1,s1_dum,I1);                        % MMP, 11/04/2026
C2op  = pi4_copvar_sop(Cb20,Cb21,[],{[],[],[]},...
                       [ny nx; 0 nxb],s1,s1_dum,I1);                        % MMP, 11/04/2026
D11op = pi4_copvar_sop(Db11,[],[],{[],[],[]},[nz nw; 0 0],s1,s1_dum,I1);    % MMP, 11/04/2026
D12op = pi4_copvar_sop(Db12,[],[],{[],[],[]},[nz nu; 0 0],s1,s1_dum,I1);    % MMP, 11/04/2026
D21op = pi4_copvar_sop(Db21,[],[],{[],[],[]},[ny nw; 0 0],s1,s1_dum,I1);    % MMP, 11/04/2026
D22op = pi4_copvar_sop(Db22,[],[],{[],[],[]},[ny nu; 0 0],s1,s1_dum,I1);    % MMP, 11/04/2026
% Former code:                                                              % MMP, 11/04/2026
% % Initialize the system 4-PI operators                                    % MMP, 11/04/2026 (was)
% opvar Top Aop C1op C2op B1op B2op D11op D12op D21op D22op TB1op TB2op;    % MMP, 11/04/2026 (was)
% Top.I = [-1 0]; Top.dim = [nx nx; nxb nxb]; Top.var1 = s1; Top.var2 = s1_dum; % MMP, 11/04/2026 (was)
% Aop.I = [-1 0]; Aop.dim = [nx nx; nxb nxb]; Aop.var1 = s1; Aop.var2 = s1_dum; % MMP, 11/04/2026 (was)
% C1op.I = [-1 0]; C1op.dim = [nz nx; 0 nxb]; C1op.var1 = s1; C1op.var2 = s1_dum; % MMP, 11/04/2026 (was)
% C2op.I = [-1 0]; C2op.dim = [ny nx; 0 nxb]; C2op.var1 = s1; C2op.var2 = s1_dum; % MMP, 11/04/2026 (was)
% B1op.I = [-1 0]; B1op.dim=[nx nw; nxb 0]; B1op.var1 = s1; B1op.var2 = s1_dum; % MMP, 11/04/2026 (was)
% B2op.I = [-1 0]; B2op.dim = [nx nu; nxb 0]; B2op.var1 = s1; B2op.var2 = s1_dum; % MMP, 11/04/2026 (was)
% D11op.I = [-1 0]; D11op.dim = [nz nw; 0 0]; D11op.var1 = s1; D11op.var2 = s1_dum; % MMP, 11/04/2026 (was)
% D12op.I = [-1 0]; D12op.dim = [nz nu; 0 0]; D12op.var1 = s1; D12op.var2 = s1_dum; % MMP, 11/04/2026 (was)
% D21op.I = [-1 0]; D21op.dim = [ny nw; 0 0]; D21op.var1 = s1; D21op.var2 = s1_dum; % MMP, 11/04/2026 (was)
% D22op.I = [-1 0]; D22op.dim = [ny nu; 0 0]; D22op.var1 = s1; D22op.var2 = s1_dum; % MMP, 11/04/2026 (was)
% TB1op.I = [-1 0]; TB1op.dim = [nx nw; nxb 0]; TB1op.var1 = s1; TB1op.var2 = s1_dum; % MMP, 11/04/2026 (was)
% TB2op.I = [-1 0]; TB2op.dim = [nx nu; nxb 0]; TB2op.var1 = s1; TB2op.var2 = s1_dum; % MMP, 11/04/2026 (was)


% % We now define the 12-opvar system representation                        % MMP, 11/04/2026 (was)
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% % MMP, 11/04/2026 (was)
% % Start with Aop                                                          % MMP, 11/04/2026 (was)
% Aop.P=Ab0;                                                                % MMP, 11/04/2026 (was)
% Aop.Q1=Ab;                                                                % MMP, 11/04/2026 (was)
% Aop.Q2=zeros(nxb,nx);                                                     % MMP, 11/04/2026 (was)
% Aop.R.R0=Itau;                                                            % MMP, 11/04/2026 (was)
% Aop.R.R1=zeros(nxb);                                                      % MMP, 11/04/2026 (was)
% Aop.R.R2=zeros(nxb);                                                      % MMP, 11/04/2026 (was)
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% % MMP, 11/04/2026 (was)
% % Next Top                                                                % MMP, 11/04/2026 (was)
% Top.P=eye(nx);                                                            % MMP, 11/04/2026 (was)
% Top.Q1=zeros(nx,nxb);                                                     % MMP, 11/04/2026 (was)
% Top.Q2=Tb0;                                                               % MMP, 11/04/2026 (was)
% Top.R.R0=zeros(nxb);                                                      % MMP, 11/04/2026 (was)
% Top.R.R1=Tba;                                                             % MMP, 11/04/2026 (was)
% Top.R.R2=Tbb;                                                             % MMP, 11/04/2026 (was)

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% % MMP, 11/04/2026 (was)
% % Next B1op                                                               % MMP, 11/04/2026 (was)
% B1op.P=Bb1;                                                               % MMP, 11/04/2026 (was)
% B1op.Q1=[];                                                               % MMP, 11/04/2026 (was)
% B1op.Q2=zeros(nxb,nw);                                                    % MMP, 11/04/2026 (was)
% B1op.R.R0=[];                                                             % MMP, 11/04/2026 (was)
% B1op.R.R1=[];                                                             % MMP, 11/04/2026 (was)
% B1op.R.R2=[];                                                             % MMP, 11/04/2026 (was)

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% % MMP, 11/04/2026 (was)
% % Next B2op                                                               % MMP, 11/04/2026 (was)
% B2op.P=Bb2;                                                               % MMP, 11/04/2026 (was)
% B2op.Q1=[];                                                               % MMP, 11/04/2026 (was)
% B2op.Q2=zeros(nxb,nu);                                                    % MMP, 11/04/2026 (was)
% B2op.R.R0=[];                                                             % MMP, 11/04/2026 (was)
% B2op.R.R1=[];                                                             % MMP, 11/04/2026 (was)
% B2op.R.R2=[];                                                             % MMP, 11/04/2026 (was)

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% % MMP, 11/04/2026 (was)
% % Next TB1op                                                              % MMP, 11/04/2026 (was)
% TB1op.P=zeros(nx,nw);                                                     % MMP, 11/04/2026 (was)
% TB1op.Q1=[];                                                              % MMP, 11/04/2026 (was)
% TB1op.Q2=Tb1;                                                             % MMP, 11/04/2026 (was)
% TB1op.R.R0=[];                                                            % MMP, 11/04/2026 (was)
% TB1op.R.R1=[];                                                            % MMP, 11/04/2026 (was)
% TB1op.R.R2=[];                                                            % MMP, 11/04/2026 (was)

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% % MMP, 11/04/2026 (was)
% % Next TB2op                                                              % MMP, 11/04/2026 (was)
% TB2op.P=zeros(nx,nu);                                                     % MMP, 11/04/2026 (was)
% TB2op.Q1=[];                                                              % MMP, 11/04/2026 (was)
% TB2op.Q2=Tb2;                                                             % MMP, 11/04/2026 (was)
% TB2op.R.R0=[];                                                            % MMP, 11/04/2026 (was)
% TB2op.R.R1=[];                                                            % MMP, 11/04/2026 (was)
% TB2op.R.R2=[];                                                            % MMP, 11/04/2026 (was)

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% % MMP, 11/04/2026 (was)
% % Next C1op                                                               % MMP, 11/04/2026 (was)
% C1op.P=Cb10;                                                              % MMP, 11/04/2026 (was)
% C1op.Q1=Cb11;                                                             % MMP, 11/04/2026 (was)
% C1op.Q2=[];                                                               % MMP, 11/04/2026 (was)
% C1op.R.R0=[];                                                             % MMP, 11/04/2026 (was)
% C1op.R.R1=[];                                                             % MMP, 11/04/2026 (was)
% C1op.R.R2=[];                                                             % MMP, 11/04/2026 (was)

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% % MMP, 11/04/2026 (was)
% % Next C2op                                                               % MMP, 11/04/2026 (was)
% C2op.P=Cb20;                                                              % MMP, 11/04/2026 (was)
% C2op.Q1=Cb21;                                                             % MMP, 11/04/2026 (was)
% C2op.Q2=[];                                                               % MMP, 11/04/2026 (was)
% C2op.R.R0=[];                                                             % MMP, 11/04/2026 (was)
% C2op.R.R1=[];                                                             % MMP, 11/04/2026 (was)
% C2op.R.R2=[];                                                             % MMP, 11/04/2026 (was)

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% % MMP, 11/04/2026 (was)
% % Next D11op                                                              % MMP, 11/04/2026 (was)
% D11op.P=Db11;                                                             % MMP, 11/04/2026 (was)
% D11op.Q1=[];                                                              % MMP, 11/04/2026 (was)
% D11op.Q2=[];                                                              % MMP, 11/04/2026 (was)
% D11op.R.R0=[];                                                            % MMP, 11/04/2026 (was)
% D11op.R.R1=[];                                                            % MMP, 11/04/2026 (was)
% D11op.R.R2=[];                                                            % MMP, 11/04/2026 (was)

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% % MMP, 11/04/2026 (was)
% % Next D12op                                                              % MMP, 11/04/2026 (was)
% D12op.P=Db12;                                                             % MMP, 11/04/2026 (was)
% D12op.Q1=[];                                                              % MMP, 11/04/2026 (was)
% D12op.Q2=[];                                                              % MMP, 11/04/2026 (was)
% D12op.R.R0=[];                                                            % MMP, 11/04/2026 (was)
% D12op.R.R1=[];                                                            % MMP, 11/04/2026 (was)
% D12op.R.R2=[];                                                            % MMP, 11/04/2026 (was)

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% % MMP, 11/04/2026 (was)
% % Next D21op                                                              % MMP, 11/04/2026 (was)
% D21op.P=Db21;                                                             % MMP, 11/04/2026 (was)
% D21op.Q1=[];                                                              % MMP, 11/04/2026 (was)
% D21op.Q2=[];                                                              % MMP, 11/04/2026 (was)
% D21op.R.R0=[];                                                            % MMP, 11/04/2026 (was)
% D21op.R.R1=[];                                                            % MMP, 11/04/2026 (was)
% D21op.R.R2=[];                                                            % MMP, 11/04/2026 (was)

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% % MMP, 11/04/2026 (was)
% % Next D22op                                                              % MMP, 11/04/2026 (was)
% D22op.P=Db22;                                                             % MMP, 11/04/2026 (was)
% D22op.Q1=[];                                                              % MMP, 11/04/2026 (was)
% D22op.Q2=[];                                                              % MMP, 11/04/2026 (was)
% D22op.R.R0=[];                                                            % MMP, 11/04/2026 (was)
% D22op.R.R1=[];                                                            % MMP, 11/04/2026 (was)
% D22op.R.R2=[];                                                            % MMP, 11/04/2026 (was)


% % correction of empty dimensions                                          % MMP, 11/04/2026 (was)
% Top.dim = [nx nx; nxb nxb];                                               % MMP, 11/04/2026 (was)
% Aop.dim = [nx nx; nxb nxb];                                               % MMP, 11/04/2026 (was)
% C1op.dim = [nz nx; 0 nxb];                                                % MMP, 11/04/2026 (was)
% C2op.dim = [ny nx; 0 nxb];                                                % MMP, 11/04/2026 (was)
% B1op.dim=[nx nw; nxb 0];                                                  % MMP, 11/04/2026 (was)
% B2op.dim = [nx nu; nxb 0];                                                % MMP, 11/04/2026 (was)
% D11op.dim = [nz nw; 0 0];                                                 % MMP, 11/04/2026 (was)
% D12op.dim = [nz nu; 0 0];                                                 % MMP, 11/04/2026 (was)
% D21op.dim = [ny nw; 0 0];                                                 % MMP, 11/04/2026 (was)
% D22op.dim = [ny nu; 0 0];                                                 % MMP, 11/04/2026 (was)
% TB1op.dim = [nx nw; nxb 0];                                               % MMP, 11/04/2026 (was)
% TB2op.dim = [nx nu; nxb 0];                                               % MMP, 11/04/2026 (was)
% % % END change MMP, 11/04/2026                                            % MMP, 11/04/2026


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
nx1=nx;nx2=nxb;
a=-1;b=0;
X = [a b];

PIE = pie_struct();
% PIE.space.vars = reshape(pvar2varname(Top.var1),1,[]);                    % MMP, 11/04/2026 (was)
PIE.space.vars = reshape(pvar2varname(s1),1,[]);  % s1 was Top.var1         % MMP, 11/04/2026
PIE.space.dom = X;                                                          % DJ, 10/29/2026
PIE.T = Top;    PIE.Tw = TB1op;     PIE.Tu = TB2op;
PIE.A = Aop;    PIE.B1 = B1op;      PIE.B2 = B2op;
PIE.C1 = C1op;  PIE.D11 = D11op;    PIE.D12 = D12op;
PIE.C2 = C2op;  PIE.D21 = D21op;    PIE.D22 = D22op;

% Make sure dimensions of operators match.
PIE = initialize(PIE);

end
