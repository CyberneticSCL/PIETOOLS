%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% convert_PIETOOLS_DDF.m     PIETOOLS 2022
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function PIE=convert_PIETOOLS_DDF(DDF,out_type)
% This a setup routine for converting DDFs into PIEs
% Input argument "out_type" is optional, and should always be set to ''pie''.
% Converts the DDF formulation to PIE formulation
% The DDF is first completed and checked (undefined terms set to zero).
% With phi_i(t,s) = r_i(t+tau(i)*s), s in [-1,0], the PIE state is
% [x; phi_s], phi_s = d/ds phi, since phi_i(t,s) = r_i(t) - int_s^0
% phi_s(t,eta) d eta and d/dt phi_i = (1/tau(i)) d/ds phi_i.
%
% A Partial Integral Equation is defined by 12 PI operators as
%
% TB1op \dot{w}(t)+TB2op \dot{u}(t)+Top \dot{x}(t)= Aop x(t) + B1op w(t) + B2op u(t)
%                                             z(t)= C1op x(t) + D11op w(t) + D12op u(t)
%                                             y(t)= C2op x(t) + D21op w(t) + D22op u(t)
%
% The formulae were defined in (12) and (15)-(16) in the paper
% ''Representation of Networks and Systems with Delay: DDEs, DDFs, ODE-PDEs and PIEs''
% reference: https://arxiv.org/abs/1910.03881
%
% INPUT
%   DDF:      'struct', the DDF with K = numel(DDF.tau) delays
%     \dot x(t) = A0 x(t) + B1 w(t) + B2 u(t) + Bv v(t)
%          z(t) = C1 x(t) + D11 w(t) + D12 u(t) + D1v v(t)
%          y(t) = C2 x(t) + D21 w(t) + D22 u(t) + D2v v(t)
%        r_i(t) = Cr{i} x(t) + Br1{i} w(t) + Br2{i} u(t) + Drv{i} v(t)
%          v(t) = sum_i Cv{i} r_i(t-tau(i)) + sum_i int_{-tau(i)}^0 Cvd{i}(s) r_i(t+s) ds
%             fields tau (1 x K 'double', delays > 0), the 'double'
%             matrices A0, ..., D2v, and 1 x K cells Cr, Br1, Br2, Drv, Cv
%             of 'double' matrices and Cvd of 'double' or 'polynomial'
%             kernels in s1; undefined fields are zero.
%   out_type: (optional) 'char', 'pie'.
%
% OUTPUT
%   PIE:  'pie_struct', initialized: the twelve 'copvar' operators T, Tw,
%         Tu, A, B1, B2, C1, D11, D12, C2, D21, D22 on R^nx x L2^nr[-1,0],
%         nr = sum_i nr{i}, spatial variable s1 (dummy s1_dum).

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
% For support, contact M. Peet, Arizona State University at mpeet@asu.edu
% Initial coding MP - 10_01_2020
%  MP - 5_30_2021; converted script to function format
% SS - 9/28, added a dimension correction step to replace 0x0 empty
% matrices with empty matrices of correct size
%  DJ - 09/01/2022: Output result as "pie_struct" object.
% DJ, 12/07/2024: Use new default vars (s1,s2) and (s1_dum,s2_dum);
% DJ, MMP, 11/21/2026: Changes since PIETOOLS 2025: the twelve operators of
%   the PIE are 'copvar' containers of 'sopvar' blocks, not 'opvar' objects;
%   the spatial variable and domain are set in PIE.space (vars, dom); help
%   and comments revised.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Convert DDF to 4-PI representation
if nargin==2 && ~strcmpi(out_type,'pie')
    error('Second argument must be ''pie'', or must be omitted. Only conversion from DDF to PIE is supported.')
end
% The following dimensions assume the DDF has been initialized, so lets
% make sure.
DDF=initialize_PIETOOLS_DDF(DDF);
 nK=length(DDF.tau);   % number of delays
 nx=size(DDF.A0,1);   % number of states
 nw=size(DDF.B1,2);   % number of disturbances
 nu=size(DDF.B2,2);   % number of inputs
 nv=size(DDF.Cv{1},1);   % dimension of v
 for i=1:nK
    nr{i}=size(DDF.Cr{i},1);   % dimension of the original delayed channel, i
 end
 nz=size(DDF.C1,1);  % number of regulated outputs
 ny=size(DDF.C2,1);  % number of sensed outputs

 pvar s1 s1_dum                                                             % DJ, 12/07/2024
% First we define the Single-Pipe PIE Representation of the TDS


%Chv{1}=Cv{1}+double(tau(1)*int(subs(Cvd{1},s,tau(1)*s),s,-1,0));

% v in terms of x, w, u and phi_s: with Chv{i} = Cv{i} + tau(i)*
% int_{-1}^0 Cvd{i}(tau(i)*s) ds and DI = (I - sum_i Chv{i}*Drv{i})^-1,
%   v = Cvx x + Dvw w + Dvu u + sum_i int_{-1}^0 CI{i}(s) phi_s,i(s) ds,
%   Cvx = DI*sum_i Chv{i}*Cr{i} (Dvw: Br1{i}, Dvu: Br2{i}),
%   CI{i}(s) = -DI*(Cv{i} + tau(i)*int_{-1}^s Cvd{i}(tau(i)*eta) d eta)
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
% DI, the coefficients Cvx, Dvw, Dvu of v, and the kernels CI{i}
DI=inv(eye(nv)-sumDI);
Cvx=DI*sumCvx;
Dvw=DI*sumDvw;
Dvu=DI*sumDvu;
for i = 1:nK
    CI{i}=-DI*(DDF.Cv{i}+DDF.tau(i)*int(subs(DDF.Cvd{i},s1,s1_dum*DDF.tau(i)),s1_dum,-1,s1));
end

%%%%%%%%%%%%%%%%%%%%%%%

% Stacked channels r = [r_1; ...; r_K] = Tb0 x + Tb1 w + Tb2 u + (the
% phi_s terms of Drv{i}*v): Tb0 rows Cr{i} + Drv{i}*Cvx, and likewise
Tb0=[];
Tb1=[];
Tb2=[];
for i = 1:nK
    Tb0=[Tb0; DDF.Cr{i}+DDF.Drv{i}*Cvx];
    Tb1=[Tb1; DDF.Br1{i}+DDF.Drv{i}*Dvw];
    Tb2=[Tb2; DDF.Br2{i}+DDF.Drv{i}*Dvu];
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% d/dt phi = Itau*phi_s, Itau = blkdiag(I/tau(i))
Itau=[];
for i=1:nK
    Itau=blkdiag(Itau,eye(nr{i})/DDF.tau(i));
end
nxb=size(Itau,2); % number of infinite-dimensional states

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%
% Kernels of phi = T*[x; phi_s] + Tw*w + Tu*u: Tbai = [CI{1} ... CI{K}],
% Tba the rows Drv{i}*Tbai in s1_dum (on [-1,s1]), Tbb = Tba - I (on
% [s1,0])
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

% Equations of x, z, y with v substituted: the multipliers (A0 + Bv*Cvx,
% ...) and the integral kernels Bv*Tbai, D1v*Tbai, D2v*Tbai of phi_s
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

% DJ, MMP, 11/21/2026: BEGIN the operators as 'copvar' containers
% Each operator [P, Q1; Q2, R]: R^p0 x L2^p1[s1] -> R^q0 x L2^q1[s1], R the % DJ, MMP, 11/21/2026
% multiplier R0 plus the integrals of R1 on [-1,s1] and R2 on [s1,0], as a  % DJ, MMP, 11/21/2026
% 'copvar' container; arguments P, Q1, Q2, {R0,R1,R2}, [q0 p0; q1 p1]       % DJ, MMP, 11/21/2026
I1 = [-1 0];                       % domain of s1                           % DJ, MMP, 11/21/2026
Top   = pi4_copvar_sop(eye(nx),zeros(nx,nxb),Tb0,{zeros(nxb),Tba,Tbb},...
                       [nx nx; nxb nxb],s1,s1_dum,I1);                      % DJ, MMP, 11/21/2026
Aop   = pi4_copvar_sop(Ab0,Ab,zeros(nxb,nx),{Itau,zeros(nxb),zeros(nxb)},...
                       [nx nx; nxb nxb],s1,s1_dum,I1);                      % DJ, MMP, 11/21/2026
B1op  = pi4_copvar_sop(Bb1,[],zeros(nxb,nw),{[],[],[]},...
                       [nx nw; nxb 0],s1,s1_dum,I1);                        % DJ, MMP, 11/21/2026
B2op  = pi4_copvar_sop(Bb2,[],zeros(nxb,nu),{[],[],[]},...
                       [nx nu; nxb 0],s1,s1_dum,I1);                        % DJ, MMP, 11/21/2026
TB1op = pi4_copvar_sop(zeros(nx,nw),[],Tb1,{[],[],[]},...
                       [nx nw; nxb 0],s1,s1_dum,I1);                        % DJ, MMP, 11/21/2026
TB2op = pi4_copvar_sop(zeros(nx,nu),[],Tb2,{[],[],[]},...
                       [nx nu; nxb 0],s1,s1_dum,I1);                        % DJ, MMP, 11/21/2026
C1op  = pi4_copvar_sop(Cb10,Cb11,[],{[],[],[]},...
                       [nz nx; 0 nxb],s1,s1_dum,I1);                        % DJ, MMP, 11/21/2026
C2op  = pi4_copvar_sop(Cb20,Cb21,[],{[],[],[]},...
                       [ny nx; 0 nxb],s1,s1_dum,I1);                        % DJ, MMP, 11/21/2026
D11op = pi4_copvar_sop(Db11,[],[],{[],[],[]},[nz nw; 0 0],s1,s1_dum,I1);    % DJ, MMP, 11/21/2026
D12op = pi4_copvar_sop(Db12,[],[],{[],[],[]},[nz nu; 0 0],s1,s1_dum,I1);    % DJ, MMP, 11/21/2026
D21op = pi4_copvar_sop(Db21,[],[],{[],[],[]},[ny nw; 0 0],s1,s1_dum,I1);    % DJ, MMP, 11/21/2026
D22op = pi4_copvar_sop(Db22,[],[],{[],[],[]},[ny nu; 0 0],s1,s1_dum,I1);    % DJ, MMP, 11/21/2026
% Former code:                                                              % DJ, MMP, 11/21/2026
% opvar Top Aop C1op C2op B1op B2op D11op D12op D21op D22op TB1op TB2op;    % DJ, MMP, 11/21/2026 (was)
% Top.I = [-1 0]; Top.dim = [nx nx; nxb nxb]; Top.var1 = s1; Top.var2 = s1_dum; % DJ, MMP, 11/21/2026 (was)
% Aop.I = [-1 0]; Aop.dim = [nx nx; nxb nxb]; Aop.var1 = s1; Aop.var2 = s1_dum; % DJ, MMP, 11/21/2026 (was)
% C1op.I = [-1 0]; C1op.dim = [nz nx; 0 nxb]; C1op.var1 = s1; C1op.var2 = s1_dum; % DJ, MMP, 11/21/2026 (was)
% C2op.I = [-1 0]; C2op.dim = [ny nx; 0 nxb]; C2op.var1 = s1; C2op.var2 = s1_dum; % DJ, MMP, 11/21/2026 (was)
% B1op.I = [-1 0]; B1op.dim=[nx nw; nxb 0]; B1op.var1 = s1; B1op.var2 = s1_dum; % DJ, MMP, 11/21/2026 (was)
% B2op.I = [-1 0]; B2op.dim = [nx nu; nxb 0]; B2op.var1 = s1; B2op.var2 = s1_dum; % DJ, MMP, 11/21/2026 (was)
% D11op.I = [-1 0]; D11op.dim = [nz nw; 0 0]; D11op.var1 = s1; D11op.var2 = s1_dum; % DJ, MMP, 11/21/2026 (was)
% D12op.I = [-1 0]; D12op.dim = [nz nu; 0 0]; D12op.var1 = s1; D12op.var2 = s1_dum; % DJ, MMP, 11/21/2026 (was)
% D21op.I = [-1 0]; D21op.dim = [ny nw; 0 0]; D21op.var1 = s1; D21op.var2 = s1_dum; % DJ, MMP, 11/21/2026 (was)
% D22op.I = [-1 0]; D22op.dim = [ny nu; 0 0]; D22op.var1 = s1; D22op.var2 = s1_dum; % DJ, MMP, 11/21/2026 (was)
% TB1op.I = [-1 0]; TB1op.dim = [nx nw; nxb 0]; TB1op.var1 = s1; TB1op.var2 = s1_dum; % DJ, MMP, 11/21/2026 (was)
% TB2op.I = [-1 0]; TB2op.dim = [nx nu; nxb 0]; TB2op.var1 = s1; TB2op.var2 = s1_dum; % DJ, MMP, 11/21/2026 (was)


% Aop.P=Ab0;                                                                % DJ, MMP, 11/21/2026 (was)
% Aop.Q1=Ab;                                                                % DJ, MMP, 11/21/2026 (was)
% Aop.Q2=zeros(nxb,nx);                                                     % DJ, MMP, 11/21/2026 (was)
% Aop.R.R0=Itau;                                                            % DJ, MMP, 11/21/2026 (was)
% Aop.R.R1=zeros(nxb);                                                      % DJ, MMP, 11/21/2026 (was)
% Aop.R.R2=zeros(nxb);                                                      % DJ, MMP, 11/21/2026 (was)
% Top.P=eye(nx);                                                            % DJ, MMP, 11/21/2026 (was)
% Top.Q1=zeros(nx,nxb);                                                     % DJ, MMP, 11/21/2026 (was)
% Top.Q2=Tb0;                                                               % DJ, MMP, 11/21/2026 (was)
% Top.R.R0=zeros(nxb);                                                      % DJ, MMP, 11/21/2026 (was)
% Top.R.R1=Tba;                                                             % DJ, MMP, 11/21/2026 (was)
% Top.R.R2=Tbb;                                                             % DJ, MMP, 11/21/2026 (was)

% B1op.P=Bb1;                                                               % DJ, MMP, 11/21/2026 (was)
% B1op.Q1=[];                                                               % DJ, MMP, 11/21/2026 (was)
% B1op.Q2=zeros(nxb,nw);                                                    % DJ, MMP, 11/21/2026 (was)
% B1op.R.R0=[];                                                             % DJ, MMP, 11/21/2026 (was)
% B1op.R.R1=[];                                                             % DJ, MMP, 11/21/2026 (was)
% B1op.R.R2=[];                                                             % DJ, MMP, 11/21/2026 (was)

% B2op.P=Bb2;                                                               % DJ, MMP, 11/21/2026 (was)
% B2op.Q1=[];                                                               % DJ, MMP, 11/21/2026 (was)
% B2op.Q2=zeros(nxb,nu);                                                    % DJ, MMP, 11/21/2026 (was)
% B2op.R.R0=[];                                                             % DJ, MMP, 11/21/2026 (was)
% B2op.R.R1=[];                                                             % DJ, MMP, 11/21/2026 (was)
% B2op.R.R2=[];                                                             % DJ, MMP, 11/21/2026 (was)

% TB1op.P=zeros(nx,nw);                                                     % DJ, MMP, 11/21/2026 (was)
% TB1op.Q1=[];                                                              % DJ, MMP, 11/21/2026 (was)
% TB1op.Q2=Tb1;                                                             % DJ, MMP, 11/21/2026 (was)
% TB1op.R.R0=[];                                                            % DJ, MMP, 11/21/2026 (was)
% TB1op.R.R1=[];                                                            % DJ, MMP, 11/21/2026 (was)
% TB1op.R.R2=[];                                                            % DJ, MMP, 11/21/2026 (was)

% TB2op.P=zeros(nx,nu);                                                     % DJ, MMP, 11/21/2026 (was)
% TB2op.Q1=[];                                                              % DJ, MMP, 11/21/2026 (was)
% TB2op.Q2=Tb2;                                                             % DJ, MMP, 11/21/2026 (was)
% TB2op.R.R0=[];                                                            % DJ, MMP, 11/21/2026 (was)
% TB2op.R.R1=[];                                                            % DJ, MMP, 11/21/2026 (was)
% TB2op.R.R2=[];                                                            % DJ, MMP, 11/21/2026 (was)

% C1op.P=Cb10;                                                              % DJ, MMP, 11/21/2026 (was)
% C1op.Q1=Cb11;                                                             % DJ, MMP, 11/21/2026 (was)
% C1op.Q2=[];                                                               % DJ, MMP, 11/21/2026 (was)
% C1op.R.R0=[];                                                             % DJ, MMP, 11/21/2026 (was)
% C1op.R.R1=[];                                                             % DJ, MMP, 11/21/2026 (was)
% C1op.R.R2=[];                                                             % DJ, MMP, 11/21/2026 (was)

% C2op.P=Cb20;                                                              % DJ, MMP, 11/21/2026 (was)
% C2op.Q1=Cb21;                                                             % DJ, MMP, 11/21/2026 (was)
% C2op.Q2=[];                                                               % DJ, MMP, 11/21/2026 (was)
% C2op.R.R0=[];                                                             % DJ, MMP, 11/21/2026 (was)
% C2op.R.R1=[];                                                             % DJ, MMP, 11/21/2026 (was)
% C2op.R.R2=[];                                                             % DJ, MMP, 11/21/2026 (was)

% D11op.P=Db11;                                                             % DJ, MMP, 11/21/2026 (was)
% D11op.Q1=[];                                                              % DJ, MMP, 11/21/2026 (was)
% D11op.Q2=[];                                                              % DJ, MMP, 11/21/2026 (was)
% D11op.R.R0=[];                                                            % DJ, MMP, 11/21/2026 (was)
% D11op.R.R1=[];                                                            % DJ, MMP, 11/21/2026 (was)
% D11op.R.R2=[];                                                            % DJ, MMP, 11/21/2026 (was)

% D12op.P=Db12;                                                             % DJ, MMP, 11/21/2026 (was)
% D12op.Q1=[];                                                              % DJ, MMP, 11/21/2026 (was)
% D12op.Q2=[];                                                              % DJ, MMP, 11/21/2026 (was)
% D12op.R.R0=[];                                                            % DJ, MMP, 11/21/2026 (was)
% D12op.R.R1=[];                                                            % DJ, MMP, 11/21/2026 (was)
% D12op.R.R2=[];                                                            % DJ, MMP, 11/21/2026 (was)

% D21op.P=Db21;                                                             % DJ, MMP, 11/21/2026 (was)
% D21op.Q1=[];                                                              % DJ, MMP, 11/21/2026 (was)
% D21op.Q2=[];                                                              % DJ, MMP, 11/21/2026 (was)
% D21op.R.R0=[];                                                            % DJ, MMP, 11/21/2026 (was)
% D21op.R.R1=[];                                                            % DJ, MMP, 11/21/2026 (was)
% D21op.R.R2=[];                                                            % DJ, MMP, 11/21/2026 (was)

% D22op.P=Db22;                                                             % DJ, MMP, 11/21/2026 (was)
% D22op.Q1=[];                                                              % DJ, MMP, 11/21/2026 (was)
% D22op.Q2=[];                                                              % DJ, MMP, 11/21/2026 (was)
% D22op.R.R0=[];                                                            % DJ, MMP, 11/21/2026 (was)
% D22op.R.R1=[];                                                            % DJ, MMP, 11/21/2026 (was)
% D22op.R.R2=[];                                                            % DJ, MMP, 11/21/2026 (was)


% Top.dim = [nx nx; nxb nxb];                                               % DJ, MMP, 11/21/2026 (was)
% Aop.dim = [nx nx; nxb nxb];                                               % DJ, MMP, 11/21/2026 (was)
% C1op.dim = [nz nx; 0 nxb];                                                % DJ, MMP, 11/21/2026 (was)
% C2op.dim = [ny nx; 0 nxb];                                                % DJ, MMP, 11/21/2026 (was)
% B1op.dim=[nx nw; nxb 0];                                                  % DJ, MMP, 11/21/2026 (was)
% B2op.dim = [nx nu; nxb 0];                                                % DJ, MMP, 11/21/2026 (was)
% D11op.dim = [nz nw; 0 0];                                                 % DJ, MMP, 11/21/2026 (was)
% D12op.dim = [nz nu; 0 0];                                                 % DJ, MMP, 11/21/2026 (was)
% D21op.dim = [ny nw; 0 0];                                                 % DJ, MMP, 11/21/2026 (was)
% D22op.dim = [ny nu; 0 0];                                                 % DJ, MMP, 11/21/2026 (was)
% TB1op.dim = [nx nw; nxb 0];                                               % DJ, MMP, 11/21/2026 (was)
% TB2op.dim = [nx nu; nxb 0];                                               % DJ, MMP, 11/21/2026 (was)
% DJ, MMP, 11/21/2026: END the operators as 'copvar' containers


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% The PIE: spatial variable s1 on [-1,0] and the twelve operators
nx1=nx;nx2=nxb;
a=-1;b=0;
X = [a b];

PIE = pie_struct();
PIE.space.vars = reshape(pvar2varname(s1),1,[]);                            % DJ, MMP, 11/21/2026
PIE.space.dom = X;                                                          % DJ, MMP, 11/21/2026
PIE.T = Top;    PIE.Tw = TB1op;     PIE.Tu = TB2op;
PIE.A = Aop;    PIE.B1 = B1op;      PIE.B2 = B2op;
PIE.C1 = C1op;  PIE.D11 = D11op;    PIE.D12 = D12op;
PIE.C2 = C2op;  PIE.D21 = D21op;    PIE.D22 = D22op;

% Make sure dimensions of operators match.
PIE = initialize(PIE);

end
