%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% convert_PIETOOLS_DDE2PIE.m     PIETOOLS 2022
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function PIE=convert_PIETOOLS_DDE2PIE_legacy(DDE)
% NAIVE Conversion of the DDE formulation to PIE formulation
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
% Initial coding MP - 7_01_2020
%  MP - 6_1_2021:  made a function. Changed name. Added warning
%  SS - 9/28: added correction to fix 0x0 empty matrices by replacing with empty matrices of correct size 
%  DJ - 09/01/2022: Output result as "pie_struct" object.
% DJ, 12/07/2024: Use new default vars s1 and s1_dum;
% MMP, 11/04/2026: Build the twelve operators as 'copvar' containers of
%                   'sopvar' blocks (PI4_COPVAR_SOP) instead of 'opvar'
%                   objects; no 'opvar' is constructed, as the class is
%                   being removed from PIETOOLS. The input dimension of D12,
%                   D22 and TB2 is nu, not nw: with nu ~= nw the former dims
%                   gave operators whose u spaces disagree, so INITIALIZE
%                   rejected the PIE (12 of the 26 library DDEs). The PIE
%                   is now initialized here, as by convert_PIETOOLS_DDF, and
%                   its variables are set in PIE.space (DJ, 10/29/2026).

disp('WARNING: convert_PIETOOLS_DDE2PIE is no longer supported! Direct conversion of DDE to PIE is NOT recommended and may result in very large computation times!!! Please use minimize_PIETOOLS_DDE2PIE instead!')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% The following dimensions should have been defined in the
% initialize_PIETOOLS_DDE script
DDE=initialize_PIETOOLS_DDE(DDE);
nx=size(DDE.A0,1);   % number of states
nw=size(DDE.B1,2);   % number of disturbances
nu=size(DDE.B2,2);   % number of inputs
nz=size(DDE.C1,1);  % number of regulated outputs
ny=size(DDE.C2,1);  % number of sensed outputs
nK=length(DDE.tau);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

pvar s1 s1_dum                                                              % DJ, 12/07/2024
% First we define the Single-Pipe PIE Representation of the DDE
Itau=[];
Tb0=[];
Tb1=[];
Tb2=[];
for i=1:nK
    Itau=blkdiag(Itau,eye(nx+nw+nu)/DDE.tau(i)); % double check this
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
Ab0=DDE.A0;
Ab=[];
Bb1=DDE.B1;
Bb2=DDE.B2;
Cb10=DDE.C1;
Cb20=DDE.C2;
Cb11=[];
Cb21=[];
Db11=DDE.D11;
Db12=DDE.D12;
Db21=DDE.D21;
Db22=DDE.D22;
for i = 1:nK
  XAi=[DDE.Ai{i} DDE.B1i{i} DDE.B2i{i}]+DDE.tau(i)*int(subs([DDE.Adi{i} DDE.B1di{i} DDE.B2di{i}],s1,DDE.tau(i)*s1),s1,-1,s1);
  XC1i=[DDE.C1i{i} DDE.D11i{i} DDE.D12i{i}]+DDE.tau(i)*int(subs([DDE.C1di{i} DDE.D11di{i} DDE.D12di{i}],s1,DDE.tau(i)*s1),s1,-1,s1);
  XC2i=[DDE.C2i{i} DDE.D21i{i} DDE.D22i{i}]+DDE.tau(i)*int(subs([DDE.C2di{i} DDE.D21di{i} DDE.D22di{i}],s1,DDE.tau(i)*s1),s1,-1,s1);
  Ab0 = Ab0 + DDE.Ai{i} + int(DDE.tau(i)*subs(DDE.Adi{i},s1,DDE.tau(i)*s1),s1,-1,0);
  Ab=[Ab -XAi];
  Bb1 = Bb1 + DDE.B1i{i} + int(DDE.tau(i)*subs(DDE.B1di{i},s1,DDE.tau(i)*s1),s1,-1,0);
  Bb2 = Bb2 + DDE.B2i{i} + int(DDE.tau(i)*subs(DDE.B2di{i},s1,DDE.tau(i)*s1),s1,-1,0);
  Db11 = Db11 + DDE.D11i{i} + int(DDE.tau(i)*subs(DDE.D11di{i},s1,DDE.tau(i)*s1),s1,-1,0);
  Db12 = Db12 + DDE.D12i{i} + int(DDE.tau(i)*subs(DDE.D12di{i},s1,DDE.tau(i)*s1),s1,-1,0);
  Db21 = Db21 + DDE.D21i{i} + int(DDE.tau(i)*subs(DDE.D21di{i},s1,DDE.tau(i)*s1),s1,-1,0);
  Db22 = Db22 + DDE.D22i{i} + int(DDE.tau(i)*subs(DDE.D22di{i},s1,DDE.tau(i)*s1),s1,-1,0);
  Cb11=[Cb11 -XC1i];
  Cb21=[Cb21 -XC2i];
  Cb10 = Cb10 + DDE.C1i{i} + int(DDE.tau(i)*subs(DDE.C1di{i},s1,DDE.tau(i)*s1),s1,-1,0);
  Cb20 = Cb20 + DDE.C2i{i} + int(DDE.tau(i)*subs(DDE.C2di{i},s1,DDE.tau(i)*s1),s1,-1,0);
  clear XAi XC1i XC2i
end
nxb=(nx+nw+nu)*nK; % number of infinite-dimensional states

% % % BEGIN change MMP, 11/04/2026: the operators as containers.            % MMP, 11/04/2026
% Each operator [P, Q1; Q2, R]: R^p0 x L2^p1[s1] -> R^q0 x L2^q1[s1], R the % MMP, 11/04/2026
% multiplier R0 plus the integrals of R1 on [-1,s1] and R2 on [s1,0], is    % MMP, 11/04/2026
% built by PI4_COPVAR_SOP from the components of the former opvar, block    % MMP, 11/04/2026
% for block as OP2COPVAR_SOP converted that opvar in INITIALIZE.            % MMP, 11/04/2026
% Arguments: P, Q1, Q2, {R0,R1,R2}, dim = [q0 p0; q1 p1] as the opvar dim.  % MMP, 11/04/2026
% The u channel has nu columns in D12, D22 and TB2; the former dims gave    % MMP, 11/04/2026
% nw (D12, D22 when declared, D22 and TB2 when corrected).                  % MMP, 11/04/2026
% Former code, kept as comments after the new: the opvar declarations and   % MMP, 11/04/2026
% fields (MP 7_01_2020; variables DJ 12/07/2024), dims (SS 9/28).           % MMP, 11/04/2026
I1 = [-1 0];                       % domain of s1, the former opvar field I % MMP, 11/04/2026
Top   = pi4_copvar_sop(eye(nx),zeros(nx,nxb),Tb0,{zeros(nxb),Tba,Tbb},...
                       [nx nx; nxb nxb],s1,s1_dum,I1);                      % MMP, 11/04/2026
TB1op = pi4_copvar_sop(zeros(nx,nw),[],Tb1,{[],[],[]},...
                       [nx nw; nxb 0],s1,s1_dum,I1);                        % MMP, 11/04/2026
TB2op = pi4_copvar_sop(zeros(nx,nu),[],Tb2,{[],[],[]},...
                       [nx nu; nxb 0],s1,s1_dum,I1);                        % MMP, 11/04/2026
Aop   = pi4_copvar_sop(Ab0,Ab,zeros(nxb,nx),{Itau,zeros(nxb),zeros(nxb)},...
                       [nx nx; nxb nxb],s1,s1_dum,I1);                      % MMP, 11/04/2026
C1op  = pi4_copvar_sop(Cb10,Cb11,[],{[],[],[]},...
                       [nz nx; 0 nxb],s1,s1_dum,I1);                        % MMP, 11/04/2026
C2op  = pi4_copvar_sop(Cb20,Cb21,[],{[],[],[]},...
                       [ny nx; 0 nxb],s1,s1_dum,I1);                        % MMP, 11/04/2026
B1op  = pi4_copvar_sop(Bb1,[],zeros(nxb,nw),{[],[],[]},...
                       [nx nw; nxb 0],s1,s1_dum,I1);                        % MMP, 11/04/2026
B2op  = pi4_copvar_sop(Bb2,[],zeros(nxb,nu),{[],[],[]},...
                       [nx nu; nxb 0],s1,s1_dum,I1);                        % MMP, 11/04/2026
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
% D12op.I = [-1 0]; D12op.dim = [nz nw; 0 0]; D12op.var1 = s1; D12op.var2 = s1_dum; % MMP, 11/04/2026 (was)
% D21op.I = [-1 0]; D21op.dim = [ny nw; 0 0]; D21op.var1 = s1; D21op.var2 = s1_dum; % MMP, 11/04/2026 (was)
% D22op.I = [-1 0]; D22op.dim = [ny nw; 0 0]; D22op.var1 = s1; D22op.var2 = s1_dum; % MMP, 11/04/2026 (was)
% TB1op.I = [-1 0]; TB1op.dim = [nx nw; nxb 0]; TB1op.var1 = s1; TB1op.var2 = s1_dum; % MMP, 11/04/2026 (was)
% TB2op.I = [-1 0]; TB2op.dim = [nx nu; nxb 0]; TB2op.var1 = s1; TB2op.var2 = s1_dum; % MMP, 11/04/2026 (was)

% % We now define the 4-PI system operators                                 % MMP, 11/04/2026 (was)
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% % MMP, 11/04/2026 (was)
% % Start with Top                                                          % MMP, 11/04/2026 (was)
% Top.P=eye(nx);                                                            % MMP, 11/04/2026 (was)
% Top.Q1=zeros(nx,nxb);                                                     % MMP, 11/04/2026 (was)
% Top.Q2=Tb0;%oneK;                                                         % MMP, 11/04/2026 (was)
% Top.R.R0=zeros(nxb);                                                      % MMP, 11/04/2026 (was)
% Top.R.R1=Tba;%zeros(nxb);                                                 % MMP, 11/04/2026 (was)
% Top.R.R2=Tbb;%-eye(nxb);                                                  % MMP, 11/04/2026 (was)
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
% % Next Aop                                                                % MMP, 11/04/2026 (was)
% Aop.P=Ab0;                                                                % MMP, 11/04/2026 (was)
% Aop.Q1=Ab;                                                                % MMP, 11/04/2026 (was)
% Aop.Q2=zeros(nxb,nx);                                                     % MMP, 11/04/2026 (was)
% Aop.R.R0=Itau;                                                            % MMP, 11/04/2026 (was)
% Aop.R.R1=zeros(nxb);                                                      % MMP, 11/04/2026 (was)
% Aop.R.R2=zeros(nxb);                                                      % MMP, 11/04/2026 (was)
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
% D22op.dim = [ny nw; 0 0];                                                 % MMP, 11/04/2026 (was)
% TB1op.dim = [nx nw; nxb 0];                                               % MMP, 11/04/2026 (was)
% TB2op.dim = [nx nw; nxb 0];                                               % MMP, 11/04/2026 (was)
% % % END change MMP, 11/04/2026                                            % MMP, 11/04/2026


PIE = pie_struct();
% PIE.vars = [Top.var1,Top.var2];                                           % MMP, 11/04/2026 (was)
% PIE.dom=[-1 0];                                                           % MMP, 11/04/2026 (was)
PIE.space.vars = reshape(pvar2varname(s1),1,[]);                            % MMP, 11/04/2026
PIE.space.dom = [-1 0];                                                     % MMP, 11/04/2026
PIE.T = Top; PIE.Tw = TB1op; PIE.Tu = TB2op;
PIE.A = Aop; PIE.B1 = B1op; PIE.B2 = B2op;
PIE.C1 = C1op; PIE.D11 = D11op; PIE.D12 = D12op;
PIE.C2 = C2op; PIE.D21 = D21op; PIE.D22 = D22op;
% Check the spaces of the operators, as convert_PIETOOLS_DDF does           % MMP, 11/04/2026
PIE = initialize(PIE);                                                      % MMP, 11/04/2026

end