%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% convert_PIETOOLS_NDS2DDF.m     PIETOOLS 2022
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function DDF=convert_PIETOOLS_NDS2DDF(NDS)
% This script takes a fully formed NDS representation and find a
% DDF representation of that NDS. 
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
%
% The use of this script is recommended as it exploits the
% structure of the system, as opposed to convert_PIETOOLS_DDE. 
% Specifically, in that script, all states and inputs are delayed,
% forming a rather large-dimensional delay channel
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Convert NDS to a nominal DDF representation
%
%For each delay in tau, the converter will create a delayed channel
%consisting of 
%
% ri(t)=Rs[x(t);
%          w(t);
%          u(t);
%          \dot x(t)];

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


% re-run initialize to make sure everything is defined.
NDS=initialize_PIETOOLS_NDS(NDS);
nK=length(NDS.tau);   % number of delays
nx=size(NDS.A0,1);   % number of states
nw=size(NDS.B1,2);   % number of disturbances
nu=size(NDS.B2,2);   % number of inputs
nz=size(NDS.C1,1);  % number of regulated outputs
ny=size(NDS.C2,1);  % number of sensed outputs

% Construct the DDF representation
% 
% First: A0, B1, B2, C1, D11, D12, C2, D21, D21 are unchanged from the NDS
DDF.tau=NDS.tau;DDF.A0=NDS.A0;DDF.B1=NDS.B1;DDF.B2=NDS.B2;DDF.C1=NDS.C1;DDF.C2=NDS.C2;DDF.D11=NDS.D11;DDF.D12=NDS.D12;DDF.D21=NDS.D21;DDF.D22=NDS.D22;

% However, we must now define:
% Bv, D1v, D2v, Cr{i}, Br1{i}, Br2{i}, Drv{i}, Cv{i}, Cvd{i} 

% The first step is to form the vertical concatenation matrices for each
% delay, tau(i)
%%% ADD case logic for nK=0
if nK==0
    DDF=convert_PIETOOLS_DDE(NDS);
else
    
    nv=nx+nz+ny; nt=nx+nw+nu;
    for i=1:nK
        nr{i}=nt+nx;
        DDF.Cv{i}=[NDS.Ai{i} NDS.B1i{i} NDS.B2i{i} NDS.Ei{i};
               NDS.C1i{i}  NDS.D11i{i} NDS.D12i{i} NDS.C1ei{i};
               NDS.C2i{i}  NDS.D21i{i} NDS.D22i{i} NDS.C2ei{i}];
        Cvd{i}=[NDS.Adi{i} NDS.B1di{i} NDS.B2di{i} NDS.Edi{i};
            NDS.C1di{i}  NDS.D11di{i} NDS.D12di{i} NDS.C1dei{i};
            NDS.C2di{i}  NDS.D21di{i} NDS.D22di{i} NDS.C2dei{i}];
        DDF.Cr{i}=[eye(nx);zeros(nw,nx);zeros(nu,nx);NDS.A0];
        DDF.Br1{i}=[zeros(nx,nw);eye(nw);zeros(nu,nw);NDS.B1];
        DDF.Br2{i}=[zeros(nx,nu);zeros(nw,nu);eye(nu);NDS.B2];
        DDF.Drv{i}=[zeros(nt,nx) zeros(nt,(nz+ny));eye(nx) zeros(nx,(nz+ny))];
    end
    DDF.Bv=[eye(nx) zeros(nx,nz) zeros(nx,ny)];
    DDF.D1v=[zeros(nz,nx) eye(nz) zeros(nz,ny)];
    DDF.D2v=[zeros(ny,nx) zeros(ny,nz) eye(ny)];
% The DDF formulation is now complete. We can call
% convert_PIETOOLS_DDF (initialize is not necessary)
% We should probably also call minimize_PIETOOLS_DDF to reduce the state
% dimension

DDF_out=DDF;
end
