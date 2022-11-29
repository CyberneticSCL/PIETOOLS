%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% minimize_PIETOOLS_DDF.m     PIETOOLS 2022
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function DDF_out=minimize_PIETOOLS_DDF(DDF)
%
% This script takes a fully formed DDF representation and find the minimal
% DDF realization of that DDF.
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
% Convert DDF to a minimal DDF representation
%
%For each delay in tau, the converter will create a delayed channel
%consisting of
%
% ri(t)=Rs[x(t);
%          w(t);
%          u(t)];
%
% The converter will retain the original variables in x,w,u,z,y,v. However,
% the dimensions nri will change, as will the physical meaning of the channels r_i

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


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Construct the DDF representation
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% The following dimensions should have been defined in the
% initialize_PIETOOLS_DDE. But, lets make sure the DDF structure is
% initialized.
DDF=initialize_PIETOOLS_DDF(DDF);
 nK=length(DDF.tau);   % number of delays
 nx=size(DDF.A0,1);   % number of states
 nw=size(DDF.B1,2);   % number of disturbances
 nu=size(DDF.B2,2);   % number of inputs
 nv=size(DDF.Cv{1},1);   % number of inputs
 for i=1:nK
    nr{i}=size(DDF.Cr{i},1);   % dimension of the original delayed channel, i
 end
 nz=size(DDF.C1,1);  % number of regulated outputs
 ny=size(DDF.C2,1);  % number of sensed outputs
 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% First: A0, B1, B2, C1, C2, D11, D12, D21, D21, Bv, D1v, D2v, are
% unchanged from the original DDF. 

DDF_out=DDF;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% However, we will re-define:
% Cr{i}, Br1{i}, Br2{i}, Drv{i}, Cv{i}, Cvd{i}

% The first step is to form the vertical concatenation matrices for each
% delay, tau(i)
%%% ADD case logic for nK=0
if nK==0
    disp('no delay channels detected')
else
    %    nv=nx+nz+ny;
    for i=1:nK
        T=DDF.Cv{i}*[DDF.Cr{i} DDF.Br1{i} DDF.Br2{i} DDF.Drv{i}];
        [ZthT,HT] = poly_decomp_1D_econ(DDF.Cvd{i}'); % so that Cvd{i}'=HT*ZthT or Cvd{i}=ZthT'*HT'
        Tdh=HT'*[DDF.Cr{i} DDF.Br1{i} DDF.Br2{i} DDF.Drv{i}];
        Zth=ZthT';
        [U,D,V]=svd([T;Tdh],'econ'); % so that [T;Tdh]=U*D*V'
        rT=rank(D);  % number of non-zero singular values
        Us=U(:,1:rT);  % First rT columns of U
        Vs=V(:,1:rT);  % First rT columns of V
        Dd=diag(D);
        Ds=diag(Dd(1:rT)); % Diagonal Matrix with rT singular values
        LsF=Us;   % Left side of UV docomposition
        RsF=Ds*Vs';     % Right side of UV decomposition
        Ls=LsF(1:nv,:); % breaking the LHS into upper and lower part - upper part
        Lsd=LsF((nv+1):end,:);   % Lower part
        DDF_out.Cv{i}=Ls;     % permanently re-assign Cv{i}
        DDF_out.Cvd{i}=Zth*Lsd;  % permanently re-assign Cvd{i}
        % At this point, we have the channel r_i=Rs*[x;w;u;v]
        % with output Ls{i}r_i(t-\tau_i)+int_{-\tau_i}^0 Zth{i}*Lsd{i}r_i(t+s)ds
        %nr{i}=size(RsF,1);  % permanently re-assign nr{i}
        DDF_out.Cr{i}=RsF(:,1:nx);  % permanently re-assign Cr{i}
        DDF_out.Br1{i}=RsF(:,(nx+1):(nx+nw));   % permanently re-assign Br1{i}
        DDF_out.Br2{i}=RsF(:,(nx+nw+1):(nx+nw+nu));  % permanently re-assign Br2{i}
        DDF_out.Drv{i}=RsF(:,(nx+nw+nu+1):end);  % permanently re-assign Drv{i}
        %       clear Pb Pbd UT DT Vt rT UsT VsT DdT DsT LsT
    end
    % The DDF re-formulation is now complete. We can call
    % convert_PIETOOLS_DDF (initialize is not necessary)
    %    convert_PIETOOLS_DDF; % This constructs the PIE representation.
end
end
