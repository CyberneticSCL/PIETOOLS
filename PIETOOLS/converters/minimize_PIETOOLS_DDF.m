%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% minimize_PIETOOLS_DDF.m     PIETOOLS 2022
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function DDF_out=minimize_PIETOOLS_DDF(DDF)
%
% This function takes a fully formed DDF representation and finds the minimal
% DDF realization of that DDF.
%
% For each delay tau(i), let Cvd{i}(s) = Zth(s)*H with Zth(s) a monomial
% basis, M = [Cr{i} Br1{i} Br2{i} Drv{i}], and [Cv{i}*M; H*M] = [Ls; Lsd]*Rs
% a rank factorization (SVD). The delay channel is replaced by
%   r_i(t) = Rs*[x(t); w(t); u(t); v(t)],    nr{i} = rank([Cv{i}*M; H*M]),
% with Cv{i} = Ls and Cvd{i}(s) = Zth(s)*Lsd: the output of channel i,
% hence v, is unchanged.
%
% The formulae were defined in (12) and (15)-(16) in the paper
% ''Representation of Networks and Systems with Delay: DDEs, DDFs, ODE-PDEs and PIEs''
% reference: https://arxiv.org/abs/1910.03881
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Convert DDF to a minimal DDF representation
%
%For each delay in tau, the converter will create a delayed channel
%consisting of
%
% ri(t)=Rs[x(t);
%          w(t);
%          u(t);
%          v(t)];
%
% The converter will retain the original variables in x,w,u,z,y,v. However,
% the dimensions nri will change, as will the physical meaning of the channels r_i
%
% INPUT
%   DDF:      'struct', the DDF with K = numel(DDF.tau) delays
%     \dot x(t) = A0 x(t) + B1 w(t) + B2 u(t) + Bv v(t)
%          z(t) = C1 x(t) + D11 w(t) + D12 u(t) + D1v v(t)
%          y(t) = C2 x(t) + D21 w(t) + D22 u(t) + D2v v(t)
%        r_i(t) = Cr{i} x(t) + Br1{i} w(t) + Br2{i} u(t) + Drv{i} v(t)
%          v(t) = sum_i Cv{i} r_i(t-tau(i)) + sum_i int_{-tau(i)}^0 Cvd{i}(s) r_i(t+s) ds
%             fields tau (1 x K 'double'), the 'double' matrices A0, ...,
%             D2v, and 1 x K cells Cr, Br1, Br2, Drv, Cv of 'double'
%             matrices and Cvd of 'double' or 'polynomial' kernels;
%             undefined fields are zero.
%
% OUTPUT
%   DDF_out:  'struct', the minimal DDF: the fields of DDF, with Cr, Br1,
%             Br2, Drv, Cv, Cvd replaced as above; DDF itself when K = 0.

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
% MMP, 11/21/2026: Changes since PIETOOLS 2025: help and comments revised.


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Construct the DDF representation
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Complete and check the DDF (undefined terms zero), then read off the
% dimensions.
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
 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% First: A0, B1, B2, C1, C2, D11, D12, D21, D22, Bv, D1v, D2v, are
% unchanged from the original DDF.

DDF_out=DDF;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% However, we will re-define:
% Cr{i}, Br1{i}, Br2{i}, Drv{i}, Cv{i}, Cvd{i}

% The first step is to form the vertical concatenation matrices for each
% delay, tau(i)
% No delays: the DDF is returned as initialized.
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
    % The DDF re-formulation is now complete.
    %    convert_PIETOOLS_DDF; % This constructs the PIE representation.
end
end
