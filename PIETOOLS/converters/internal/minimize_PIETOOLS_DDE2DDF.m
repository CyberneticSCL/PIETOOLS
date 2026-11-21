%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% minimize_PIETOOLS_DDE2DDF.m     PIETOOLS 2022a
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function DDF_out=minimize_PIETOOLS_DDE2DDF(DDE)
%
% This function takes a fully formed DDE representation and finds the minimal
% DDF representation of that DDE.
%
% For each delay tau(i), let Pb be the discrete and Pdb(s) the distributed
% delay coefficients of tau(i) (rows: the x, z, y equations; columns: x,
% w, u), Pdb(s) = Zth(s)*H with Zth(s) a monomial basis, and
% [Pb; H] = [Ls; Lsd]*Rs a rank factorization (SVD). The delay channel is
%   r_i(t) = Rs*[x(t); w(t); u(t)],    nr{i} = rank([Pb; H]),
% with output Ls*r_i(t-tau(i)) + int_{-tau(i)}^0 Zth(s)*Lsd*r_i(t+s) ds;
% v = [v_x; v_z; v_y] collects these outputs for the x, z and y equations.
%
% The formulae were defined in (12) and (15)-(16) in the paper
% ''Representation of Networks and Systems with Delay: DDEs, DDFs, ODE-PDEs and PIEs''
% reference: https://arxiv.org/abs/1910.03881
%
% INPUT
%   DDE:      'struct', the DDE with K = numel(DDE.tau) >= 1 delays
%     \dot x(t) = A0 x(t) + sum_i Ai{i} x(t-tau(i)) + sum_i int_{-tau(i)}^0 Adi{i}(s) x(t+s) ds
%               + B1 w(t) + (B1i, B1di likewise) + B2 u(t) + (B2i, B2di likewise)
%          z(t) = C1 x(t) + D11 w(t) + D12 u(t) + (C1i, C1di, D11i, D11di, D12i, D12di likewise)
%          y(t) = C2 x(t) + D21 w(t) + D22 u(t) + (C2i, C2di, D21i, D21di, D22i, D22di likewise)
%             with 1 x K cells of 'double' matrices (Ai, ...) and of
%             'double' or 'polynomial' kernels (Adi, ...); undefined fields
%             are zero.
%
% OUTPUT
%   DDF_out:  'struct', the minimal DDF
%     \dot x(t) = A0 x(t) + B1 w(t) + B2 u(t) + Bv v(t)
%          z(t) = C1 x(t) + D11 w(t) + D12 u(t) + D1v v(t)
%          y(t) = C2 x(t) + D21 w(t) + D22 u(t) + D2v v(t)
%        r_i(t) = Cr{i} x(t) + Br1{i} w(t) + Br2{i} u(t) + Drv{i} v(t)
%          v(t) = sum_i Cv{i} r_i(t-tau(i)) + sum_i int_{-tau(i)}^0 Cvd{i}(s) r_i(t+s) ds
%             with tau, A0, ..., D22 of the DDE, [Cr{i} Br1{i} Br2{i}] = Rs,
%             Drv{i} = 0, Cv{i} = Ls, Cvd{i} = Zth*Lsd, and Bv, D1v, D2v
%             the selectors of v_x, v_z, v_y.

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
% MMP, 11/19/2026: Release comment pass (comment-only).
% MMP, 11/20/2026: No delays (tau empty): an error, not an unbounded recursion.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Convert DDE to a minimal DDF representation
% The following dimensions assume the DDE has been initialized, so lets
% make sure.
DDE=initialize_PIETOOLS_DDE(DDE);
nK=length(DDE.tau);   % number of delays
nx=size(DDE.A0,1);   % number of states
nw=size(DDE.B1,2);   % number of disturbances
nu=size(DDE.B2,2);   % number of inputs
nz=size(DDE.C1,1);  % number of regulated outputs
ny=size(DDE.C2,1);  % number of sensed outputs
%
%For each delay in tau, the converter will create a delayed channel
%consisting of
%
% ri(t)=Rs[x(t);
%          w(t);
%          u(t)];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% If you modify this code, document all changes carefully and include date
% authorship, and a brief description of modifications
%
% Initial coding MP - 10_01_2020


% Construct the DDF representation
%

% The first step is to form the vertical concatenation matrices for each
% delay, tau(i)
% No delays: not supported; an error, since convert_PIETOOLS_DDE calls      % MMP, 11/20/2026
% this function again (unbounded recursion).                                % MMP, 11/20/2026
if nK==0
%   DDE=convert_PIETOOLS_DDE(DDE);                                          % MMP, 11/20/2026 (was)
    error('The DDE has no delays (tau is empty). Define at least one delay.') % MMP, 11/20/2026
else
    % First: tau, A0, B1, B2, C1, D11, D12, C2, D21, D22 are unchanged from the DDE
    DDF.tau=DDE.tau;DDF.A0=DDE.A0;DDF.B1=DDE.B1;DDF.B2=DDE.B2;DDF.C1=DDE.C1;DDF.C2=DDE.C2;DDF.D11=DDE.D11;DDF.D12=DDE.D12;DDF.D21=DDE.D21;DDF.D22=DDE.D22;
    % However, we must now define:
    % Bv, D1v, D2v, Cr{i}, Br1{i}, Br2{i}, Drv{i}, Cv{i}, Cvd{i}
    
    nv=nx+nz+ny;
    for i=1:nK
        % Discrete (Pb) and distributed (Pdb) delay coefficients of tau(i):
        % rows the x, z, y equations, columns x, w, u
        Pb=[DDE.Ai{i} DDE.B1i{i} DDE.B2i{i};
            DDE.C1i{i}  DDE.D11i{i} DDE.D12i{i};
            DDE.C2i{i}  DDE.D21i{i} DDE.D22i{i}];
        Pdb=[DDE.Adi{i} DDE.B1di{i} DDE.B2di{i};
            DDE.C1di{i}  DDE.D11di{i} DDE.D12di{i};
            DDE.C2di{i}  DDE.D21di{i} DDE.D22di{i}];
        % Pdb(s) = Zth(s)*HT', Zth(s) a monomial basis
        [ZthT,HT] = poly_decomp_1D_econ(Pdb');
        Zth=ZthT';
        % Rank rT of [Pb; HT'] from its SVD; error if delay i has no terms
        [UT,DT,VT]=svd([Pb;HT'],'econ');
        rT=rank(DT);
        if rT==0
            error(['Delay ', int2str(i), ' has no associated coefficients. Remove the delay or correct the coefficients.'])
        end
        % Rank factorization [Pb; HT'] = LsT*Rs, LsT = [Ls; Lsd]
        UsT=UT(:,1:rT);
        VsT=VT(:,1:rT);
        DdT=diag(DT);
        DsT=diag(DdT(1:rT));
        LsT=UsT*DsT;
        Rs=VsT';
        Ls=LsT(1:(nx+nz+ny),:);
        Lsd=LsT((nx+nz+ny+1):end,:);
        % At this point, we have the channel r_i=Rs*[x;w;u]
        % with output Ls{i}r_i(t-\tau_i)+int_{-\tau_i}^0 Zth{i}*Lsd{i}r_i(t+s)ds
        nr{i}=size(Rs,1);
        DDF.Cr{i}=Rs(:,1:nx);
        DDF.Br1{i}=Rs(:,(nx+1):(nx+nw));
        DDF.Br2{i}=Rs(:,(nx+nw+1):end);
        DDF.Drv{i}=zeros(nr{i},nv);
        DDF.Cv{i}=Ls;
        DDF.Cvd{i}=Zth*Lsd;
        %       clear Pb Pbd UT DT Vt rT UsT VsT DdT DsT LsT
    end
    % v = [v_x; v_z; v_y]: Bv, D1v, D2v select the delayed terms of the
    % x, z and y equations
    DDF.Bv=[eye(nx) zeros(nx,nz) zeros(nx,ny)];
    DDF.D1v=[zeros(nz,nx) eye(nz) zeros(nz,ny)];
    DDF.D2v=[zeros(ny,nx) zeros(ny,nz) eye(ny)];
    % The DDF formulation is now complete.
    %    PIE=convert_PIETOOLS_DDF(DDF);
    DDF_out=DDF;
end
end