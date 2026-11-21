%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% convert_PIETOOLS_NDS2DDF.m     PIETOOLS 2022
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function DDF=convert_PIETOOLS_NDS2DDF(NDS)
% This function takes a fully formed NDS representation and finds a
% DDF representation of that NDS.
%
% Each delay tau(i) gets the channel r_i = [x; w; u; \dot x], with
% \dot x = A0 x + B1 w + B2 u + v_x, and the output
%   Cv{i} = [Ai{i}  B1i{i}  B2i{i}  Ei{i};
%            C1i{i} D11i{i} D12i{i} C1ei{i};
%            C2i{i} D21i{i} D22i{i} C2ei{i}],
% where v = [v_x; v_z; v_y] collects the delayed terms of the x, z and y
% equations (Bv, D1v, D2v select them). The kernels Cvd{i} hold the
% distributed delay coefficients (Adi{i}, B1di{i}, ..., C2dei{i}) in the
% same layout.
%
% The formulae were defined in (12) and (15)-(16) in the paper
% ''Representation of Networks and Systems with Delay: DDEs, DDFs, ODE-PDEs and PIEs''
% reference: https://arxiv.org/abs/1910.03881
%
% INPUT
%   NDS:      'struct', the neutral delay system with K = numel(NDS.tau) >= 1 delays
%     \dot x(t) = A0 x(t) + sum_i Ai{i} x(t-tau(i)) + sum_i int_{-tau(i)}^0 Adi{i}(s) x(t+s) ds
%               + sum_i Ei{i} \dot x(t-tau(i)) + sum_i int_{-tau(i)}^0 Edi{i}(s) \dot x(t+s) ds
%               + B1 w(t) + (B1i, B1di likewise) + B2 u(t) + (B2i, B2di likewise)
%          z(t) = C1 x(t) + D11 w(t) + D12 u(t) + (C1i, C1di, C1ei, C1dei, D11i, ..., D12di likewise)
%          y(t) = C2 x(t) + D21 w(t) + D22 u(t) + (C2i, C2di, C2ei, C2dei, D21i, ..., D22di likewise)
%             undefined fields are zero.
%
% OUTPUT
%   DDF:      'struct', the DDF
%     \dot x(t) = A0 x(t) + B1 w(t) + B2 u(t) + Bv v(t)
%          z(t) = C1 x(t) + D11 w(t) + D12 u(t) + D1v v(t)
%          y(t) = C2 x(t) + D21 w(t) + D22 u(t) + D2v v(t)
%        r_i(t) = Cr{i} x(t) + Br1{i} w(t) + Br2{i} u(t) + Drv{i} v(t)
%          v(t) = sum_i Cv{i} r_i(t-tau(i)) + sum_i int_{-tau(i)}^0 Cvd{i}(s) r_i(t+s) ds
%             fields tau, A0, B1, B2, C1, C2, D11, D12, D21, D22 of the
%             NDS, Bv, D1v, D2v, and 1 x K cells Cr, Br1, Br2, Drv, Cv, Cvd.
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
% For support, contact M. Peet, Arizona State University at mpeet@asu.edu
% Initial coding MP - 10_01_2020
% MMP, 11/21/2026: Changes since PIETOOLS 2025: the distributed delay kernels
%   are stored in the field Cvd of the output; an NDS without delays (tau
%   empty) throws an error; help and comments revised.


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
% First: A0, B1, B2, C1, D11, D12, C2, D21, D22 are unchanged from the NDS
DDF.tau=NDS.tau;DDF.A0=NDS.A0;DDF.B1=NDS.B1;DDF.B2=NDS.B2;DDF.C1=NDS.C1;DDF.C2=NDS.C2;DDF.D11=NDS.D11;DDF.D12=NDS.D12;DDF.D21=NDS.D21;DDF.D22=NDS.D22;

% However, we must now define:
% Bv, D1v, D2v, Cr{i}, Br1{i}, Br2{i}, Drv{i}, Cv{i}, Cvd{i} 

% The first step is to form the vertical concatenation matrices for each
% delay, tau(i)
% No delays: not supported, an error                                        % MMP, 11/21/2026
if nK==0
%   DDF=convert_PIETOOLS_DDE(NDS);                                          % MMP, 11/21/2026 (was)
    error('The NDS has no delays (tau is empty). Define at least one delay.') % MMP, 11/21/2026
else

    nv=nx+nz+ny; nt=nx+nw+nu;
    for i=1:nK
        % Channel r_i = [x; w; u; \dot x]; discrete delay coefficients of
        % [x; w; u; \dot x] in the x, z, y equations
        nr{i}=nt+nx;
        DDF.Cv{i}=[NDS.Ai{i} NDS.B1i{i} NDS.B2i{i} NDS.Ei{i};
               NDS.C1i{i}  NDS.D11i{i} NDS.D12i{i} NDS.C1ei{i};
               NDS.C2i{i}  NDS.D21i{i} NDS.D22i{i} NDS.C2ei{i}];
        % Distributed delay kernels of [x; w; u; \dot x], layout of Cv{i}   % MMP, 11/21/2026
%       Cvd{i}=[NDS.Adi{i} NDS.B1di{i} NDS.B2di{i} NDS.Edi{i};              % MMP, 11/21/2026 (was)
        DDF.Cvd{i}=[NDS.Adi{i} NDS.B1di{i} NDS.B2di{i} NDS.Edi{i};          % MMP, 11/21/2026
            NDS.C1di{i}  NDS.D11di{i} NDS.D12di{i} NDS.C1dei{i};
            NDS.C2di{i}  NDS.D21di{i} NDS.D22di{i} NDS.C2dei{i}];
        % r_i in terms of x, w, u, v: \dot x = A0 x + B1 w + B2 u + v_x
        DDF.Cr{i}=[eye(nx);zeros(nw,nx);zeros(nu,nx);NDS.A0];
        DDF.Br1{i}=[zeros(nx,nw);eye(nw);zeros(nu,nw);NDS.B1];
        DDF.Br2{i}=[zeros(nx,nu);zeros(nw,nu);eye(nu);NDS.B2];
        DDF.Drv{i}=[zeros(nt,nx) zeros(nt,(nz+ny));eye(nx) zeros(nx,(nz+ny))];
    end
    % v = [v_x; v_z; v_y]: Bv, D1v, D2v select the delayed terms of the
    % x, z and y equations
    DDF.Bv=[eye(nx) zeros(nx,nz) zeros(nx,ny)];
    DDF.D1v=[zeros(nz,nx) eye(nz) zeros(nz,ny)];
    DDF.D2v=[zeros(ny,nx) zeros(ny,nz) eye(ny)];
% The DDF formulation is now complete; each channel has dimension
% 2*nx+nw+nu (no rank reduction).

DDF_out=DDF;
end
