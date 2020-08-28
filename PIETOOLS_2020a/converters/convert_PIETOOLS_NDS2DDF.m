%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% convert_PIETOOLS_NDS2DDF.m     PIETOOLS 2020a
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
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
% The following dimensions should have been defined in the
% initialize_PIETOOLS_NDS script
% nx=size(A0,1);   % number of states
% nw=size(B1,2);   % number of disturbances
% nu=size(B2,2);   % number of inputs
% nz=size(C10,1);  % number of regulated outputs
% ny=size(C20,1);  % number of sensed outputs
%
%For each delay in tau, the converter will create a delayed channel
%consisting of 
%
% ri(t)=Rs[x(t);
%          w(t);
%          u(t);
%          \dot x(t)];

% Construct the DDF representation
% 
% First: A0, B1, B2, C1, D11, D12, C2, D21, D21 are unchanged from the DDE
% However, we must now define:
% Bv, D1v, D2v, Cr{i}, Br1{i}, Br2{i}, Drv{i}, Cv{i}, Cvd{i} 

% The first step is to form the vertical concatenation matrices for each
% delay, tau(i)
%%% ADD case logic for nK=0
if nK==0
    convert_PIETOOLS_DDE;
else
    
    nv=nx+nz+ny; nt=nx+nw+nu;
    for i=1:nK
        nr{i}=nt+nx;
        Cv{i}=[Ai{i} B1i{i} B2i{i} Ei{i};
               C1i{i}  D11i{i} D12i{i} C1ei{i};
               C2i{i}  D21i{i} D22i{i} C2ei{i}];
        Cvd{i}=[Adi{i} B1di{i} B2di{i} Edi{i};
            C1di{i}  D11di{i} D12di{i} C1dei{i};
            C2di{i}  D21di{i} D22di{i} C2dei{i}];
        Cr{i}=[eye(nx);zeros(nw,nx);zeros(nu,nx);A0];
        Br1{i}=[zeros(nx,nw);eye(nw);zeros(nu,nw);B1];
        Br2{i}=[zeros(nx,nu);zeros(nw,nu);eye(nu);B2];
        Drv{i}=[zeros(nt,nx) zeros(nt,(nz+ny));eye(nx) zeros(nx,(nz+ny))];
    end
    Bv=[eye(nx) zeros(nx,nz) zeros(nx,ny)];
    D1v=[zeros(nz,nx) eye(nz) zeros(nz,ny)];
    D2v=[zeros(ny,nx) zeros(ny,nz) eye(ny)];
% The DDF formulation is now complete. We can call
% convert_PIETOOLS_DDF (initialize is not necessary)
% We should probably also call minimize_PIETOOLS_DDF to reduce the state
% dimension
end
