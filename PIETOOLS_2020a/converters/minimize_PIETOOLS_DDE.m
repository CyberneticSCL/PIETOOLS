%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% minimize_PIETOOLS_DDE.m     PIETOOLS 2020a
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% This script takes a fully formed DDE representation and find the minimal
% DDF representation of that DDE. 
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
% Convert DDE to a minimal DDF representation
% The following dimensions should have been defined in the
% initialize_PIETOOLS_DDE script
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
%          u(t)];

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
    
    nv=nx+nz+ny;
    for i=1:nK
        Pb=[Ai{i} B1i{i} B2i{i};
            C1i{i}  D11i{i} D12i{i};
            C2i{i}  D21i{i} D22i{i}];
        Pdb=[Adi{i} B1di{i} B2di{i};
            C1di{i}  D11di{i} D12di{i};
            C2di{i}  D21di{i} D22di{i}];
        [ZthT,HT] = poly_decomp_1D_econ(Pdb');
        Zth=ZthT';
        [UT,DT,VT]=svd([Pb;HT'],'econ');
        rT=rank(DT);
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
        Cr{i}=Rs(:,1:nx);
        Br1{i}=Rs(:,(nx+1):(nx+nw));
        Br2{i}=Rs(:,(nx+nw+1):end);
        Drv{i}=zeros(nr{i},nv);
        Cv{i}=Ls;
        Cvd{i}=Zth*Lsd;
        %       clear Pb Pbd UT DT Vt rT UsT VsT DdT DsT LsT
    end
    Bv=[eye(nx) zeros(nx,nz) zeros(nx,ny)];
    D1v=[zeros(nz,nx) eye(nz) zeros(nz,ny)];
    D2v=[zeros(ny,nx) zeros(ny,nz) eye(ny)];
% The DDF formulation is now complete. We can call
% convert_PIETOOLS_DDF (initialize is not necessary)
    convert_PIETOOLS_DDF;
end
