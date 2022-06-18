%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% minimize_PIETOOLS_DDE.m     PIETOOLS 2021b
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function DDF_out=minimize_PIETOOLS_DDE2DDF(DDE)
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
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% If you modify this code, document all changes carefully and include date
% authorship, and a brief description of modifications
%
% Initial coding MP - 10_01_2020
%  MP - 5_30_2021; converted script to function format

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Convert DDE to a minimal DDF representation
% The following dimensions assume the DDE has been initialized, so lets
% make sure.
DDE=initialize_PIETOOLS_DDE(DDE);
nK=length(DDE.tau);   % number of states
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
%%% ADD case logic for when there are no delays
if nK==0
    DDE=convert_PIETOOLS_DDE(DDE);
else
    % First: tau, A0, B1, B2, C1, D11, D12, C2, D21, D21 are unchanged from the DDE
    DDF.tau=DDE.tau;DDF.A0=DDE.A0;DDF.B1=DDE.B1;DDF.B2=DDE.B2;DDF.C1=DDE.C1;DDF.C2=DDE.C2;DDF.D11=DDE.D11;DDF.D12=DDE.D12;DDF.D21=DDE.D21;DDF.D22=DDE.D22;
    % However, we must now define:
    % Bv, D1v, D2v, Cr{i}, Br1{i}, Br2{i}, Drv{i}, Cv{i}, Cvd{i}
    
    nv=nx+nz+ny;
    for i=1:nK
        Pb=[DDE.Ai{i} DDE.B1i{i} DDE.B2i{i};
            DDE.C1i{i}  DDE.D11i{i} DDE.D12i{i};
            DDE.C2i{i}  DDE.D21i{i} DDE.D22i{i}];
        Pdb=[DDE.Adi{i} DDE.B1di{i} DDE.B2di{i};
            DDE.C1di{i}  DDE.D11di{i} DDE.D12di{i};
            DDE.C2di{i}  DDE.D21di{i} DDE.D22di{i}];
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
        DDF.Cr{i}=Rs(:,1:nx);
        DDF.Br1{i}=Rs(:,(nx+1):(nx+nw));
        DDF.Br2{i}=Rs(:,(nx+nw+1):end);
        DDF.Drv{i}=zeros(nr{i},nv);
        DDF.Cv{i}=Ls;
        DDF.Cvd{i}=Zth*Lsd;
        %       clear Pb Pbd UT DT Vt rT UsT VsT DdT DsT LsT
    end
    DDF.Bv=[eye(nx) zeros(nx,nz) zeros(nx,ny)];
    DDF.D1v=[zeros(nz,nx) eye(nz) zeros(nz,ny)];
    DDF.D2v=[zeros(ny,nx) zeros(ny,nz) eye(ny)];
    % The DDF formulation is now complete. We can call
    % convert_PIETOOLS_DDF (initialize is not necessary)
    %    PIE=convert_PIETOOLS_DDF(DDF);
    DDF_out=DDF;
end
end