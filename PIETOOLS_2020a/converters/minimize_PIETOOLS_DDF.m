%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% minimize_PIETOOLS_DDF.m     PIETOOLS 2020a
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
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
% The following dimensions should have been defined in the
% initialize_PIETOOLS_DDE script
% nx=size(A0,1);   % number of states
% nw=size(B1,2);   % number of disturbances
% nu=size(B2,2);   % number of inputs
% nz=size(C10,1);  % number of regulated outputs
% ny=size(C20,1);  % number of sensed outputs
% nv=size(Bv,2);   % dimension of the output from the delayed channels
% nri{i}=size(Cri{i},1) % dimension of the original delayed channel, i
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
%
% Construct the DDF representation
% 
% First: A0, B1, B2, C1, C2, D11, D12, D21, D21, Bv, D1v, D2v, are
% unchanged from the original DDF. 
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
        T=Cv{i}*[Cr{i} Br1{i} Br2{i} Drv{i}];
        [ZthT,HT] = poly_decomp_1D_econ(Cvd{i}'); % so that Cvd{i}'=HT*ZthT or Cvd{i}=ZthT'*HT'
        Tdh=HT'*[Cr{i} Br1{i} Br2{i} Drv{i}];
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
        Cv{i}=Ls;     % permanently re-assign Cv{i}
        Cvd{i}=Zth*Lsd;  % permanently re-assign Cvd{i}
        % At this point, we have the channel r_i=Rs*[x;w;u;v]
        % with output Ls{i}r_i(t-\tau_i)+int_{-\tau_i}^0 Zth{i}*Lsd{i}r_i(t+s)ds
        nr{i}=size(RsF,1);  % permanently re-assign nr{i}
        Cr{i}=RsF(:,1:nx);  % permanently re-assign Cr{i}
        Br1{i}=RsF(:,(nx+1):(nx+nw));   % permanently re-assign Br1{i}
        Br2{i}=RsF(:,(nx+nw+1):(nx+nw+nu));  % permanently re-assign Br2{i}
        Drv{i}=RsF(:,(nx+nw+nu+1):end);  % permanently re-assign Drv{i}
        %       clear Pb Pbd UT DT Vt rT UsT VsT DdT DsT LsT
    end
% The DDF re-formulation is now complete. We can call
% convert_PIETOOLS_DDF (initialize is not necessary)
    convert_PIETOOLS_DDF; % This constructs the PIE representation. 
end
