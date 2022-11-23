%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% examples_DDE_library_PIETOOLS.m     PIETOOLS 2022a
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This library file contains the definition of some common DDE systems
% drawn from the literature. To use one of the examples, simply uncomment
% the variable definitions for that example. Please make sure no other
% example is uncommented, or you will have problems.
%
% The examples are grouped into various problem types. These types are as
% follows:
% 1. Stability Analysis Tests
% 2. Hinf-gain Analysis
% 3. Optimal Control Problems
% 4. Optimal Estimator Design Problems
%
% The examples are typically called using a line in the main PIETOOLS_DDE.m
% file. Simply save your changes to the library file, uncomment the line
% examples_DDE_library_PIETOOLS.m in PIETOOLS_DDE.m and run PIETOOLS_DDE.m.
% Of course, you can also run the library file directly.
%
% Most examples in each problem type include the option to solve the LPI
% according to their problem type. Hence, it is recommended that you also
% uncomment the toggle as well. e.g. stability=1, Hinf_gain=1, etc.
%
% When relevant, we also include citations for each example, indicating the
% sources. The bibtex for each citation is included at the end of the
% library file.
%
% If you wish to include a new example in our library, please send it to us
% and we will include it in the next release. Please also include the
% citation information, if available.
%
% Parts of a DDE: tau DDE.Ai A0 B1 B2 C1 C1i D11 D12 D21 D22 Dop
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% If you modify this code, document all changes carefully and include date
% authorship, and a brief description of modifications
%
% Initial coding MP - 10_01_2020
% updated to DDE data structure 5_3_2021


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% STABILITY TEST EXAMPLES %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %%% Distributed delay example
% stability=1
% stability_dual=1
% DDE.A0=[0.2 0;0.2 0.1]; DDE.Ai{1}=[0 0;0 0]; DDE.Adi{1}=[-1 0;-1 -1];
% DDE.tau(1)=2.0;%    h=1/DDE.tau(1);
% the analytical range is [0.2, 2.04]
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %%% Distributed delay from [1], page 270
% stability=1
% stability_dual=1
% Ap1=[2 2.25;0 -.5];Ap2=[-1 0;0 -1];
% DDE.A0=[-1.5 0;0.5 -1];  DDE.Adi{1}=Ap1-Ap2; DDE.Adi{2}=Ap2;
% tauN=1.8;
% DDE.tau(1)=tauN/2;DDE.tau(2)=tauN;%
% % [1] says the stable range is tauNmax<2
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % % Example B from [1], page 193
% stability=1
% stability_dual=1
% DDE.A0=[0 1;-2 .1];
% DDE.Ai{1}=[0 0;1 0];
% DDE.tau(1)=1.5;
% % %%%% answer - taumin=.100174
% % %%%% answer - taumax=1.71785
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % % % Example from [1], page 192
% stability=1
% stability_dual=1
%    DDE.A0=[-2 0;0 -.9];%
%    DDE.Ai{1}=[-1 0; -1 -1];%
%    DDE.tau(1) = 6;
% % %%%% [1] says stable range is tau < 6.17258
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % % Example A from [2]
% stability=1
% stability_dual=1
% % % \dot x=-x(t-tau)+w(t)
% % % z(t)=x(t)+w(t)
%  DDE.A0=[0];
%  DDE.Ai{1}=[-1];
%  DDE.tau(1)=1.56;
% % %%%% answer - taumax=1.5707
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % % Example C from [2] and from [3]
% stability=1
% stability_dual=1
% DDE.A0=[-2];
% DDE.Ai{1}=[2.99];
% DDE.Ai{2}=[-1]
% DDE.tau(1)=1;DDE.tau(2)=2;
% % %%%% answer - bmax=3
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % % % Example D from [2]
%  stability=1
%  stability_dual=1
%   tauN=1.3;
%   DDE.A0=[0 1;-1 .1];
%   DDE.Ai{1}=[0 0; -1 0];
%   DDE.Ai{2}=[0 0;1 0];
%   DDE.tau(1)=tauN/2;DDE.tau(2)=tauN;
% % % %%%% answer - tauNmax=1.3722
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % % % % Example E from [2], taken from [9]
%  stability=1
%  stability_dual=1
%   tauN=1.37;
%   A=[0 0 1 0; 0 0 0 1;-10 10 0 0;5 -15 0 -.25];
%   B=[0;0;1;0]; C=[1 0 0 0];K=1;
%   DDE.A0=A-B*K*C;
%   DDE.Ai{1}=B*K*C;
%   DDE.tau(1)=1.1;
% % % %%%% answer - tauNmax\in [3.9,3.95] (Note:requires veryheavy settings)
%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Hinf gain %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % Example 1 (unsourced)
% % \dot x=-x(t-tau)+w(t)
% % z(t)=x(t)+w(t)
% Hinf_gain=1
% Hinf_gain_dual=1
%  DDE.A0=[0];
%  DDE.Ai{1}=[-1];
%  DDE.B1=[1];
%  DDE.C1=[1];
%  DDE.D11=[1];
%  DDE.tau(1)=1.5;
% %%%% answer - 2 @ tau=.5
% %%%% answer - 25.907 @ tau=1.5
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % % % Example A.1 from [3]
%    Hinf_gain=1
%    Hinf_gain_dual=1
    DDE.A0=[-2 0;0 -.9];%
    DDE.Ai{1}=[-1 0; -1 -1];%
    DDE.B1=[-.5;1];
    DDE.C1=[1 0];
% %   DDE.D12=[0;.1];
    DDE.tau(1) = .846;
% % %%%% answer - .2364 (tau=.846)
% %
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % % % Example A.2 from [3]
%    Hinf_gain=1
%    Hinf_gain_dual=1
%     DDE.A0=[0 1;-2 .1];%
%     DDE.Ai{1}=[0 0; 1 0];%
%     DDE.B1=[1 0;0 1];
%     DDE.C1=[0 1];
%     DDE.tau(1) = .846;%.846; %1.5707
% % %%%% answer - 2.87 (tau=.846)
% %
% % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % % % % Example D, derived from [2]
% Hinf_gain=1
% Hinf_gain_dual=1
%  tauN=1.0;
%  DDE.A0=[0 1;-1 .1];
%  DDE.Ai{1}=[0 0; -1 0];
%  DDE.Ai{2}=[0 0;1 0];
%     DDE.B1=[-.5;1];
%     DDE.C1=[1 0]; %DDE.C1i{1}=[0 0;0 0];
%  DDE.tau(1)=tauN/2;DDE.tau(2)=tauN;
% % % %%%% answer - gamma=4.9784


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Hinf optimal control EXAMPLES %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %
% % % % % Example 3 from [4], taken from [8], [7]
% %
% Hinf_control=1
%   DDE.A0=[0 0;0 1];%
%   DDE.Ai{1}=[-1 -1; 0 -.9];%
%   DDE.B1=[1;1]; DDE.B2=[0;1];
%   DDE.C1=[0 1;0 0]; %DDE.C1i{1}=[0 0;0 0];
%   %D11=[0;0];
%   DDE.D12=[0;.1];
%    DDE.tau(1) = 2;
% %
% % % % %%%% answer - .1000 (tau=.99)
% % % % %%%% answer - 1.339 (tau=2)

% % % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % % % Example B.2 from [4], adapted from [7]
% % % % \dot x=.1x(t)-x(t-tau)+w(t)+u(t)
% % % % % z(t)=x(t)+w(t)
% Hinf_control=1
%     DDE.A0=[2 1;0 -1];%
%     DDE.Ai{1}=[-1 0; -1 1];%
%     DDE.B1=[-.5;1]; DDE.B2=[3;1];
%     DDE.C1=[1 -.5;0 0];
%     DDE.D12=[0;1];
%     DDE.tau(1) = .3;
% %%% answer - .3953

% % % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % % Example B.3 from [4], adapted from [6]
% Hinf_control=1
%    DDE.A0=[-1 2;0 1];%
%    DDE.Ai{1}=[.6 -.4; 0 0];%
%    DDE.Ai{2}=[0 0; 0 -.5];%
%    DDE.B1=[1;1]; DDE.B2=[0;1];
%    DDE.C1=[1 0;0 1;0 0];
%    DDE.D12=[0;0;.1];
%    DDE.tau(1) = 1;DDE.tau(2)=2;
% %%% answer - .6104
% % % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % % Scalability Example from [4]
% % delaymat=[1 2 3 5 10]
% % dimmat=[1 2 3 5 10]
% %  for iii=5:5 %dimmat=[1 2 3 5 10] n
% %      for jjj=1:5 %delaymat=[1 2 3 5 10] K
% %         clear prog tau Ai A0 B1 B2 C1 C1i D11 D12 D21 D22 Dop
% %  nndim=dimmat(iii);nndelay=delaymat(jjj);
%
% Hinf_control=1
%  nndim=3;nndelay=3;
% 
% DDE.A0=zeros(nndim);
% DDE.B1=ones(nndim,1);
% DDE.B2=ones(nndim,1);
% DDE.C1=[ones(1,nndim);zeros(1,nndim)];
% DDE.D11=[0;0]; DDE.D12=[0;1];
% for j=1:nndelay
%   DDE.Ai{j}=-eye(nndim)/nndelay;
% %  DDE.C1i{j}=zeros(2,nndim);
%   DDE.tau(j) = j/nndelay;
% end
% Answer [3,3] - .9630
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % % % % %%% Spring-Mass Network delay problem
% % First node is referenced to wall. Delay between node i and i+1 is
% % tau(i+1). Control input is force at first node and is delayed by tau(n+1). disturbance input is
% % force at last node. regulated output is position of last node plus
% % control effort. sensed output is position of last node.
% n=3;
% a=1;b=.2;
% DDE.A0=[zeros(n) eye(n); -2*a*eye(n) -b*eye(n)];
% DDE.A0(2*n,n)=-a;DDE.A0(2*n,2*n)=-b;
% DDE.A0(n+1,1)=-a;DDE.A0(n+1,n+1)=-b;
% DDE.Ai{1}=zeros(2*n,2*n);DDE.Ai{1}(n+1,1)=-a;DDE.Ai{1}(n+1,n+1)=-b;
% h=.05;
% DDE.tau(1)=h;
% for i=2:n
%     DDE.Ai{i}=zeros(2*n,2*n);
%     DDE.Ai{i}(n+i,i-1)=a;DDE.Ai{i}(n+i-1,i)=a;DDE.Ai{i}(n+i,n+i-1)=b;DDE.Ai{i}(n+i-1,n+i)=b;
%     DDE.tau(i)=h*i;
% end
% %DDE.tau(n+1)=.1;
% DDE.B2=zeros(2*n,1);B2(n+1,1)=1;
% %DDE.B2i{n+1}=zeros(2*n,1);DDE.B2i{n+1}(n+1,1)=1;
% DDE.B1=zeros(2*n,1);DDE.B1(2*n,1)=1;
% DDE.C1=zeros(1,2*n);DDE.C1(1,n)=1;
% DDE.C2=zeros(1,2*n);DDE.C2(1,n)=1;D12=1;
% stability=1;
% %Hinf_gain=1;
% %Hinf_gain=1;
% %Hinf_estimator=1;
% Hinf_control=1;
% % stable at h=.08, unstable at h=1
% % for n=3: full size: d=24, reduced size: d=5 
% % for n=5: full size: d=60, reduced size: d=9 % minimal comp: 154.95s
% % for n=6: full size: d=84, reduced size: d=11 
% % for n=10: full size: d=220, reduced size: d=19 
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % % % %%% Multiple Showering People from [2] - tracking with integral control
% % % This is the DDE implementation
% Hinf_control=1
% n=3; %This is the number of users and can be changed
% shower_ex=2;
% nndelay=n;
% for i=1:nndelay
%     DDE.tau(i)=n;
% end
% alpha=1;
% alphav=alpha*ones(n,1);
% gamma=1/n;
% gammaM=gamma*ones(n,n);
% zn=zeros(n);
%  DDE.A0=[zn eye(n);zn zn];
% for i=1:n
%     gammav=gammaM(:,i);
%     Atemp=zn;
%     Atemp(:,i)=alphav(i)*gammav;
%     Atemp(i,i)=-alphav(i);
%     DDE.Ai{i}=[zn zn;zn Atemp];
% end
% Gam=zn;
% for i=1:nndelay
%     for j=1:n
%         if i~=j
%             Gam(i,j)=alphav(j)*gammaM(i,j);
%         end
%     end
% end
% DDE.B1=[-eye(n);diag(alphav)-Gam];
% DDE.B2=[zn;eye(n)];
% DDE.C1=[eye(n) zn;zn zn];
% DDE.D12=[zn;.1*eye(n)];
% % % %gam_guess; % min in [.32 .38]
% % % n=3 -> .4092
% % % n=10 -> .5964

% for n=3: full size: d=36, reduced size: d=3 
% for n=5: full size: d=5, reduced size: d=100 
% for n=10: full size: d=400, reduced size: d=10 

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % % % % Wind Tunnel Problem from [4], adapted from [5]
% Hinf_control=1
% aM=1/1.964; kM=-.0117;zM=.8;wM=6;
% DDE.A0=[-aM 0 0; 0 0 1; 0 -wM^2 -2*zM*wM];
%   Ad{1}=[0 kM*aM 0; 0 0 0; 0 0 0];
%  DDE.C0=[1 0 0; 0 1 0;0 0 0]; DDE.C{1}=[0 0 0;0 0 0;0 0 0];
%  DDE.B2=[0;0;wM^2]; %B in paper
%  DDE.B1=[1 0;0 0;0 10]; % D in paper
%  DDE.tau(1)=.33;
%  DDE.D1=[0 0;0 0;0 0];
%  DDE.D2=[0;0;.1];
% % %  gam=1.9369   @6,12, eps1=eps2=.000001, tau=.33 IPM=6.61
% % %  gam=1.964   @4,8, eps1=eps2=.000001, tau=.33 IPM=1.98
% % %  gam=1.964   @2,4, eps1=eps2=.000001, tau=.33 IPM=.743
% % %  % C=[0 10]

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Hinf optimal ESTIMATOR SYNTHESIS EXAMPLES %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % % % Example 1 from [12], adapted from [11]
% Hinf_estimator=1;
%    DDE.A0=[-3 4;2 0];%
%    DDE.Ai{1}=[0 0; 1 0];%
%    DDE.B1=[1 0;0 1];
%    DDE.C1=[1 0;0 1];
%    DDE.C2=[0 7];
% 
%    DDE.tau(1) = .3;
% % % %%%% answer - .2357

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % Example 2 from [12], adapted from [13]
%    Hinf_estimator=1
%    DDE.A0=[0 0;0 1];%
%    DDE.Ai{1}=[-1 -1; 0 .9];%
%    DDE.B1=[1 0;0 1];
%    DDE.C1=[1 0];
%    DDE.C2=[0 1];
% 
%    DDE.tau(1) = 1;
% %%%% answer - 2.327

% % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % Example 3 from [12]
% Hinf_estimator=1
%    DDE.A0=[0 0;0 1];%
%    DDE.Ai{1}=[-1 -1; 0 .9]/2;%
%    DDE.Ai{2}=[-1 -1; 0 .9]/2;%
%    DDE.B1=[1 0;0 1];
%    DDE.C1=[1 0];
%    DDE.C2=[0 1];
% 
%    DDE.tau(1)=.5; DDE.tau(2) = 1;
% %%%% answer - 1.3501
% % % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % % Example 1 from [10]
% Hinf_estimator=1;
%    DDE.A0=[0 0;0 1];%
%    DDE.Ai{1}=[-1 -1; 0 .9];%
%    DDE.B1=[1 0;0 1];
%    DDE.C1=[1 0];
%     DDE.C2i{1}=[1 10];
%     DDE.D21i{1}=[0 5];
%     DDE.tau(1)=1;

% % %%%% answer - 1.808
% % % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % % % Example 2 from [10]
% Hinf_estimator=1;
%    DDE.A0=[0 3;-4 -5];%
%    DDE.Ai{1}=[-.1 0; .2 -.2];%
%    DDE.Ai{2}=[0 .1; -.2 -.3];%
%    DDE.B1=[-.4545 0;0 .9090];
%    DDE.C1=[0 100];
%    DDE.C2=[0 100];
%    DDE.C2i{1}=[0 10];
%    DDE.C2i{2}=[0 2];
%    DDE.D21=[1 1];
%    DDE.tau(1)=.3; DDE.tau(2)=.5;
% % %
% % % %%%% answer - 1.0448
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% RANDOM EXAMPLE GENERATOR %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% seed='default';%234235;
% maxn=3; pvar_deg=1;
% stability=0;
% stability_dual=0;
% Hinf_gain=0;
% Hinf_gain_dual=0;
% Hinf_control=0;
% Hinf_estimator=1;
% 
% %%%%%%%% Begin generating case
% nK=randi(maxn+1)-1; nx =randi(maxn+1)-1; nw=randi(maxn+1)-1; nu=randi(maxn+1)-1; nz=randi(maxn+1)-1; ny=randi(maxn+1)-1;
% %nK=2; nx =2; nw=1; nu=2; nz=1; ny=2;
% rng(seed)
% DDE.tau=rand([1 nK]);
% DDE.A0=rand([nx nx]);
% DDE.B1=rand([nx nw]);
% DDE.B2=rand([nx nu]);
% DDE.C1=rand([nz nx]);
% DDE.C2=rand([ny nx]);
% DDE.D11=rand([nz nw]);
% DDE.D12=rand([nz nu]);
% DDE.D21=rand([ny nw]);
% DDE.D22=rand([ny nu]);
% for i=1:nK
%     DDE.Ai{i}=rand([nx nx]);
%     DDE.B1i{i}=rand([nx nw]);
%     DDE.B2i{i}=rand([nx nu]);
%     DDE.C1i{i}=rand([nz nx]);
%     DDE.C2i{i}=rand([ny nx]);
%     DDE.D11i{i}=rand([nz nw]);
%     DDE.D12i{i}=rand([nz nu]);
%     DDE.D21i{i}=rand([ny nw]);
%     DDE.D22i{i}=rand([ny nu]);
%     DDE.Adi{i}=rand([nx nx]);
%     DDE.B1di{i}=rand([nx nw]);
%     DDE.B2di{i}=rand([nx nu]);
%     DDE.C1di{i}=rand([nz nx]);
%     DDE.C2di{i}=rand([ny nx]);
%     DDE.D11di{i}=rand([nz nw]);
%     DDE.D12di{i}=rand([nz nu]);
%     DDE.D21di{i}=rand([ny nw]);
%     DDE.D22di{i}=rand([ny nu]);
%     for j=1:pvar_deg
%         DDE.Adi{i}=DDE.Adi{i}+rand([nx nx])*s^j;
%         DDE.B1di{i}=DDE.B1di{i}+rand([nx nw])*s^j;
%         DDE.B2di{i}=DDE.B2di{i}+rand([nx nu])*s^j;
%         DDE.C1di{i}=DDE.C1di{i}+rand([nz nx])*s^j;
%         DDE.C2di{i}=DDE.C2di{i}+rand([ny nx])*s^j;
%         DDE.D11di{i}=DDE.D11di{i}+rand([nz nw])*s^j;
%         DDE.D12di{i}=DDE.D12di{i}+rand([nz nu])*s^j;
%         DDE.D21di{i}=DDE.D21di{i}+rand([ny nw])*s^j;
%         DDE.D22di{i}=DDE.D22di{i}+rand([ny nu])*s^j;
%     end
%     
% end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% References %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% [1] -
%@book{gu2003stability,
%  title={Stability of time-delay systems},
%  author={Gu, Keqin and Chen, Jie and Kharitonov, Vladimir L},
%  year={2003},
%  publisher={Springer Science \& Business Media}
%}
%
% [2] -
% @article{peet2018dual,
%   title={A dual to Lyapunov's second method for linear systems with multiple delays and implementation using SOS},
%   author={Peet, Matthew M},
%   journal={IEEE Transactions on Automatic Control},
%   volume={64},
%   number={3},
%   pages={944--959},
%   year={2018},
%   publisher={IEEE}
% }
%
% [3] -
% @article{egorov2014necessary,
%   title={Necessary stability conditions for linear delay systems},
%   author={Egorov, Alexey V and Mondi{\'e}, Sabine},
%   journal={Automatica},
%   volume={50},
%   number={12},
%   pages={3204--3208},
%   year={2014},
%   publisher={Elsevier}
% }
%
% [4] -
% @article{peet2020convex,
%   title={A Convex Solution of the H\_\infty-Optimal Controller Synthesis Problem for Multidelay Systems},
%   author={Peet, Matthew M},
%   journal={SIAM Journal on Control and Optimization},
%   volume={58},
%   number={3},
%   pages={1547--1578},
%   year={2020},
%   publisher={SIAM}
% }
%
% [5] -
% @article{manitius1984feedback,
%   title={Feedback controllers for a wind tunnel model involving a delay: Analytical design and numerical simulation},
%   author={Manitius, A},
%   journal={IEEE Transactions on Automatic Control},
%   volume={29},
%   number={12},
%   pages={1058--1068},
%   year={1984},
%   publisher={IEEE}
% }
%
% [6] -
% @article{cao1998delay,
%   title={Delay-dependent robust stabilization of uncertain systems with multiple state delays},
%   author={Cao, Yong-Yan and Sun, You-Xian and Cheng, Chuwang},
%   journal={IEEE Transactions on Automatic Control},
%   volume={43},
%   number={11},
%   pages={1608--1612},
%   year={1998},
%   publisher={IEEE}
% }
%
% [7] -
% @article{fridman1998h,
%   title={H_\infty-state-feedback control of linear systems with small state delay},
%   author={Fridman, E and Shaked, U},
%   journal={Systems \& control letters},
%   volume={33},
%   number={3},
%   pages={141--150},
%   year={1998},
%   publisher={Elsevier}
% }
%
% [8] -
% @article{li1997criteria,
%   title={Criteria for robust stability and stabilization of uncertain linear systems with state delay},
%   author={Li, Xi and De Souza, Carlos E},
%   journal={Automatica},
%   volume={33},
%   number={9},
%   pages={1657--1662},
%   year={1997},
%   publisher={Elsevier}
% }
%
% [9] -
% @inproceedings{seuret2014complete,
%   title={Complete quadratic Lyapunov functionals using Bessel-Legendre inequality},
%   author={Seuret, Alexandre and Gouaisbaut, Fr{\'e}d{\'e}ric},
%   booktitle={2014 European Control Conference (ECC)},
%   pages={448--453},
%   year={2014},
%   organization={IEEE}
% }
%
% [10] -
% @article{wu2020h,
%   title={H$\backslash$infty-Optimal Observer Design for Linear Systems with Delays in States, Outputs and Disturbances},
%   author={Wu, Shuangshuang and Shivakumar, Sachin and Peet, Matthew M and Hua, Changchun},
%   journal={arXiv preprint arXiv:2004.04482},
%   year={2020}
% }
%
% [11] -
% @incollection{fattouh1998robust,
%   title={Robust observer design for time-delay systems: A riccati equation approach},
%   author={Fattouh, Anas and Sename, Olivier and Dion, Jean-Michel},
%   booktitle={Theory and Practice of Control and Systems},
%   pages={432--436},
%   year={1998},
%   publisher={World Scientific}
% }
%
% [12] -
% @inproceedings{peet2019sos,
%   title={SOS for systems with multiple delays: Part 2. H∞-optimal estimation},
%   author={Peet, Matthew M and Gu, Keqin},
%   booktitle={2019 American Control Conference (ACC)},
%   pages={3870--3876},
%   year={2019},
%   organization={IEEE}
% }
%
% [12] -
% @article{fridman2001new,
%   title={A new H/sub/spl infin//filter design for linear time delay systems},
%   author={Fridman, Emilia and Shaked, Uri},
%   journal={IEEE Transactions on Signal Processing},
%   volume={49},
%   number={11},
%   pages={2839--2843},
%   year={2001},
%   publisher={IEEE}
% }


