% examples_DDE_library_PIETOOLS.m lines 321-335, one comment level removed (lb_blocks)
 % % % % Wind Tunnel Problem from [4], adapted from [5]
 Hinf_control=1
 aM=1/1.964; kM=-.0117;zM=.8;wM=6;
 DDE.A0=[-aM 0 0; 0 0 1; 0 -wM^2 -2*zM*wM];
  DDE.Adi{1}=[0 kM*aM 0; 0 0 0; 0 0 0];
  DDE.C0=[1 0 0; 0 1 0;0 0 0]; DDE.C{1}=[0 0 0;0 0 0;0 0 0];
  DDE.B2=[0;0;wM^2]; %B in paper
  DDE.B1=[1 0;0 0;0 10]; % D in paper
  DDE.tau(1)=.33;
  DDE.D1=[0 0;0 0;0 0];
  DDE.D2=[0;0;.1];
 % %  gam=1.9369   @6,12, eps1=eps2=.000001, tau=.33 IPM=6.61
 % %  gam=1.964   @4,8, eps1=eps2=.000001, tau=.33 IPM=1.98
 % %  gam=1.964   @2,4, eps1=eps2=.000001, tau=.33 IPM=.743
 % %  % C=[0 10]
