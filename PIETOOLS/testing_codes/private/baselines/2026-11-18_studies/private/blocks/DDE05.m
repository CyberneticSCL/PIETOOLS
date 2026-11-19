% examples_DDE_library_PIETOOLS.m lines 82-90, one comment level removed (lb_blocks)
 % % Example A from [2]
 stability=1
 stability_dual=1
 % % \dot x=-x(t-tau)+w(t)
 % % z(t)=x(t)+w(t)
  DDE.A0=[0];
  DDE.Ai{1}=[-1];
  DDE.tau(1)=1.56;
 % %%%% answer - taumax=1.5707
