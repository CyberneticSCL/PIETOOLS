% examples_DDE_library_PIETOOLS.m lines 127-139, one comment level removed (lb_blocks)
 % Example 1 (unsourced)
 % \dot x=-x(t-tau)+w(t)
 % z(t)=x(t)+w(t)
 Hinf_gain=1
 Hinf_gain_dual=1
  DDE.A0=[0];
  DDE.Ai{1}=[-1];
  DDE.B1=[1];
  DDE.C1=[1];
  DDE.D11=[1];
  DDE.tau(1)=1.5;
 %%%% answer - 2 @ tau=.5
 %%%% answer - 25.907 @ tau=1.5
