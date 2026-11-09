% examples_DDE_library_PIETOOLS.m lines 164-174, one comment level removed (lb_blocks)
 % % % % Example D, derived from [2]
 Hinf_gain=1
 Hinf_gain_dual=1
  tauN=1.0;
  DDE.A0=[0 1;-1 .1];
  DDE.Ai{1}=[0 0; -1 0];
  DDE.Ai{2}=[0 0;1 0];
     DDE.B1=[-.5;1];
     DDE.C1=[1 0]; %DDE.C1i{1}=[0 0;0 0];
  DDE.tau(1)=tauN/2;DDE.tau(2)=tauN;
 % % %%%% answer - gamma=4.9784
