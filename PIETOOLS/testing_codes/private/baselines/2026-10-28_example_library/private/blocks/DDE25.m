% examples_DDE_library_PIETOOLS.m lines 407-418, one comment level removed (lb_blocks)
 % % Example 1 from [13]
 H2_norm=1;
 H2_norm_dual=1;
  a=1;
  b=2;
  c=2;
  tau=0.5;
 DDE.Ai{1}=-a;
 DDE.B1=b;
 DDE.C1=c;
 DDE.tau = tau;
 % %%%% answer - 3.6724
