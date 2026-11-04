% examples_DDE_library_PIETOOLS.m lines 101-109, one comment level removed (lb_blocks)
 % % % Example D from [2]
  stability=1
  stability_dual=1
   tauN=1.3;
   DDE.A0=[0 1;-1 .1];
   DDE.Ai{1}=[0 0; -1 0];
   DDE.Ai{2}=[0 0;1 0];
   DDE.tau(1)=tauN/2;DDE.tau(2)=tauN;
 % % %%%% answer - tauNmax=1.3722
