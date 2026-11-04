% examples_DDE_library_PIETOOLS.m lines 111-120, one comment level removed (lb_blocks)
 % % % % Example E from [2], taken from [9]
  stability=1
  stability_dual=1
   tauN=1.37;
   A=[0 0 1 0; 0 0 0 1;-10 10 0 0;5 -15 0 -.25];
   B=[0;0;1;0]; C=[1 0 0 0];K=1;
   DDE.A0=A-B*K*C;
   DDE.Ai{1}=B*K*C;
   DDE.tau(1)=1.1;
 % % %%%% answer - tauNmax\in [3.9,3.95] (Note:requires veryheavy settings)
