% examples_DDE_library_PIETOOLS.m lines 92-99, one comment level removed (lb_blocks)
 % % Example C from [2] and from [3]
 stability=1
 stability_dual=1
 DDE.A0=[-2];
 DDE.Ai{1}=[2.99];
 DDE.Ai{2}=[-1]
 DDE.tau(1)=1;DDE.tau(2)=2;
 % %%%% answer - bmax=3
