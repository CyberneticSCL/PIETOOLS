% examples_NDSDDF_library_PIETOOLS.m lines 103-112, one comment level removed (lb_blocks)
 % This is an example from [7], Zhao, 2018
 stability=1
 stability_dual=1
 NDS.A0=[-2 0;0 -.9];
 NDS.Ai{2}=[-1 0;-1 -1];
 NDS.Ei{1}=[.1 0;0 .1];
 h=1.58; %max at 1.5804? (stability verified for h<2.04)
 NDS.tau=[h 3*h];
 % % dims = 8, 2
 % % cpusecs = 22.56, .332
