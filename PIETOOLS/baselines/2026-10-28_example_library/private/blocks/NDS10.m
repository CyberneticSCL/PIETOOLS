% examples_NDSDDF_library_PIETOOLS.m lines 179-185, one comment level removed (lb_blocks)
%%%%% Source Uncertain
 stability=1
 stability_dual=1
 NDS.A0=[-4 -1;0 -3];
 NDS.Ai{1}=[2 0;1 1];
 NDS.Ei{1}=+[.2 0;-.1 -.2];
 NDS.tau(1)=.5;
