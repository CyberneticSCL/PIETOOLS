% examples_NDSDDF_library_PIETOOLS.m lines 88-99, one comment level removed (lb_blocks)
 % This is the version used in [2](park and Won, 1999a, Ex. 2.3), as well as [1];
 stability=1
 stability_dual=1
 
 NDS.A0=[-2 1;0 -1];
 NDS.Ai{1}=[0 .3;-.3 0];
 NDS.Ai{2}=[.1 -.05;.05 .1];
 NDS.Ei{1}=[0 -.1;-.1 0];
 NDS.Ei{2}=[.05 0;0 .05];
 NDS.tau(1)= 2.262; 
 NDS.tau(2)= .5;
% Said to be stable for any h1>0, h2<.552
