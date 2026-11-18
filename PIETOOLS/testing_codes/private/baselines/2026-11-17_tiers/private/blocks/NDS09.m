% examples_NDSDDF_library_PIETOOLS.m lines 157-176, one comment level removed (lb_blocks)
 % An Example from Ramirez and Sipahi, 2020
 stability=1
 stability_dual=1
 NDS.A0=[-2 .2 -.3 0 -.4;
     .2 -3.8 0 .7 0;
     .8 0 -1.6 0 0;
     0 .8 -.6 -2 .3;
     -1 -.1 -1.5 0 -1.8];
 NDS.Ai{1}=[-2.2 0 0 1 0;
     1.6 -2.2 1.6 0 0;
     -0.2 -0.2 -0.2 -0.2 -0.2;
     0 0.4 -1.4 -3.4 1;
     -0.2 0.4 -0.1 -1.1 -3.3];
 NDS.Ei{1}=[ 0.40888 0.00888 0.20888 -0.09112 -0.29112;
     0 0.2 0 0 0.6;
     -0.1 -0.4 0 -0.8 0;
     0 0 -0.1 0 0;
     0 0 0 -0.2 -0.1];
  NDS.tau(1)=.60; 
% stable for tau <.603
