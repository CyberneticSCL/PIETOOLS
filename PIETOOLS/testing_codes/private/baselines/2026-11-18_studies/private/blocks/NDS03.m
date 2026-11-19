% examples_NDSDDF_library_PIETOOLS.m lines 73-84, one comment level removed (lb_blocks)
 % This is the version used in [2] (park and Won, 1999a, Ex. 2.2), as well as [1]; 
 stability=1
 stability_dual=1
 
 NDS.A0=[-2 .1;0 -2];
 NDS.Ai{1}=[.1 .1;-.1 .2];
 NDS.Ai{2}=[.1 0;0 .1];
 NDS.Ei{1}=[.1 .05;.02 .1];
 NDS.Ei{2}=[.05 0;0 .05];
 NDS.tau(1)= .3; 
 NDS.tau(2)= .6;
% Said to be stable for any tau1=.3, tau2=.6 (verified by pietools)
