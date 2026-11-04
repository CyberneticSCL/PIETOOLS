% examples_NDSDDF_library_PIETOOLS.m lines 61-69, one comment level removed (lb_blocks)
 % % This is the version used in [4] - Park and Won, 2000, Ex. 3
 stability=1
 stability_dual=1
 
 NDS.A0=[-2 0;0 -1];
 NDS.Ai{1}=[0 .5;.5 0];
 NDS.Ei{1}=[.2 0;0 .2];
  NDS.tau(1)= 20; 
% supposedly stable for tau<.7516 (verified by pietools. indeed, stable out to at least tau=20)
