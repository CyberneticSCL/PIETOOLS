% examples_NDSDDF_library_PIETOOLS.m lines 138-154, one comment level removed (lb_blocks)
 % An Example from Olgac and Sipahi, 2004
 stability=1
 stability_dual=1
 NDS.A0=[0 1 -1 0; 
 -3.346 -2.715 2.075 -2.007;
 -4 0 2 0;
 -3 0 0 6];
 NDS.Ai{1}=[-1 2 2 -1;
     3 3 -2 0;
     1 2 -1 1;
     2 3 1 -3];
 NDS.Ei{1}=[.2 -.1 .5 -.1 ;
     -.3 .09 -.15 -.027;
     -3.333 .1 .2 1;
     -1 2 .5 1];
   NDS.tau(1)=.0275; 
% stable for tau <.02755 
