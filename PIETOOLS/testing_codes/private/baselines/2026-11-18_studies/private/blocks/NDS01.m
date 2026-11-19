% examples_NDSDDF_library_PIETOOLS.m lines 48-57, one comment level removed (lb_blocks)
 % This is the version used in [1], as well as Han, 2005; Fridman, 2001; Fridman,Shaked
 % 2002, Han, 2004ab, [6] He et al. 2004; Lien and Chen, 2003[3]; Niculescu, 2000;
 stability=1
 stability_dual=1
 
 NDS.A0=[-.9 .2;.1 -.9];
 NDS.Ai{1}=[-1.1 -.2;-.1 -1.1];
 NDS.Ei{1}=[-.2 0;.2 -.1];
 NDS.tau(1)= 2.2; 
%% analytic max is supposedly 2.2255 (Han 2005) although no citation is given
