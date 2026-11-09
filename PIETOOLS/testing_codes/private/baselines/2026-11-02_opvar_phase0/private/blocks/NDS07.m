% examples_NDSDDF_library_PIETOOLS.m lines 126-133, one comment level removed (lb_blocks)
 % % This is an example from [7], Zhao, 2018
 stability=1
 stability_dual=1
 NDS.A0=[-1];
 NDS.Ai{2}=[-2];
 NDS.Ei{1}=[.2];
 h=.43; %max at .5? (stability verified for h<.603)
 NDS.tau=[h 2*h];
