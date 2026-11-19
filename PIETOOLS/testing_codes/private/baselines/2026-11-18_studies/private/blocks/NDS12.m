% examples_NDSDDF_library_PIETOOLS.m lines 206-215, one comment level removed (lb_blocks)
 %%% A Neutral-Type System: 
 %%% A 2-delay Neutral-Type System: 
 %%% \dot x(t)- \sum_i Ei{i}\dot x(t-\tau(i))=A0 x(t)+\sum_i Ai{i}x(t-\tau_i)
 stability=1
 stability_dual=1
 
 NDS.Ei{1}=.2;NDS.Ei{2}=0;
 NDS.A0=-1;NDS.Ai{1}=0;NDS.Ai{2}=-2;
 h=.4; % max at .4448?
 NDS.tau=[h 2*h];
