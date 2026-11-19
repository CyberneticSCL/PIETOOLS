% examples_NDSDDF_library_PIETOOLS.m lines 239-243, one comment level removed (lb_blocks)
 %%% A pure Difference Equation r_1(t)=Drv{1}*Cv{1}*r_1(t-tau(1))+Drv{2}*Cv{2}*r_2(t-tau(2))
 stability=1;
 stability_dual=1;
 DDF.Drv{1}=.5; DDF.Cv{1}=1;DDF.Drv{2}=.25; DDF.Cv{2}=1;
 DDF.tau(1)=1; DDF.tau(2)=2;
