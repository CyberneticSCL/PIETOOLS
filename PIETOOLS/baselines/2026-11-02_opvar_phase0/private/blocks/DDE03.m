% examples_DDE_library_PIETOOLS.m lines 64-71, one comment level removed (lb_blocks)
 % % Example B from [1], page 193
 stability=1
 stability_dual=1
 DDE.A0=[0 1;-2 .1];
 DDE.Ai{1}=[0 0;1 0];
 DDE.tau(1)=1.5;
 % %%%% answer - taumin=.100174
 % %%%% answer - taumax=1.71785
