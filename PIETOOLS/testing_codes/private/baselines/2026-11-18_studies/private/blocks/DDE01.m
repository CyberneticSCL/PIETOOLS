% examples_DDE_library_PIETOOLS.m lines 47-52, one comment level removed (lb_blocks)
 %%% Distributed delay example
 stability=1
 stability_dual=1
 DDE.A0=[0.2 0;0.2 0.1]; DDE.Ai{1}=[0 0;0 0]; DDE.Adi{1}=[-1 0;-1 -1];
 DDE.tau(1)=2.0;%    h=1/DDE.tau(1);
% the analytical range is [0.2, 2.04]
