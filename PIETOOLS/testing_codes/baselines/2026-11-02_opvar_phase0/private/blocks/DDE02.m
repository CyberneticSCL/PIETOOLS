% examples_DDE_library_PIETOOLS.m lines 54-61, one comment level removed (lb_blocks)
 %%% Distributed delay from [1], page 270
 stability=1
 stability_dual=1
 Ap1=[2 2.25;0 -.5];Ap2=[-1 0;0 -1];
 DDE.A0=[-1.5 0;0.5 -1];  DDE.Adi{1}=Ap1-Ap2; DDE.Adi{2}=Ap2;
 tauN=1.8;
 DDE.tau(1)=tauN/2;DDE.tau(2)=tauN;%
 % [1] says the stable range is tauNmax<2
