% examples_DDE_library_PIETOOLS.m lines 209-218, one comment level removed (lb_blocks)
 % % Example B.3 from [4], adapted from [6]
 Hinf_control=1
    DDE.A0=[-1 2;0 1];%
    DDE.Ai{1}=[.6 -.4; 0 0];%
    DDE.Ai{2}=[0 0; 0 -.5];%
    DDE.B1=[1;1]; DDE.B2=[0;1];
    DDE.C1=[1 0;0 1;0 0];
    DDE.D12=[0;0;.1];
    DDE.tau(1) = 1;DDE.tau(2)=2;
 %%% answer - .6104
