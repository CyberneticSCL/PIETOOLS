% examples_DDE_library_PIETOOLS.m lines 141-150, one comment level removed (lb_blocks)
 % % % Example A.1 from [3]
    Hinf_gain=1
    Hinf_gain_dual=1
     DDE.A0=[-2 0;0 -.9];%
     DDE.Ai{1}=[-1 0; -1 -1];%
     DDE.B1=[-.5;1];
     DDE.C1=[1 0];
 % %   DDE.D12=[0;.1];
     DDE.tau(1) = .846;
 % %%%% answer - .2364 (tau=.846)
