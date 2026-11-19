% examples_DDE_library_PIETOOLS.m lines 153-161, one comment level removed (lb_blocks)
 % % % Example A.2 from [3]
    Hinf_gain=1
    Hinf_gain_dual=1
     DDE.A0=[0 1;-2 .1];%
     DDE.Ai{1}=[0 0; 1 0];%
     DDE.B1=[1 0;0 1];
     DDE.C1=[0 1];
     DDE.tau(1) = .846;%.846; %1.5707
 % %%%% answer - 2.87 (tau=.846)
