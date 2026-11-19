% examples_DDE_library_PIETOOLS.m lines 341-350, one comment level removed (lb_blocks)
 % % % Example 1 from [12], adapted from [11]
 Hinf_estimator=1;
    DDE.A0=[-3 4;2 0];%
    DDE.Ai{1}=[0 0; 1 0];%
    DDE.B1=[1 0;0 1];
    DDE.C1=[1 0;0 1];
    DDE.C2=[0 7];
 
    DDE.tau(1) = .3;
 % % %%%% answer - .2357
