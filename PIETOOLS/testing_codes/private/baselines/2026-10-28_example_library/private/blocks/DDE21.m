% examples_DDE_library_PIETOOLS.m lines 353-362, one comment level removed (lb_blocks)
 % Example 2 from [12], adapted from [13]
    Hinf_estimator=1
    DDE.A0=[0 0;0 1];%
    DDE.Ai{1}=[-1 -1; 0 .9];%
    DDE.B1=[1 0;0 1];
    DDE.C1=[1 0];
    DDE.C2=[0 1];
    
    DDE.tau(1) = 1;
 %%%% answer - 2.327
