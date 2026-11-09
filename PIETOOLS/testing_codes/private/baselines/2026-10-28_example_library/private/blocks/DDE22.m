% examples_DDE_library_PIETOOLS.m lines 365-375, one comment level removed (lb_blocks)
 % Example 3 from [12]
 Hinf_estimator=1
    DDE.A0=[0 0;0 1];%
    DDE.Ai{1}=[-1 -1; 0 .9]/2;%
    DDE.Ai{2}=[-1 -1; 0 .9]/2;%
    DDE.B1=[1 0;0 1];
    DDE.C1=[1 0];
    DDE.C2=[0 1];
 
    DDE.tau(1)=.5; DDE.tau(2) = 1;
 %%%% answer - 1.3501
