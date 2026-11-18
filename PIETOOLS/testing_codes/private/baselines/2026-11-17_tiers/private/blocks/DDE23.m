% examples_DDE_library_PIETOOLS.m lines 377-387, one comment level removed (lb_blocks)
 % % Example 1 from [10]
 Hinf_estimator=1;
    DDE.A0=[0 0;0 1];%
    DDE.Ai{1}=[-1 -1; 0 .9];%
    DDE.B1=[1 0;0 1];
    DDE.C1=[1 0];
     DDE.C2i{1}=[1 10];
     DDE.D21i{1}=[0 5];
     DDE.tau(1)=1;

 % %%%% answer - 1.808
