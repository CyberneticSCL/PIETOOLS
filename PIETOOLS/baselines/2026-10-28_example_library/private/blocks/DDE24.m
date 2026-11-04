% examples_DDE_library_PIETOOLS.m lines 389-401, one comment level removed (lb_blocks)
 % % % Example 2 from [10]
 Hinf_estimator=1;
    DDE.A0=[0 3;-4 -5];%
    DDE.Ai{1}=[-.1 0; .2 -.2];%
    DDE.Ai{2}=[0 .1; -.2 -.3];%
    DDE.B1=[-.4545 0;0 .9090];
    DDE.C1=[0 100];
    DDE.C2=[0 100];
    DDE.C2i{1}=[0 10];
    DDE.C2i{2}=[0 2];
    DDE.D21=[1 1];
    DDE.tau(1)=.3; DDE.tau(2)=.5;
 %
