% examples_DDE_library_PIETOOLS.m lines 181-193, one comment level removed (lb_blocks)
 % % % % Example 3 from [4], taken from [8], [7]
 %
 Hinf_control=1
   DDE.A0=[0 0;0 1];%
   DDE.Ai{1}=[-1 -1; 0 -.9];%
   DDE.B1=[1;1]; DDE.B2=[0;1];
   DDE.C1=[0 1;0 0]; %DDE.C1i{1}=[0 0;0 0];
   %D11=[0;0];
   DDE.D12=[0;.1];
    DDE.tau(1) = 2;
 %
 % % % %%%% answer - .1000 (tau=.99)
 % % % %%%% answer - 1.339 (tau=2)
