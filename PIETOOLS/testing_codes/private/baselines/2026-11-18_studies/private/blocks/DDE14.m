% examples_DDE_library_PIETOOLS.m lines 196-206, one comment level removed (lb_blocks)
 % % % Example B.2 from [4], adapted from [7]
 % % % \dot x=.1x(t)-x(t-tau)+w(t)+u(t)
 % % % % z(t)=x(t)+w(t)
 Hinf_control=1
     DDE.A0=[2 1;0 -1];%
     DDE.Ai{1}=[-1 0; -1 1];%
     DDE.B1=[-.5;1]; DDE.B2=[3;1];
     DDE.C1=[1 -.5;0 0];
     DDE.D12=[0;1];
     DDE.tau(1) = .3;
 %%% answer - .3953
