% examples_DDE_library_PIETOOLS.m lines 243-276, one comment level removed (lb_blocks)
 % % % % %%% Spring-Mass Network delay problem
 % First node is referenced to wall. Delay between node i and i+1 is
 % tau(i+1). Control input is force at first node and is delayed by tau(n+1). disturbance input is
 % force at last node. regulated output is position of last node plus
 % control effort. sensed output is position of last node.
 n=3;
 a=1;b=.2;
 DDE.A0=[zeros(n) eye(n); -2*a*eye(n) -b*eye(n)];
 DDE.A0(2*n,n)=-a;DDE.A0(2*n,2*n)=-b;
 DDE.A0(n+1,1)=-a;DDE.A0(n+1,n+1)=-b;
 DDE.Ai{1}=zeros(2*n,2*n);DDE.Ai{1}(n+1,1)=-a;DDE.Ai{1}(n+1,n+1)=-b;
 h=.05;
 DDE.tau(1)=h;
 for i=2:n
     DDE.Ai{i}=zeros(2*n,2*n);
     DDE.Ai{i}(n+i,i-1)=a;DDE.Ai{i}(n+i-1,i)=a;DDE.Ai{i}(n+i,n+i-1)=b;DDE.Ai{i}(n+i-1,n+i)=b;
     DDE.tau(i)=h*i;
 end
 %DDE.tau(n+1)=.1;
 DDE.B2=zeros(2*n,1);B2(n+1,1)=1;
 %DDE.B2i{n+1}=zeros(2*n,1);DDE.B2i{n+1}(n+1,1)=1;
 DDE.B1=zeros(2*n,1);DDE.B1(2*n,1)=1;
 DDE.C1=zeros(1,2*n);DDE.C1(1,n)=1;
 DDE.C2=zeros(1,2*n);DDE.C2(1,n)=1;D12=1;
 stability=1;
 %Hinf_gain=1;
 %Hinf_gain=1;
 %Hinf_estimator=1;
 Hinf_control=1;
 % stable at h=.08, unstable at h=1
 % for n=3: full size: d=24, reduced size: d=5 
 % for n=5: full size: d=60, reduced size: d=9 % minimal comp: 154.95s
 % for n=6: full size: d=84, reduced size: d=11 
 % for n=10: full size: d=220, reduced size: d=19 
