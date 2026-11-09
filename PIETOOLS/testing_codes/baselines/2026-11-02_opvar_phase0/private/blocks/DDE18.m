% examples_DDE_library_PIETOOLS.m lines 278-318, one comment level removed (lb_blocks)
 % % % %%% Multiple Showering People from [2] - tracking with integral control
 % % This is the DDE implementation
 Hinf_control=1
 n=3; %This is the number of users and can be changed
 shower_ex=2;
 nndelay=n;
 for i=1:nndelay
     DDE.tau(i)=n;
 end
 alpha=1;
 alphav=alpha*ones(n,1);
 gamma=1/n;
 gammaM=gamma*ones(n,n);
 zn=zeros(n);
  DDE.A0=[zn eye(n);zn zn];
 for i=1:n
     gammav=gammaM(:,i);
     Atemp=zn;
     Atemp(:,i)=alphav(i)*gammav;
     Atemp(i,i)=-alphav(i);
     DDE.Ai{i}=[zn zn;zn Atemp];
 end
 Gam=zn;
 for i=1:nndelay
     for j=1:n
         if i~=j
             Gam(i,j)=alphav(j)*gammaM(i,j);
         end
     end
 end
 DDE.B1=[-eye(n);diag(alphav)-Gam];
 DDE.B2=[zn;eye(n)];
 DDE.C1=[eye(n) zn;zn zn];
 DDE.D12=[zn;.1*eye(n)];
 % % %gam_guess; % min in [.32 .38]
 % % n=3 -> .4092
 % % n=10 -> .5964

% for n=3: full size: d=36, reduced size: d=3 
% for n=5: full size: d=5, reduced size: d=100 
% for n=10: full size: d=400, reduced size: d=10 
