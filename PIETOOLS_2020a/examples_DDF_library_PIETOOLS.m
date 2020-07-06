%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% examples_DDF_library_PIETOOLS.m     PIETOOLS 2020a
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This library file contains the definition of some common DDF systems
% drawn from the literature. To use one of the examples, simply uncomment
% the variable definitions for that example. Please make sure no other
% example is uncommented, or you will have problems.
%
% The examples are grouped into various problem types. These types are as
% follows:
% 1. Stability Analysis Tests
% 2. Hinf-gain Analysis
% 3. Optimal Control Problems
% 4. Optimal Estimator Design Problems
%
% The examples are typically called using a line in the main PIETOOLS_DDF.m
% file. Simply save your changes to the library file, uncomment the line
% examples_DDF_library_PIETOOLS.m in PIETOOLS_DDF.m and run PIETOOLS_DDF.m.
% Of course, you can also run the library file directly.
%
% Most examples in each problem type include the option to solve the LPI
% according to their problem type. Hence, it is recommended that you also
% uncomment the toggle as well. e.g. stability=1, Hinf_gain=1, etc.
%
% When relevant, we also include citations for each example, indicating the
% sources. The bibtex for each citation is included at the end of the
% library file.
%
% If you wish to include a new example in our library, please send it to us
% and we will include it in the next release. Please also include the
% citation information, if available.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% STABILITY TEST EXAMPLES %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% A pure Difference Equation r_1(t)=Drv{1}*Cv{1}*r_1(t-tau(1))+Drv{2}*Cv{2}*r_2(t-tau(2))
stability=1;
Drv{1}=.5; Cv{1}=1;Drv{2}=.25; Cv{2}=1;
tau(1)=1; tau(2)=2;
% 

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% CONTROLLER SYNTHESIS EXAMPLES
% % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % % % %%% Multiple Showering People - tracking with integral control
% % % This is the DDF implementation
% Hinf_control=1
% n=5; %This is the number of users and can be changed
% shower_ex=2;
% nndelay=n;
% % the delay for user i is tau(i)=i;
% for i=1:nndelay
%     tau(i)=n;
% end
% % set all alphas to 1
% alpha=1;
% alphav=alpha*ones(n,1);
% % set all gammas to 1/n
% gamma=1/n;
% gammaM=gamma*ones(n,n);
% zn=zeros(n);
% 
% A0=[zn eye(n);zn zn];
% 
% Gam=zn;
% for i=1:n
%     for j=1:n
%         if i~=j
%             Gam(i,j)=alphav(j)*gammaM(i,j);
%         else
%             Gam(i,j)=-alphav(i);
%         end
%         
%     end
% end
% B1=[-eye(n);-Gam];
% B2=[zn;eye(n)];
% 
% %There are 3 option for the number of outputs in this example:
% %%%% Full outputs, Full disturbances
% % B1=[-eye(n);-Gam];
% % C1=[eye(n) zn;zn zn];
% % D11=[zn;zn];
% % D12=[zn;.1*eye(n)];
% 
% %%%%% 2 outputs, Full disturbances
% B1=[-eye(n);-Gam];
% C1=[ones(1,n) zeros(1,n);zeros(1,2*n)];
% D11=[zeros(2,n)];
% D12=[zeros(1,n);.1*ones(1,n)];
% % N=5 - light - gam = .714, IPM=17.4
% % N=30 - extreme - gam = 5.37, IPM=35,620
% 
% % Now for new DDF terms
% Bv=[zn;Gam];
% %D1v=zn
% 
% for i=1:nndelay
%     Cr{i}=[zeros(1,n+i-1) 1 zeros(1,n-i)];
% % Br1i=Br2i=Drvi=0;
% Cv{i}=zeros(n,1);Cv{i}(i,1)=1;
% %Cvdi=0;
% end
%%%gam_guess; % min in [.32 .38]
