% examples_NDSDDF_library_PIETOOLS.m showering block (active in the library), verbatim (11/19/2026: option 2 output; tau(i) = i)
% % % % %%% Multiple Showering People - tracking with integral control
% % % This is the DDF implementation
Hinf_control=1
n=2; %This is the number of users and can be changed
shower_ex=2;
nndelay=n;
% the delay for user i is tau(i)=i;
for i=1:nndelay
%   DDF.tau(i)=n;                                                           % MMP, 11/19/2026 (was)
    DDF.tau(i)=i;   % as the comment above and the source, eq. (15)         % MMP, 11/19/2026
end
% set all alphas to 1
alpha=1;
alphav=alpha*ones(n,1);
% set all gammas to 1/n
gamma=1/n;
gammaM=gamma*ones(n,n);
zn=zeros(n);

DDF.A0=[zn eye(n);zn zn];

Gam=zn;
for i=1:n
    for j=1:n
        if i~=j
            Gam(i,j)=alphav(j)*gammaM(i,j);
        else
            Gam(i,j)=-alphav(i);
        end

    end
end
DDF.B1=[-eye(n);-Gam];
DDF.B2=[zn;eye(n)];

% %There are 3 option for the number of outputs in this example:
% %%%% 1) Full outputs, Full disturbances
% DDF.B1=[-eye(n);-Gam];
% DDF.C1=[eye(n) zn;zn zn];
% DDF.D12=[zn;.1*eye(n)];

% %%%%% 2) 2 outputs, Full disturbances
% DDF.B1=[-eye(n);-Gam];
% Option 2 is the regulated output of the source: z = [sum_i T1i;           % MMP, 11/19/2026
% .1*sum_i ui], eq. (15) of Peet, Automatica 2021 (arXiv 1910.03881).       % MMP, 11/19/2026
% Active: with every option commented out, Hinf_control had no z.           % MMP, 11/19/2026
% The source sets tau(i) = i, as the loop above does since 11/19/2026.      % MMP, 11/19/2026
% DDF.C1=[ones(1,n) zeros(1,n);zeros(1,2*n)];                               % MMP, 11/19/2026 (was)
DDF.C1=[ones(1,n) zeros(1,n);zeros(1,2*n)];                                 % MMP, 11/19/2026
% DDF.D11=[zeros(2,n)];                                                     % MMP, 11/19/2026 (was)
DDF.D11=[zeros(2,n)];                                                       % MMP, 11/19/2026
% DDF.D12=[zeros(1,n);.1*ones(1,n)];                                        % MMP, 11/19/2026 (was)
DDF.D12=[zeros(1,n);.1*ones(1,n)];                                          % MMP, 11/19/2026
% % N=5 - light - gam = .714, IPM=17.4
% % N=30 - extreme - gam = 5.37, IPM=35,620
% 
% % 3) Now for new DDF terms
DDF.Bv=[zn;Gam];

for i=1:nndelay
    DDF.Cr{i}=[zeros(1,n+i-1) 1 zeros(1,n-i)];
DDF.Cv{i}=zeros(n,1);DDF.Cv{i}(i,1)=1;
end
% st: example-specific settings (11/19/2026); closed-loop gain              % MMP, 11/19/2026
%     0.439636, 4,673 decision variables (heavy as it is: 0.440633 at       % MMP, 11/19/2026
%     8,465)                                                                % MMP, 11/19/2026
st = lpisettings('heavy');                                                  % MMP, 11/19/2026
st.lpi.pos.len_sos{1} = [4;2;2];                                            % MMP, 11/19/2026
st.lpi.ineq.cap_inc = 'auto';                                               % MMP, 11/19/2026
st.lpi.var.deg{1} = [1 2 2];                                                % MMP, 11/19/2026
%%%gam_guess; % min in [.32 .38]
