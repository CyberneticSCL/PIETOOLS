%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PIESIM_int_bdf.m     PIETOOLS 2021b
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% % Implicit time-advancement. 

% BDF (backward-difference formulation) is used
% Unconditionally stable for diffusive
% problems, not occasionally stable for hyperbolic problems
% Input:
% Atotal
% dt
% Norder
% Nsteps
% acheb_f0
% nf
% ns
% Nsize
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% If you modify this code, document all changes carefully and include date
% authorship, and a brief description of modifications
%
% Initial coding YP  - 5_26_2021


% Define BDF coefficients and matrices
%

if (Norder>=1) 
bdf1=[1 -1];
A_impl_1=eye(Nsize)*bdf1(1)-dt*Atotal;
A_inv_impl_1=inv(A_impl_1);
end
if (Norder>=2)
bdf2=[3/2 -2 1/2];
A_impl_2=eye(Nsize)*bdf2(1)-dt*Atotal;
A_inv_impl_2=inv(A_impl_2); 
end
if (Norder>=3)
bdf3=[11/6 -3 3/2 -1/3];
A_impl_3=eye(Nsize)*bdf3(1)-dt*Atotal;
A_inv_impl_3=inv(A_impl_3);
end
if (Norder==4)
bdf4=[25/12 -4 3 -4/3 1/4];
A_impl_4=eye(Nsize)*bdf4(1)-dt*Atotal;
A_inv_impl_4=inv(A_impl_4);
end

% achebi denotes Chebyshev coefficients on the fundamental states computed with
% (implicit) BDF method

achebi=acheb_f0;
achebi_p=acheb_f0;
achebi_pp=acheb_f0;
achebi_ppp=acheb_f0;
achebi_pppp=acheb_f0;

% Implicit solution with BDF
t=0;


   for n=1:Nsteps
     t=n*dt;
     inhom=0;
     
     if (nw>0)
     if (isempty(B1cheb)==0&any(B1cheb,'all')~=0)
     wvec(:,1)=double(subs(uinput.w(:),t));
     inhom=inhom+Mcheb_inv*B1cheb*wvec;
     end
     if (isempty(Twcheb)==0&any(Twcheb,'all')~=0)
     wdotvec(:,1)=double(subs(uinput.wdot(:),t));
     inhom=inhom-Mcheb_inv*Twcheb*wdotvec;
     end
     end
     
     if (nu>0)
     if (isempty(B2cheb)==0&any(B2cheb,'all')~=0)
     uvec(:,1)=double(subs(uinput.u(:),t));
     inhom=inhom+Mcheb_inv*B2cheb*uvec;
     end
     if (isempty(Tucheb)==0&any(Tucheb,'all')~=0)
     udotvec(:,1)=double(subs(uinput.udot(:),t));
     inhom=inhom-Mcheb_inv*Tucheb*udotvec;
     end
     end
     
%      bcinf_1=B1op_1*bcvec-Twop_1*bcvecder;
%      bcinf_s=B1op_s*bcvec-Twop_s*bcvecder;
          
     
%     index=1;
%     for i=1:ns;
%     bcinput=bcinput+Mcheb_inv(:,index)*bcinf_1(i)+Mcheb_inv(:,index+1)*bcinf_s(i);
%     index=index+N-p(i)+1;
%     end
%     end 
%     
    force_contr=0;
    for i=1:nf
        force_contr=force_contr+Mcheb_inv*double(subs(tforce(i),t))*acheb_force(:,i);
    end
    
     inhom=inhom+force_contr;
    
    % Time advancement for first several time steps needs to be done with
    % lower order schemes since not enough previous solutions are yet
    % available
    
    
      achebi_pppp=achebi_ppp;
      achebi_ppp=achebi_pp;
      achebi_pp=achebi_p;
      achebi_p=achebi;
      
    
    switch n
        case 1 
        achebi =  A_inv_impl_1*(-bdf1(2)*achebi_p+inhom*dt);
    
        case 2
        if (Norder==1) 
        achebi =  A_inv_impl_1*(-bdf1(2)*achebi_p+inhom*dt);
        elseif (Norder>1) 
        achebi = A_inv_impl_2*(-bdf2(2)*achebi_p-bdf2(3)*achebi_pp+...
        inhom*dt);
        end
        
        case 3
        if (Norder==1) 
        achebi =  A_inv_impl_1*(-bdf1(2)*achebi_p+inhom*dt);
        elseif (Norder==2) 
        achebi = A_inv_impl_2*(-bdf2(2)*achebi_p-bdf2(3)*achebi_pp+...
        inhom*dt);
        elseif (Norder>2) 
        achebi =  A_inv_impl_3*(-bdf3(2)*achebi_p-bdf3(3)*achebi_pp-bdf3(4)*achebi_ppp+...
        inhom*dt);
        end
        
        otherwise
            
        switch Norder
        case 1    
        achebi =  A_inv_impl_1*(-bdf1(2)*achebi_p+inhom*dt);
        case 2
        achebi =  A_inv_impl_2*(-bdf2(2)*achebi_p-bdf2(3)*achebi_pp+...
        inhom*dt);    
        case 3
        achebi =  A_inv_impl_3*(-bdf3(2)*achebi_p-bdf3(3)*achebi_pp-bdf3(4)*achebi_ppp+...
       inhom*dt);
        case 4
        achebi =  A_inv_impl_4*(-bdf4(2)*achebi_p-bdf4(3)*achebi_pp-bdf4(4)*achebi_ppp-...
        bdf4(5)*achebi_pppp+inhom*dt);
        end  
        
    end
            
    % Save temporal evolution of ODE states for plotting
%     if (nx>0)
    odesol(n,1:nx)=achebi(1:nx);
    dtime(n)=t;
%     end
    pdesol(n,:) = achebi(nx+1:end);
    end % Loop on time steps
        
   % Chebyshev coefficients to output
       acheb_f=achebi;
       tf=t;