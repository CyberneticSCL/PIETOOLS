<<<<<<< Updated upstream:PIETOOLS/PIESIM/PIESIM_routines/Time_Integrators/PIESIM_int_symbolic.m
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PIESIM_int_symbolic.m    PIETOOLS 2024
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Performs time-advancement of a discretized PIE using a convolution integral in
% a symbolic form
% Convolution integral is an exact solution in time of a
% discretized in space PIE (similar to an analytical solution of
% system of ODEs)
% Provides the most accurate solution but slow for large N (due to the use
% of symbolic integration)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Inputs:
% psize - size of the PIE problem: nw, nu, no 
% opts - options for temporal scheme parameters
% uinput   - user-defined boundary inputs
% coeff - Chebyshev coefficients of initial conditions and forcing terms, if any
% Dop - discretized PIE operators
%
% % Output:
% solcoeff - contains:
% 1) solcoeff.acheb_f - Chebyshev coefficients of the final
% solution, ODE+PDE states
% 2) solcoeff.tf - final time of the solution
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% If you modify this code, document all changes carefully and include date
% authorship, and a brief description of modifications
%
% Initial coding YP  - 5_26_2021

function solcoeff=PIESIM_int_symbolic(psize, opts, uinput, coeff, Dop)
        syms tt intsym st;

        nw=psize.nw;
        nu=psize.nu;
        no=psize.no;        
        tf=opts.tf;
        
        Twcheb=Dop.Twcheb;
        Tucheb=Dop.Tucheb;
        B1cheb=Dop.B1cheb;
        B2cheb=Dop.B2cheb;
        Acheb=Dop.Acheb;
        Tcheb_inv=Dop.Tcheb_inv;
        Atotal=Dop.Atotal;
        V=Dop.V;
        D=Dop.D;
        Nsize=size(Atotal,1);
        acheb_f0=coeff.acheb_f0;
        
           inhom=zeros(Nsize,Nsize);
           for j=1:Nsize
             lambda=D(j,j);
             if (nw>0)
                 
    % Contribution to inhomogeneous term due to boundary disturbances  
             if (isempty(B1cheb)==0&any(B1cheb,'all')~=0)
             for i=1:size(uinput.w,2);      
             intsym=int(exp(lambda*(tt-st))*uinput.w(i),st,[0,tt]);
             wvec(i,1)=double(subs(intsym,tt,tf));
             end        
             end 
             
    % Contribution to inhomogeneous term due to time derivative of boundary disturbances             
             if (isempty(Twcheb)==0&any(Twcheb,'all')~=0)
             for i=1:size(uinput.wdot,2);      
             intsym=int(exp(lambda*(tt-st))*uinput.wdot(i),st,[0,tt]);
             wdotvec(i,1)=double(subs(intsym,tt,tf));
             end 
             end
            
             
             if (isempty(B1cheb)==0&any(B1cheb,'all')~=0)
             inhom(j,:)=inhom(j,:)+(B1cheb*coeff.w*wvec).';
             end
             if (isempty(Twcheb)==0&any(Twcheb,'all')~=0)
             inhom(j,:)=inhom(j,:)-(Twcheb*coeff.w*wdotvec).';
             end
             end % nw>0
             
              if (nu>0)
                 
    % Contribution to inhomogeneous term due to boundary inputs  
             if (isempty(B2cheb)==0&any(B2cheb,'all')~=0)
             for i=1:size(uinput.u,2);      
             intsym=int(exp(lambda*(tt-st))*uinput.u(i),st,[0,tt]);
             uvec(i,1)=double(subs(intsym,tt,tf));
             end        
             end 
             
    % Contribution to inhomogeneous term due to time derivative of boundary inputs             
             if (isempty(Tucheb)==0&any(Tucheb,'all')~=0)
             for i=1:size(uinput.udot,2);      
             intsym=int(exp(lambda*(tt-st))*uinput.udot(i),st,[0,tt]);
             udotvec(i,1)=double(subs(intsym,tt,tf));
             end 
             end
            
             
             if (isempty(B2cheb)==0&any(B2cheb,'all')~=0)
             inhom(j,:)=inhom(j,:)+(B2cheb*coeff.u*uvec).';
             end
             if (isempty(Tucheb)==0&any(Tucheb,'all')~=0)
             inhom(j,:)=inhom(j,:)-(Tucheb*coeff.u*udotvec).';
             end
             end % nu>0
  
           end    % for j loop
       
   % Adding all contributions to the time integral. Only considering the
   % columns of inhom array that are non-zero
           
       addition_symbolic=0;
       if (isempty(inhom)==0) 
           [row,col]=find(inhom);
           C=unique(col);
       for k=1:size(C,1)
        i=C(k);
       addition_symbolic=addition_symbolic+diag(inhom(:,i))*inv(V)*Tcheb_inv(:,i);
       end
       end

     
        exp_lambda=exp(diag(D*tf));
        solcoeff.final=V*diag(exp_lambda)/V*acheb_f0+V*addition_symbolic;
        solcoeff.tf=tf;
=======
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PIESIM_int_symbolic.m    PIETOOLS 2024
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Performs time-advancement of a discretized PIE using a convolution integral in
% a symbolic form
% Convolution integral is an exact solution in time of a
% discretized in space PIE (similar to an analytical solution of
% system of ODEs)
% Provides the most accurate solution but slow for large N (due to the use
% of symbolic integration)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Inputs:
% psize - size of the PIE problem: nw, nu, no 
% opts - options for temporal scheme parameters
% uinput   - user-defined boundary inputs
% coeff - Chebyshev coefficients of initial conditions and forcing terms, if any
% Dop - discretized PIE operators
%
% % Output:
% solcoeff - contains:
% 1) solcoeff.acheb_f - Chebyshev coefficients of the final
% solution, ODE+PDE states
% 2) solcoeff.tf - final time of the solution
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% If you modify this code, document all changes carefully and include date
% authorship, and a brief description of modifications
%
% Initial coding YP  - 5_26_2021

function solcoeff=PIESIM_int_symbolic(psize, opts, uinput, coeff, Dop)
        syms tt intsym st;

        nw=psize.nw;
        nu=psize.nu;
        no=psize.no;        
        tf=opts.tf;
        
        Twcheb=Dop.Twcheb;
        Tucheb=Dop.Tucheb;
        B1cheb=Dop.B1cheb;
        B2cheb=Dop.B2cheb;
        Acheb=Dop.Acheb;
        Tcheb_inv=Dop.Tcheb_inv;
        Atotal=Dop.Atotal;
        V=Dop.V;
        D=Dop.D;
        Nsize=size(Atotal,1);
        acheb_f0=coeff.acheb_f0;
        
           inhom=zeros(Nsize,Nsize);
           for j=1:Nsize
             lambda=D(j,j);
             if (nw>0)
                 
    % Contribution to inhomogeneous term due to boundary disturbances  
             if (isempty(B1cheb)==0&any(B1cheb,'all')~=0)
             for i=1:size(uinput.w,2);      
             intsym=int(exp(lambda*(tt-st))*uinput.w(i),st,[0,tt]);
             wvec(i,1)=double(subs(intsym,tt,tf));
             end        
             end 
             
    % Contribution to inhomogeneous term due to time derivative of boundary disturbances             
             if (isempty(Twcheb)==0&any(Twcheb,'all')~=0)
             for i=1:size(uinput.wdot,2);      
             intsym=int(exp(lambda*(tt-st))*uinput.wdot(i),st,[0,tt]);
             wdotvec(i,1)=double(subs(intsym,tt,tf));
             end 
             end
            
             
             if (isempty(B1cheb)==0&any(B1cheb,'all')~=0)
             inhom(j,:)=inhom(j,:)+(B1cheb*coeff.w*wvec).';
             end
             if (isempty(Twcheb)==0&any(Twcheb,'all')~=0)
             inhom(j,:)=inhom(j,:)-(Twcheb*coeff.w*wdotvec).';
             end
             end % nw>0
             
              if (nu>0)
                 
    % Contribution to inhomogeneous term due to boundary inputs  
             if (isempty(B2cheb)==0&any(B2cheb,'all')~=0)
             for i=1:size(uinput.u,2);      
             intsym=int(exp(lambda*(tt-st))*uinput.u(i),st,[0,tt]);
             uvec(i,1)=double(subs(intsym,tt,tf));
             end        
             end 
             
    % Contribution to inhomogeneous term due to time derivative of boundary inputs             
             if (isempty(Tucheb)==0&any(Tucheb,'all')~=0)
             for i=1:size(uinput.udot,2);      
             intsym=int(exp(lambda*(tt-st))*uinput.udot(i),st,[0,tt]);
             udotvec(i,1)=double(subs(intsym,tt,tf));
             end 
             end
            
             
             if (isempty(B2cheb)==0&any(B2cheb,'all')~=0)
             inhom(j,:)=inhom(j,:)+(B2cheb*coeff.u*uvec).';
             end
             if (isempty(Tucheb)==0&any(Tucheb,'all')~=0)
             inhom(j,:)=inhom(j,:)-(Tucheb*coeff.u*udotvec).';
             end
             end % nu>0
  
           end    % for j loop
       
   % Adding all contributions to the time integral. Only considering the
   % columns of inhom array that are non-zero
           
       addition_symbolic=0;
       if (isempty(inhom)==0) 
           [row,col]=find(inhom);
           C=unique(col);
       for k=1:size(C,1)
        i=C(k);
       addition_symbolic=addition_symbolic+diag(inhom(:,i))*inv(V)*Tcheb_inv(:,i);
       end
       end

     
        exp_lambda=exp(diag(D*tf));
        solcoeff.final=V*diag(exp_lambda)/V*acheb_f0+V*addition_symbolic;
        solcoeff.tf=tf;
>>>>>>> Stashed changes:PIESIM/PIESIM_routines/Time_Integrators/PIESIM_int_symbolic.m
        