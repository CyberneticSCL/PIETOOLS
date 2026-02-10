%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PIESIM_int_bdf.m    PIETOOLS 2025
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Performs time-advancement of a discretized PIE with implicit
% backward-difference formula (BDF) - unconditionally stable for diffusive
% problems, occasionally not stable for hyperbolic problems
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
% 2) solcoeff.w and solcoeff.u - Chebyshev coefficients of a spatial
% component of forcing functions introduced through uinput.w and uinput.u
% (time-independent, needed for reconstruction of the solution) 
% 3) solcoeff.tf - final time of the solution 
% 4) solcoeff.timedep.ode -
% time-dependent Chebyshev coefficients of the ODE states of PIE system for
% output (for ODE states, Chebyshev coefficients are equal to the solution
% of the states, since the ODE states are not spatially dependent)
% 5) solcoeff.timedep.pde -
% time-dependent Chebyshev coefficients of the PDE states of PIE system for output 
% (for PDE states, an inverse Fourier transform needs to be performed to recover physical solution from its )
% Chebyshev coefficients) 
% 6)solcoeff.timedep.dtime - temporal stamps (discrete time values) of the time-dependent solution
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% If you modify this code, document all changes carefully and include date
% authorship, and a brief description of modifications
%
% Initial coding YP  - 5_26_2021

% YP 8_18_2024 - added support for spatially-variant multiple disturbances
% across multiple states (thrugh coeff.u, coeff.w)
% Added solcoeff,u, solcoeff.w to solcoeff structure
% YP 12/31/2025 - modified the treatment of disturbances and control inputs via functon
% handles

function solcoeff=PIESIM_int_bdf(psize, opts, uinput, coeff, Dop)
nw=psize.nw;
nu=psize.nu;
no=psize.no;
dt=opts.dt;
Nsteps=opts.Nsteps;
Norder=opts.Norder;

Twcheb=Dop.Twcheb;
Tucheb=Dop.Tucheb;
B1cheb=Dop.B1cheb;
B2cheb=Dop.B2cheb;
Acheb=Dop.Acheb;
Tcheb_inv=Dop.Tcheb_inv;
Atotal=Dop.Atotal;


Nsize=size(Atotal,1);
acheb_f0=coeff.acheb_f0;

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
     t=(n)*dt;
     inhom=0;
     
     if (nw>0)
    % Contribution to inhomogeneous term due to boundary disturbances 
     if (isempty(B1cheb)==0&any(B1cheb,'all')~=0)
         for k = 1:numel(uinput.w)
        wvec(k,1)=uinput.w{k}(t);
        end
     inhom=inhom+Tcheb_inv*B1cheb*coeff.w*wvec;
     end
    % Contribution to inhomogeneous term due to time derivative of boundary disturbances  
     if (isempty(Twcheb)==0&any(Twcheb,'all')~=0)
          for k = 1:numel(uinput.w)
          wdotvec(k,1)=uinput.wdot{k}(t);
          end
     inhom=inhom-Tcheb_inv*Twcheb*coeff.w*wdotvec;
     end
     end
     
     if (nu>0)
    % Contribution to inhomogeneous term due to boundary inputs  
     if (isempty(B2cheb)==0&any(B2cheb,'all')~=0) 
          for k = 1:numel(uinput.u)
          uvec(k,1)=uinput.u{k}(t);
          end
     inhom=inhom+Tcheb_inv*B2cheb*coeff.u*uvec;
     end
    % Contribution to inhomogeneous term due to time derivative of boundary inputs  
     if (isempty(Tucheb)==0&any(Tucheb,'all')~=0)
          for k = 1:numel(uinput.u)
          udotvec(k,1)=uinput.udot{k}(t);
          end
     inhom=inhom-Tcheb_inv*Tucheb*coeff.u*udotvec;
     end
     end
    
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
    if (no>0)
    solcoeff.timedep.ode(n,:)=achebi(1:no);
    end
    
    % Save temporal evolution of PDE states for plotting
    
    solcoeff.timedep.pde(n,:) = achebi(no+1:end);
    
    solcoeff.timedep.dtime(n)=t;
    
    solcoeff.timedep.coeff(:,n)=achebi;
    
    end % Loop on time steps
        
   % Chebyshev coefficients to output
       solcoeff.final=achebi;
       solcoeff.tf=t;

       solcoeff.w=coeff.w;
       solcoeff.u=coeff.u;

       
       