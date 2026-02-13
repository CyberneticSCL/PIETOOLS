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
% grid - physical grids of solution states
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
% YP 2/12/2026 - added handling of non-separable disturbances

function solcoeff=PIESIM_int_bdf(psize, opts, uinput, coeff, grid, Dop)
syms st;
nw=psize.nw;
if (psize.dim==2)
N=opts.N;
Np=[1, N(1)+1, N(2)+1, prod(N+1)];

nw_groups=[psize.nw0, psize.nwx, psize.nwy, psize.nw2];
nw_points=nw_groups.*Np;
cum_nwg = cumsum(nw_groups);
cum_nwp = cumsum(nw_points);

nu_groups=[psize.nu0, psize.nux, psize.nuy, psize.nu2];
nu_points=nu_groups.*Np;
cum_nug = cumsum(nu_groups);
cum_nup = cumsum(nu_points);
end

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
             if (uinput.wsep{k})
        wvec(k,1)=uinput.w{k}(t);
             else
        uinput.wspace{k}=subs(uinput.w{k},st,t);
        % Relate the disturbance index to its type
        group = find(k <= cum_nwg, 1);
        index=cum_nwp(group-1) +(k-cum_nwg(group-1)-1)*Np(group)+1;  
             switch group
             case 2
        coeff.w(index:index+N(1),k)=PIESIM_NonPoly2Mat_cheb(N(1), uinput.wspace{k}, 0, grid.x);
             case 3
        coeff.w(index:index+N(2),k)=PIESIM_NonPoly2Mat_cheb(N(2), uinput.wspace{k}, 0, grid.y);
             case 4
        coeff.w(index:index+prod(N+1)-1,k)=PIESIM_NonPoly2Mat_cheb_2D(N, uinput.wspace{k}, 0, grid);
             end
        wvec(k,1)=1;
             end % if uinput.wsep{k}
         end % k 
     inhom=inhom+Tcheb_inv*B1cheb*coeff.w*wvec;
     end

    % Contribution to inhomogeneous term due to time derivative of boundary disturbances  
     if (isempty(Twcheb)==0&any(Twcheb,'all')~=0)
          for k = 1:numel(uinput.w)
                if (uinput.wsep{k})
          wdotvec(k,1)=uinput.wdot{k}(t);
                else
        uinput.wspace{k}=subs(uinput.wdot{k},st,t);
        % Relate the disturbance index to its type
        group = find(k <= cum_nwg, 1);
        index=cum_nwp(group-1) +(k-cum_nwg(group-1)-1)*Np(group)+1;  
             switch group
             case 2
        coeff.w(index:index+N(1),k)=PIESIM_NonPoly2Mat_cheb(N(1), uinput.wspace{k}, 0, grid.x);
             case 3
        coeff.w(index:index+N(2),k)=PIESIM_NonPoly2Mat_cheb(N(2), uinput.wspace{k}, 0, grid.y);
             case 4
        coeff.w(index:index+prod(N+1)-1,k)=PIESIM_NonPoly2Mat_cheb_2D(N, uinput.wspace{k}, 0, grid);
             end
          wdotvec(k,1)=1;
                end % (uinput.wsep{k})
          end % k 
     inhom=inhom-Tcheb_inv*Twcheb*coeff.w*wdotvec;
     end
     end
     
     if (nu>0)
    % Contribution to inhomogeneous term due to boundary inputs  
     if (isempty(B2cheb)==0&any(B2cheb,'all')~=0) 
          for k = 1:numel(uinput.u)
                if (uinput.usep{k})
          uvec(k,1)=uinput.u{k}(t);
                else  uinput.uspace{k}=subs(uinput.u{k},st,t);
        % Relate the control input index to its type
        group = find(k <= cum_nug, 1);
        index=cum_nup(group-1) +(k-cum_nug(group-1)-1)*Np(group)+1;  
             switch group
             case 2
        coeff.u(index:index+N(1),k)=PIESIM_NonPoly2Mat_cheb(N(1), uinput.uspace{k}, 0, grid.x);
             case 3
        coeff.u(index:index+N(2),k)=PIESIM_NonPoly2Mat_cheb(N(2), uinput.uspace{k}, 0, grid.y);
             case 4
        coeff.u(index:index+prod(N+1)-1,k)=PIESIM_NonPoly2Mat_cheb_2D(N, uinput.uspace{k}, 0, grid);
             end
        uvec(k,1)=1;
             end % if uinput.usep{k}
         end % k 
     inhom=inhom+Tcheb_inv*B2cheb*coeff.u*uvec;
     end
    % Contribution to inhomogeneous term due to time derivative of boundary inputs  
     if (isempty(Tucheb)==0&any(Tucheb,'all')~=0)
          for k = 1:numel(uinput.u)
                if (uinput.usep{k})
          udotvec(k,1)=uinput.udot{k}(t);
                else
 uinput.uspace{k}=subs(uinput.udot{k},st,t);
        % Relate the disturbance index to its type
        group = find(k <= cum_nug, 1);
        index=cum_nup(group-1) +(k-cum_nug(group-1)-1)*Np(group)+1;  
             switch group
             case 2
        coeff.u(index:index+N(1),k)=PIESIM_NonPoly2Mat_cheb(N(1), uinput.uspace{k}, 0, grid.x);
             case 3
        coeff.u(index:index+N(2),k)=PIESIM_NonPoly2Mat_cheb(N(2), uinput.uspace{k}, 0, grid.y);
             case 4
        coeff.u(index:index+prod(N+1)-1,k)=PIESIM_NonPoly2Mat_cheb_2D(N, uinput.uspace{k}, 0, grid);
             end
          udotvec(k,1)=1;
                end % (uinput.usep{k})
          end % k 
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

       
       