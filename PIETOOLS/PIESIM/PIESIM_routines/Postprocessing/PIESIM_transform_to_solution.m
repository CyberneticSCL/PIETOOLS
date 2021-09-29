%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PIESIM_transform_to_solution.m     PIETOOLS 2021d
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function solution=PIESIM_transform_to_solution(psize, Tu, Tw, Mcheb_nonsquare, uinput, grid, solcoeff, opts);
% This routine transforms solution from the Chebyshev coefficient space of the fundamental states to
% the physical solution of the primary states
% It first performs a transform of Chebshev coefficients of fundamental
% states to Chebyshev coefficients of primary states using Mcheb_nonsaquare
% operator
% It then transforms Chebysehv coefficients of the primary states to the
% physical solution of the preimary states
% Inputs: 
% psize - size of the PIE problem: nw, nu, nf, nx 
% Tu
% Tw
% Mcheb_nonsquare - nonsquare matrix operator that transforms Chebyshev
% cpefficients between fundamental and primary states
% uinput -  user-defined boundary inputs
% grid
% solcoeff - 
% opts - options for simulation parameters
% 
% Outputs: solution. This field contains
% AVAILABLE FOR ALL OPTS.INTSCHEME OPTIONS
% solution.tf - actual final time of the solution
% solution.final.pde - pde (distributed state) solution at a final time : matrix of
% dimension (N+1) x ns, ns=n0+n1+n2
% solution.final.ode - ode solution at a final time : array of dimensions
% nx x 1

% AVAILABLE ONLY FOR OPTS.INTSCHEME=1 (BDF) OPTION
% solution.timedep.ode - array of size nx x Nsteps - time-dependent solution of nx ODE states
% solution.timedep.pde - array of size (N+1) x ns x Nsteps - time-dependent
% solution of ns PDE (distribited) states of the primary PDE system
% solution.timedep.dtime - array of size 1 x Nsteps - array of temporal stamps (discrete time values) of the time-dependent solution
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Requires: ifcht (inverse Chebyshev transform)
%
% If you modify this code, document all changes carefully and include date
% authorship, and a brief description of modifications
%
% Initial coding YP  - 9_19_2021

%----------------------------------------   
% We first transform the final solution, since it is available for all
% integration schemes
%----------------------------------------

     N=psize.N;
     nx=psize.nx;
     ns=sum(psize.n);
   
     if (solcoeff.tf~=opts.tf) 
     disp('Warning: solution final time does not match user input final time');
     disp('Defaulting to solution final time');
     end
     tf=solcoeff.tf;
     solution.tf=tf;
%      

 % Reconstruction of the coefficients of the primary states
     acheb_p = Mcheb_nonsquare*solcoeff.final;
     
     
 % Reconstuctution of primary states in the physical space    
     
     % Define inhomogeneous contributions due to boundary disturbances
     if (psize.nw>0)
     wvec(:,1)=double(subs(uinput.w(:),tf));
     bc=Tw.Q2*wvec;
     
     if isa(bc,'polynomial')         
     coeff=zeros(bc.maxdeg+1,size(bc,1));
     for nterms=1:bc.nterms 
         coeff(bc.degmat(nterms)+1,:)=bc.coefficient(nterms,:);
     end
     coeff=flipud(coeff);
     for n=1:size(coeff,2);
     bcw_input(:,n)=polyval(coeff(:,n),grid.comp);
     end
     else
     bcw_input=0;
     end
     end % psize.nw>0
     
     
     % Define inhomogeneous contributions due to boundary inputs
     if (psize.nu>0)
     uvec(:,1)=double(subs(uinput.u(:),tf));
     bc=Tu.Q2*uvec;
     if isa(bc,'polynomial') 
     coeff=zeros(bc.maxdeg+1,size(bc,1));
     for nterms=1:bc.nterms 
       coeff(bc.degmat(nterms)+1,:)=bc.coefficient(nterms,:);
     end
     coeff=flipud(coeff);
     for n=1:size(coeff,2);
     bcu_input(:,n)=polyval(coeff(:,n),grid.comp);
     end
     else
     bcu_input=0;
     end
     end
     
     % Reconstruct solution using inverse Chebyshev transform (ifcht
     % function)
     for n=1:ns
     acheb_p_local=acheb_p(nx+1+(n-1)*(N+1):nx+n*(N+1));
     solution.final.pde(:,n) = ifcht(acheb_p_local);
      % Add contribution to solution due to inhomogeneous boundary inputs
      % and disturbances
     if(psize.nw>0 & bcw_input~=0)
         solution.final.pde(:,n)=solution.final.pde(:,n)+bcw_input(:,n);
     end
     if(psize.nu>0 & bcu_input~=0)
     solution.final.pde(:,n)=solution.final.pde(:,n)+bcu_input(:,n);
     end
     end
     solution.final.ode=solcoeff.final(1:nx);
    
     
     
%----------------------------------------   
% We now transform time-dependent solution, which is available for BDF
% scheme only (opts.IntScheme=1)
%----------------------------------------
          
if (opts.intScheme==1)
    
    % Define ODE solution and temporal stamps array
     solution.timedep.dtime=solcoeff.timedep.dtime;
     solution.timedep.ode=solcoeff.timedep.coeff(1:nx,:);
     
     
 %---------------------------------------------    
 % Reconstruct PDE solution for every time step
 %--------------------------------------------
     for ntime=1:size(solution.timedep.dtime,2);
         
 % Reconstruction of the coefficients of the primary states
         acheb_p = Mcheb_nonsquare*solcoeff.timedep.coeff(:,ntime);
         
         tt=solution.timedep.dtime(ntime);
         
         
 % Reconstuctution of primary states in the physical space  
     
     % Define inhomogeneous contributions due to boundary disturbances
     if (psize.nw>0)
     wvec(:,1)=double(subs(uinput.w(:),tt));
     bc=Tw.Q2*wvec;
     if isa(bc,'polynomial') 
     coeff=zeros(bc.maxdeg+1,size(bc,1));
     for nterms=1:bc.nterms 
         coeff(bc.degmat(nterms)+1,:)=bc.coefficient(nterms,:);
     end
     coeff=flipud(coeff);
     for n=1:size(coeff,2);
     bcw_input(:,n)=polyval(coeff(:,n),grid.comp);
     end
     else
     bcw_input=0;
     end
     end
     
     % Define inhomogeneous contributions due to boundary inputs
     if (psize.nu>0)
     uvec(:,1)=double(subs(uinput.u(:),tt));
     bc=Tu.Q2*uvec;
     if isa(bc,'polynomial') 
     coeff=zeros(bc.maxdeg+1,size(bc,1));
     for nterms=1:bc.nterms 
       coeff(bc.degmat(nterms)+1,:)=bc.coefficient(nterms,:);
     end
     coeff=flipud(coeff);
     for n=1:size(coeff,2);
     bcu_input(:,n)=polyval(coeff(:,n),grid.comp);
     end
     else
     bcu_input=0;
     end
     end
     
     % Reconstruct solution using inverse Chebyshev transform (ifcht
     % function)
     for n=1:ns
     acheb_p_local=acheb_p(nx+1+(n-1)*(N+1):nx+n*(N+1));
     solution.timedep.pde(:,n,ntime) = ifcht(acheb_p_local);
      % Add contribution to solution due to inhomogeneous boundary inputs
      % and disturbances
%     if(psize.nw>0 & bcw_input~=0)        
     if(psize.nw>0 & bcw_input~=0)
         solution.timedep.pde(:,n,ntime)=solution.timedep.pde(:,n,ntime)+bcw_input(:,n);
     end
     if(psize.nu>0 & bcu_input~=0)
     solution.timedep.pde(:,n,ntime)=solution.timedep.pde(:,n,ntime)+bcu_input(:,n);
     end
     end
     end
     
     
end
%  
