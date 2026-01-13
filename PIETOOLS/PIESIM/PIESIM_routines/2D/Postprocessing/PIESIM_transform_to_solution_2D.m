%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PIESIM_transform_to_solution.m     PIETOOLS 2025
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function solution=PIESIM_transform_to_solution_2D(psize, PIE, Dop, uinput, grid, solcoeff, opts);
% This routine transforms solution from the Chebyshev coefficient space of the fundamental states to
% the physical solution of the primary states in 2D
% It first performs a transform of Chebshev coefficients of fundamental
% states to Chebyshev coefficients of primary states using Tcheb_nonsaquare
% operator
% It then transforms Chebysehv coefficients of the primary states to the
% physical solution of the preimary states

% Inputs: 
% 1) psize - size of the PIE problem: all variables defining the size of the PIE problem
% 2) PIE structure - needed for reconstruction of non-zero boundary terms,
%  observed and regulated outputs
% 3) Dop- discretized PIE operators
% 4) uinput -  user-defined boundary inputs, forcing functions and initial conditions
% 5) grid - physical and computational grid for states differentiable up to order zero (corresponding to a primary = PDE state discretization)
% 6) solcoeff - Chebyshev coefficients of the time-dependent PDE solution
% 7) opts - options for simulation parameters
% 
% Outputs: 
% 1) solution 
% solution is a structure with the following fields
% --- solution.tf - scalar - actual final time of the solution
% --- solution.final.pde{1,2} - array containing all pde solutions at a final time
% --- solution.final.pde{1} - array containing the solution for states that are only the functions of one variable - 
% it is array of size (N+1) x (nx+ny), nx - number of states depending only on s1, 
%                                      ny - number of states depending only on s2 
% --- solution.final.pde{2} - array containing the solution for states that are the functions of two variables - 
% it is array of size (N+1) x (N+1) x ns, ns - number of states depending on both s1 and s2
% --- solution.final.ode - array of size no - ode solution at a final time 
% --- solution.final.observed - array of size noo  - final value of observed outputs
% --- solution.final.regulated - array of size nro  - final value of regulated outputs

% IF OPTS.INTSCHEME=1 (BDF) OPTION, there are additional outputs
% --- solution.timedep.dtime - array of size 1 x Nsteps - array of temporal stamps (discrete time values) of the time-dependent solution
% --- solution.timedep.pde{1,2} - array containing all time-dependent pde solutions 
% --- solution.timedep.pde{1} - array containing the solution for states that are only the functions of one variable - 
% it is array of size (N+1) x (nx+ny) x Nsteps, nx - number of states depending only on s1, 
%                                      ny - number of states depending only on s2 
% --- solution.timedep.pde{2} - array containing the solution for states that are the functions of two variables - 
% it is array of size (N+1) x (N+1) x ns x Nsteps, ns - number of states depending on both s1 and s2
% --- solution.timedep.ode - array of size no x Nsteps - time-dependent solution of no ODE states
% --- solution.timedep.observed - array of size noo x Nsteps -
%     time-dependent value of observed outputs
% --- solution.timedep.regulated - array of size nro x Nsteps -
%     time-dependent value of regulated outputs
%-------------------------------------
%
% Requires: ifcht (inverse Chebyshev transform)
%
% If you modify this code, document all changes carefully and include date
% authorship, and a brief description of modifications
%
% Initial coding YP  - 4_16_2024
% YP 12/31/2025 - modified treatment of disturbances and control inputs

%----------------------------------------   
% We first transform the final solution, since it is available for all
% integration schemes
%----------------------------------------

% Define local variables
     N=psize.N;
     no=psize.no;
     nx=sum(psize.nx);
     ny=sum(psize.ny);
     ns=sum(psize.n);
     Tcheb_nonsquare=Dop.Tcheb_nonsquare; 
     Tucheb_nonsquare=Dop.Tucheb_nonsquare;
     Twcheb_nonsquare=Dop.Twcheb_nonsquare;
     if isfield(Dop,'Tchebmap_nonsquare') 
         Tcheb_nonsquare=Dop.Tchebmap_nonsquare;
     end
     if isfield(Dop,'Tuchebmap_nonsquare') 
         Tucheb_nonsquare=Dop.Tuchebmap_nonsquare;
     end
      if isfield(Dop,'Twchebmap_nonsquare') 
         Twcheb_nonsquare=Dop.Twchebmap_nonsquare;
     end
     C1op=Dop.C1cheb;
     C2op=Dop.C2cheb;
     D11op=Dop.D11cheb;
     D12op=Dop.D12cheb;
     D21op=Dop.D21cheb;
     D22op=Dop.D22cheb;

     if (solcoeff.tf~=opts.tf) 
     disp('Warning: solution final time does not match user input final time');
     disp('Defaulting to solution final time');
     end
     tf=solcoeff.tf;
     solution.tf=tf;
%      

 % Reconstruction of the coefficients of the primary states
   acheb_p = Tcheb_nonsquare*solcoeff.final;
          
 % Reconstuctution of primary states in the physical space    
     
     % Define inhomogeneous contributions due to boundary disturbances 
     if (psize.nw>0)
      for k = 1:numel(uinput.w)
        wvec(k,1)=uinput.w{k}(tf); 
      end
     acheb_p=acheb_p+Twcheb_nonsquare*solcoeff.w*wvec;
     end % psize.nw>0


     % Define inhomogeneous contributions due to boundary inputs
     if (psize.nu>0)
        for k = 1:numel(uinput.u)
        uvec(k,1)=uinput.u{k}(tf); 
        end
     acheb_p=acheb_p+Tucheb_nonsquare*solcoeff.u*uvec;
     end % psize.nu>0
     
     
     % Reconstruct solution using inverse Chebyshev transform 

     % Reconstruct 1D states

     for n=1:nx
         n1=no+1+(n-1)*(N+1);
         n2=no+n*(N+1);
     acheb_p_local=acheb_p(n1:n2);
     solution.final.pde{1}(:,n) = ifcht(acheb_p_local);
     end % nx

     for n=1:ny
         ng=nx+n;
         n1=no+nx*(N+1)+1+(n-1)*(N+1);
         n2=no+nx*(N+1)+n*(N+1);
     acheb_p_local=acheb_p(n1:n2);
     solution.final.pde{1}(:,ng) = ifcht(acheb_p_local);
     end % ny

     % Reconstruct 2D states
     for n=1:ns
         ng=nx+ny+n;
     n1=no+(nx+ny)*(N+1)+1+(n-1)*(N+1)^2;
     n2=no+(nx+ny)*(N+1)+n*(N+1)^2;
     acheb_p_local_2D=acheb_p(n1:n2);
     acheb_local_2D=reshape(acheb_p_local_2D,N+1,N+1);
     solution.final.pde{2}(:,:,n)=fcgltran2d(acheb_local_2D,0);
     end % ns
     solution.final.ode=solcoeff.final(1:no);

%   Reconstruction of observed and regulated outputs

    if(psize.nro>0)
    solution.final.regulated=C1op*solcoeff.final;
    if (psize.nw>0)
    solution.final.regulated=solution.final.regulated+D11op*solcoeff.w*wvec(:,1);
    end
    if (psize.nu>0)
    solution.final.regulated=solution.final.regulated+D12op*solcoeff.u*uvec(:,1);
    end
    end

    if(psize.noo>0)
    solution.final.observed=C2op*solcoeff.final;
    if (psize.nw>0)
    solution.final.observed=solution.final.observed+D21op*solcoeff.w*wvec(:,1);
    end
    if (psize.nu>0)
    solution.final.observed=solution.final.observed+D22op*solcoeff.u*uvec(:,1);
    end
    end
     
     
%----------------------------------------   
% We now transform time-dependent solution, which is available for BDF
% scheme only (opts.IntScheme=1)
%----------------------------------------
          
if (opts.intScheme==1&opts.tf~=0)
    
    % Define ODE solution and temporal stamps array
     solution.timedep.dtime=solcoeff.timedep.dtime;
     solution.timedep.ode=solcoeff.timedep.coeff(1:no,:);


%   Reconstruction of observed and regulated outputs
    if(psize.nro>0)
    solution.timedep.regulated=C1op*solcoeff.timedep.coeff;
    end
    if(psize.noo>0)
    solution.timedep.observed=C2op*solcoeff.timedep.coeff;
    end
     
     
 %---------------------------------------------    
 % Reconstruct PDE solution for every time step
 %--------------------------------------------
     for ntime=1:size(solution.timedep.dtime,2);
         
 % Reconstruction of the coefficients of the primary states
         acheb_p = Tcheb_nonsquare*solcoeff.timedep.coeff(:,ntime);
         
         tt=solution.timedep.dtime(ntime);
         
         
 % Reconstuctution of primary states in the physical space  
     
     % Define inhomogeneous contributions due to disturbances
     if (psize.nw>0)
        for k = 1:numel(uinput.w)
        wvec(k,1)=uinput.w{k}(tt); 
     end
     acheb_p=acheb_p+Twcheb_nonsquare*solcoeff.w*wvec;

     % Add disturbances to regulated and observed outputs
     if(psize.nro>0)
     solution.timedep.regulated(:,ntime)=solution.timedep.regulated(:,ntime)+D11op*solcoeff.w*wvec(:,1);
     end
     if(psize.noo>0)
     solution.timedep.observed(:,ntime)=solution.timedep.observed(:,ntime)+D21op*solcoeff.w*wvec(:,1);
     end

     end % psize.nw>0

     % Define inhomogeneous contributions due to controlled inputs
    if (psize.nu>0)
    for k = 1:numel(uinput.u)
    uvec(k,1)=uinput.u{k}(tt); 
    end
     acheb_p=acheb_p+Tucheb_nonsquare*solcoeff.u*uvec;

      % Add control inputs to regulated and observed outputs
     if(psize.nro>0)
     solution.timedep.regulated(:,ntime)=solution.timedep.regulated(:,ntime)+D12op*solcoeff.u*uvec(:,1);
     end
     if(psize.noo>0)
     solution.timedep.observed(:,ntime)=solution.timedep.observed(:,ntime)+D22op*solcoeff.u*uvec(:,1);
     end

     end % psize.nu>0
     

     % Reconstruct solution using inverse Chebyshev transform 

     % Reconstruct 1D states

     for n=1:nx
         n1=no+1+(n-1)*(N+1);
         n2=no+n*(N+1);
     acheb_p_local=acheb_p(n1:n2);
     solution.timedep.pde{1}(:,n,ntime) = ifcht(acheb_p_local);
     end % nx

     for n=1:ny
         ng=nx+n;
         n1=no+nx*(N+1)+1+(n-1)*(N+1);
         n2=no+nx*(N+1)+n*(N+1);
     acheb_p_local=acheb_p(n1:n2);
     solution.timedep.pde{1}(:,ng,ntime) = ifcht(acheb_p_local);
     end % ny

     % Reconstruct 2D states
     for n=1:ns
     ng=nx+ny+n;
         n1=no+(nx+ny)*(N+1)+1+(n-1)*(N+1)^2;
         n2=no+(nx+ny)*(N+1)+n*(N+1)^2;
     acheb_p_local=acheb_p(n1:n2);
     acheb_local_2D=reshape(acheb_p_local,N+1,N+1);
     solution.timedep.pde{2}(:,:,n,ntime)=fcgltran2d(acheb_local_2D,0);
     end % ns

    
     end % ntime

end %  if (opts.intScheme==1&opts.tf~=0)
     
%  
