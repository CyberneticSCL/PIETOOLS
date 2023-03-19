%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PIESIM_transform_to_solution.m     PIETOOLS 2021b
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function solution=PIESIM_transform_to_solution(psize, PIE, Dop, uinput, grid, solcoeff, opts);
% This routine transforms solution from the Chebyshev coefficient space of the fundamental states to
% the physical solution of the primary states
% It first performs a transform of Chebshev coefficients of fundamental
% states to Chebyshev coefficients of primary states using Mcheb_nonsaquare
% operator
% It then transforms Chebysehv coefficients of the primary states to the
% physical solution of the preimary states

% Inputs: 
% 1) psize - size of the PIE problem: nw, nu, nf, nx 
% 2) PIE structure - needed for reconstruction of non-zero boundary terms,
%  observed and regulated outputs
% 3) Dop- discretized PIE operators
% 4) uinput -  user-defined boundary inputs, forcing functions and initial conditions
% 5) grid - physical and computational grid for n0 states
% 6) solcoeff - Chebyshev coefficients of the time-dependent PDE solution
% 7) opts - options for simulation parameters
% 
% Outputs: 
% 1) solution 
% solution is a strucutre with the following fields
% --- solution.tf - scalar - actual final time of the solution
% --- solution.final.pde - array of size (N+1) x ns, ns=n0+n1+n2 - pde (distributed state) solution at a final time
% --- solution.final.ode - array of size nx - ode solution at a final time 
% --- solution.final.observed - array of size ny  - final value of observed outputs
% --- solution.final.regulated - array of size nz  - final value of regulated outputs

% IF OPTS.INTSCHEME=1 (BDF) OPTION, there are additional outputs
% --- solution.timedep.dtime - array of size 1 x Nsteps - array of temporal stamps (discrete time values) of the time-dependent solution
% --- solution.timedep.pde - array of size (N+1) x ns x Nsteps - time-dependent
% --- solution.timedep.ode - array of size nx x Nsteps - time-dependent solution of nx ODE states
%     solution of ns PDE (distributed) states of the primary PDE system
% --- solution.timedep.observed - array of size ny x Nsteps -
%     time-dependent value of observed outputs
% --- solution.timedep.regulated - array of size nz x Nsteps -
%     time-dependent value of regulated outputs
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Requires: ifcht (inverse Chebyshev transform)
%
% If you modify this code, document all changes carefully and include date
% authorship, and a brief description of modifications
%
% Initial coding YP  - 9_19_2021
% YP 6/16/2022. Added reconstruction of observed and regulated outputs.
% Added a support of case when LHS PDE operator is not the same as the
% state map operator. Note: input arguments changed.

%----------------------------------------   
% We first transform the final solution, since it is available for all
% integration schemes
%----------------------------------------

% Define local variables
     N=psize.N;
     nx=psize.nx;
     ns=sum(psize.n);
     Mcheb_nonsquare=Dop.Mcheb_nonsquare; 
     if isfield(Dop,'Mcheb0_nonsquare') 
         Mcheb_nonsquare=Dop.Mcheb0_nonsquare;
     end
     C1op=Dop.C1cheb;
     C2op=Dop.C2cheb;
     D11op=PIE.D11.P;
     D12op=PIE.D12.P;
     D21op=PIE.D21.P;
     D22op=PIE.D22.P;
     Tu=PIE.Tu;
     Tw=PIE.Tw;
   
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
     if(~isempty(bc.degmat))
         for nterms=1:bc.nterms 
         coeff(bc.degmat(nterms)+1,:)=bc.coefficient(nterms,:);
         end
     else
         coeff(:)=bc.coefficient(:);
     end
     coeff=flipud(coeff);
     for n=1:size(coeff,2);
     bcw_input(:,n)=polyval(coeff(:,n),grid.comp);
     end
     else
     for n=1:length(bc)
     bcw_input(1:length(grid.comp),n)=bc(n);
     end
     end % isa(bc,'polynomial')
     end % psize.nw>0
     
     
     % Define inhomogeneous contributions due to boundary inputs
     if (psize.nu>0)
     uvec(:,1)=double(subs(uinput.u(:),tf));
     bc=Tu.Q2*uvec;

     if isa(bc,'polynomial') 
     coeff=zeros(bc.maxdeg+1,size(bc,1));
     if(~isempty(bc.degmat))
     for nterms=1:bc.nterms 
       coeff(bc.degmat(nterms)+1,:)=bc.coefficient(nterms,:);
     end
     else
     coeff(:)=bc.coefficient(:);
     end

     coeff=flipud(coeff);
     for n=1:size(coeff,2);
     bcu_input(:,n)=polyval(coeff(:,n),grid.comp);
     end
     else
     for n=1:length(bc)
     bcu_input(1:length(grid.comp),n)=bc(n);
     end
     end
     end %psize.nu>0
     
     % Reconstruct solution using inverse Chebyshev transform (ifcht
     % function)
     if(~isnan(acheb_p))
     for n=1:ns
     acheb_p_local=acheb_p(nx+1+(n-1)*(N+1):nx+n*(N+1));
     solution.final.pde(:,n) = ifcht(acheb_p_local);
      % Add contribution to solution due to inhomogeneous boundary inputs
      % and disturbances

      if(psize.nw>0 & ~isempty(bcw_input))
      solution.final.pde(:,n)=solution.final.pde(:,n)+bcw_input(:,n);
      end
      if(psize.nu>0 & ~isempty(bcu_input))
      solution.final.pde(:,n)=solution.final.pde(:,n)+bcu_input(:,n);
      end
     end
     solution.final.ode=solcoeff.final(1:nx);
     else
         solution.final.ode=[];
         if (ns>0) 
             solution.final.pde=[];
         end
     end


%   Reconstruction of observed and regulated outputs

    if(psize.nz>0)
        % Don't execute if solution blew up
        if(~isnan(solcoeff.final)) 
    solution.final.regulated=C1op*solcoeff.final;
    if (psize.nw>0)
    solution.final.regulated=solution.final.regulated+D11op*wvec(:,1);
    end
    if (psize.nu>0)
    solution.final.regulated=solution.final.regulated+D12op*uvec(:,1);
    end
        else
        solution.final.regulated=[];
        end % if
    end % psize

    if(psize.ny>0)
        % Don't execute if solution blew up
        if(~isnan(solcoeff.final)) 
    solution.final.observed=C2op*solcoeff.final;
    if (psize.nw>0)
    solution.final.observed=solution.final.observed+D21op*wvec(:,1);
    end
    if (psize.nu>0)
    solution.final.observed=solution.final.observed+D22op*uvec(:,1);
    end
    else
        solution.final.observed=[];
        end % if
    end % psize
     
     
%----------------------------------------   
% We now transform time-dependent solution, which is available for BDF
% scheme only (opts.IntScheme=1)
%----------------------------------------
          
if (opts.intScheme==1)
    
    % Define ODE solution and temporal stamps array
     solution.timedep.dtime=solcoeff.timedep.dtime;
     if (~isnan(solcoeff.timedep.coeff))
     solution.timedep.ode=solcoeff.timedep.coeff(1:nx,:);


%   Reconstruction of observed and regulated outputs
    if(psize.nz>0)
    solution.timedep.regulated=C1op*solcoeff.timedep.coeff;
    end
    if(psize.ny>0)
    solution.timedep.observed=C2op*solcoeff.timedep.coeff;
    end
    else
    solution.timedep.ode=[];
    solution.timedep.regulated=[];
    solution.timedep.observed=[];
    end
          
     
     
 %---------------------------------------------    
 % Reconstruct PDE solution for every time step
 %--------------------------------------------
     for ntime=1:size(solution.timedep.dtime,2);
         
 % Reconstruction of the coefficients of the primary states
         acheb_p = Mcheb_nonsquare*solcoeff.timedep.coeff(:,ntime);
         
         tt=solution.timedep.dtime(ntime);
         
         
 % Reconstuctution of primary states in the physical space  
     
     % Define inhomogeneous contributions due to disturbances
     if (psize.nw>0)
     wvec(:,1)=double(subs(uinput.w(:),tt));
     bc=Tw.Q2*wvec;

     if isa(bc,'polynomial') 
     coeff=zeros(bc.maxdeg+1,size(bc,1));
     if (~isempty(bc.degmat))
     for nterms=1:bc.nterms 
         coeff(bc.degmat(nterms)+1,:)=bc.coefficient(nterms,:);
     end
     else
     coeff(:)=bc.coefficient(:);
     end % isempty(bc.degmat)
     coeff=flipud(coeff);
     for n=1:size(coeff,2);
     bcw_input(:,n)=polyval(coeff(:,n),grid.comp);
     end % for loop
     else 
     for n=1:length(bc)
     bcw_input(1:length(grid.comp),n)=bc(n);
     end
     end % if isa(bc,'polynomial')

     % Add disturbances to regulated and observed outputs
     if(psize.nz>0 & ~isempty(solution.timedep.regulated))
     solution.timedep.regulated(:,ntime)=solution.timedep.regulated(:,ntime)+D11op*wvec(:,1);
     end
     if(psize.ny>0 & ~isempty(solution.timedep.observed))
     solution.timedep.observed(:,ntime)=solution.timedep.observed(:,ntime)+D21op*wvec(:,1);
     end

     end % if (psize.nw>0)

     % Define inhomogeneous contributions due to controlled inputs
     if (psize.nu>0)
     uvec(:,1)=double(subs(uinput.u(:),tt));
     bc=Tu.Q2*uvec;
     if isa(bc,'polynomial') 
     coeff=zeros(bc.maxdeg+1,size(bc,1));
     if (~isempty(bc.degmat))
     for nterms=1:bc.nterms 
       coeff(bc.degmat(nterms)+1,:)=bc.coefficient(nterms,:);
     end
     else
     coeff(:)=bc.coefficient(:);
     end

     coeff=flipud(coeff);
     for n=1:size(coeff,2);
     bcu_input(:,n)=polyval(coeff(:,n),grid.comp);
     end
     else
     for n=1:length(bc)
     bcu_input(1:length(grid.comp),n)=bc(n);
     end
     end
     % Add controled inputs to regulated and observed outputs
     if(psize.nz>0 & ~isempty(solution.timedep.regulated))
     solution.timedep.regulated(:,ntime)=solution.timedep.regulated(:,ntime)+D12op*uvec(:,1);
     end
     if(psize.ny>0& ~isempty(solution.timedep.observed))
     solution.timedep.observed(:,ntime)=solution.timedep.observed(:,ntime)+D22op*uvec(:,1);
     end
     end % if (psize.nu>0)
     
     % Reconstruct solution using inverse Chebyshev transform (ifcht function)
     if (~isnan(acheb_p))
     for n=1:ns
     acheb_p_local=acheb_p(nx+1+(n-1)*(N+1):nx+n*(N+1));
     solution.timedep.pde(:,n,ntime) = ifcht(acheb_p_local);
      % Add contribution to solution due to inhomogeneous boundary inputs
      % and disturbances 
     if(psize.nw>0 & ~isempty(bcw_input))
         solution.timedep.pde(:,n,ntime)=solution.timedep.pde(:,n,ntime)+bcw_input(:,n);
     end % endif
     if(psize.nu>0 & ~isempty(bcu_input))
     solution.timedep.pde(:,n,ntime)=solution.timedep.pde(:,n,ntime)+bcu_input(:,n);
     end % endif
     end % ns
     else
     solution.timepde.pde=[];
     end
     end
end
%  
