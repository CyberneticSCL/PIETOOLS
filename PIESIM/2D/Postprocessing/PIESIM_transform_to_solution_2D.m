%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PIESIM_transform_to_solution.m     PIETOOLS 2024
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function solution=PIESIM_transform_to_solution_2D(psize, PIE, Dop, uinput, grid, solcoeff, opts);
% This routine transforms solution from the Chebyshev coefficient space of the fundamental states to
% the physical solution of the primary states in 2D
% It first performs a transform of Chebshev coefficients of fundamental
% states to Chebyshev coefficients of primary states using Mcheb_nonsaquare
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
% it is array of size (N+1) x (nx+ny), nx - number of states dependending on only x, 
%                                      ny - number of states dependending on only y 
% --- solution.final.pde{2} - array containing the solution for states that are the functions of two variables - 
% it is array of size (N+1) x (N+1) x ns, ns - number of states dependending on both x and y
% --- solution.final.ode - array of size no - ode solution at a final time 
% --- solution.final.observed - array of size noo  - final value of observed outputs
% --- solution.final.regulated - array of size nro  - final value of regulated outputs

% IF OPTS.INTSCHEME=1 (BDF) OPTION, there are additional outputs
% --- solution.timedep.dtime - array of size 1 x Nsteps - array of temporal stamps (discrete time values) of the time-dependent solution
% --- solution.timedep.pde{1,2} - array containing all time-dependent pde solutions 
% --- solution.timedep.pde{1} - array containing the solution for states that are only the functions of one variable - 
% it is array of size (N+1) x (nx+ny) x Nsteps, nx - number of states dependending on only x, 
%                                      ny - number of states dependending on only y 
% --- solution.timedep.pde{2} - array containing the solution for states that are the functions of two variables - 
% it is array of size (N+1) x (N+1) x ns x Nsteps, ns - number of states dependending on both x and y
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
     Mcheb_nonsquare=Dop.Mcheb_nonsquare; 
     if isfield(Dop,'Mcheb0_nonsquare') 
         Mcheb_nonsquare=Dop.Mcheb0_nonsquare;
     end
     C1op=Dop.C1cheb;
     C2op=Dop.C2cheb;
     D11op=PIE.D11.R00;
     D12op=PIE.D12.R00;
     D21op=PIE.D21.R00;
     D22op=PIE.D22.R00;
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

 %  acheb_p=inv(Dop.Mcheb_inv)*solcoeff.final;
          
 % Reconstuctution of primary states in the physical space    
     
     % Define inhomogeneous contributions due to boundary disturbances 
     if (psize.nw>0)
     wvec(:,1)=double(subs(uinput.w(:),tf));

     if ~isempty(Tw.Rx0) 
     bc{1}=Tw.Rx0*wvec;
     else
     bc{1}=zeros(size(Tw,1));
     end

     if ~isempty(Tw.Ry0) 
     bc{2}=Tw.Ry0*wvec;
     else
     bc{2}=zeros(size(Tw,1));
     end

     if ~isempty(Tw.R20) 
     bc{3}=Tw.R20*wvec;
     else
     bc{3}=zeros(size(Tw,1));
     end
     
     nshift=[0 nx nx+ny];
     for k=1:3
     if (isa(bc{k},'polynomial') & bc{k}.maxdeg>0)   
     coeff=zeros(bc{k}.maxdeg+1,size(bc{k},1));
     if(~isempty(bc{k}.degmat))
         for nterms=1:bc{k}.nterms 
         coeff(bc{k}.degmat(nterms)+1,:)=bc{k}.coefficient(nterms,:);
         end
     else
         coeff(:)=bc{k}.coefficient(:);
     end
     coeff=flipud(coeff);
     for n=1:size(coeff,2);
     bcw_input(:,n+nshift(k))=polyval(coeff(:,n),grid.comp);
     end
     else
     for n=1:length(bc{k})
         if (k<3)
     bcw_input{n+nshift(k)}(1:length(grid.comp),:)=bc{k}(n);
         else
     bcw_input{n+nshift(k)}(1:length(grid.comp),1:length(grid.comp))=bc{k}(n);
         end
     end
     end % isa(bc{k},'polynomial')
     end %k=1:3
     end % psize.nw>0
     
     
     % Define inhomogeneous contributions due to boundary inputs
     if (psize.nu>0)
     uvec(:,1)=double(subs(uinput.u(:),tf));
     bc=Tu.Q2*uvec;

     if (isa(bc,'polynomial') & bc{k}.maxdeg>0)  
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

     
     % Reconstruct solution using inverse Chebyshev transform 

     % Reconstruct 1D states

     for n=1:nx
         n1=no+1+(n-1)*(N+1);
         n2=no+n*(N+1);
     acheb_p_local=acheb_p(n1:n2);
     solution.final.pde{1}(:,n) = ifcht(acheb_p_local);
      % Add contribution to solution due to inhomogeneous boundary inputs
      % and disturbances 
     if(psize.nw>0 && ~isempty(bcw_input) && (~isdouble(bcw_input) || ~all(all(bcw_input==0))))
         solution.final.pde{1}(:,n)=solution.final.pde{1}(:,n)+bcw_input{n};
     end % endif
     if(psize.nu>0 && ~isempty(bcu_input) && (~isdouble(bcu_input) || ~all(all(bcu_input==0))))
     solution.final.pde{1}(:,n)=solution.final.pde{1}(:,n)+bcu_input(:,n);
     end % endif
     end % ns

     for n=1:ny
         ng=nx+n;
         n1=no+nx*(N+1)+1+(n-1)*(N+1);
         n2=no+nx*(N+1)+n*(N+1);
     acheb_p_local=acheb_p(n1:n2);
     solution.final.pde{1}(:,ng) = ifcht(acheb_p_local);
      % Add contribution to solution due to inhomogeneous boundary inputs
      % and disturbances 
     if(psize.nw>0 && ~isempty(bcw_input) && (~isdouble(bcw_input) || ~all(all(bcw_input==0))))
         solution.final.pde{1}(:,ng)=solution.final.pde{1}(:,ng)+bcw_input{ng};
     end % endif
     if(psize.nu>0 && ~isempty(bcu_input) && (~isdouble(bcu_input) || ~all(all(bcu_input==0))))
     solution.final.pde{1}(:,ng)=solution.final.pde{1}(:,ng)+bcu_input(:,ng);
     end % endif
     end % ns

     % Reconstruct 2D states
     for n=1:ns
         ng=nx+ny+n;
     n1=no+(nx+ny)*(N+1)+1+(n-1)*(N+1)^2;
     n2=no+(nx+ny)*(N+1)+n*(N+1)^2;
     acheb_p_local_2D=acheb_p(n1:n2);
     acheb_local_2D=reshape(acheb_p_local_2D,N+1,N+1);
     solution.final.pde{2}(:,:,n)=fcgltran2d(acheb_local_2D,0);
      % Add contribution to solution due to inhomogeneous boundary inputs
      % and disturbances

      if(psize.nw>0 && ~isempty(bcw_input) && (~isdouble(bcw_input) || ~all(all(bcw_input==0))))
      solution.final.pde{2}(:,:,n)=solution.final.pde{2}(:,:,n)+bcw_input{ng};
      end
      if(psize.nu>0 && ~isempty(bcu_input) && (~isdouble(bcu_input) || ~all(all(bcu_input==0))))
      solution.final.pde{2}(:,:,n)=solution.final.pde{2}(:,:,n)+bcu_input(:,ng);
      end
     end
     solution.final.ode=solcoeff.final(1:no);

%   Reconstruction of observed and regulated outputs

    if(psize.nro>0)
    solution.final.regulated=C1op*solcoeff.final;
    if (psize.nw>0)
    solution.final.regulated=solution.final.regulated+D11op*wvec(:,1);
    end
    if (psize.nu>0)
    solution.final.regulated=solution.final.regulated+D12op*uvec(:,1);
    end
    end

    if(psize.noo>0)
    solution.final.observed=C2op*solcoeff.final;
    if (psize.nw>0)
    solution.final.observed=solution.final.observed+D21op*wvec(:,1);
    end
    if (psize.nu>0)
    solution.final.observed=solution.final.observed+D22op*uvec(:,1);
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
         acheb_p = Mcheb_nonsquare*solcoeff.timedep.coeff(:,ntime);

   %      acheb_p=inv(Dop.Mcheb_inv)*solcoeff.timedep.coeff(:,ntime);
         
         tt=solution.timedep.dtime(ntime);
         
         
 % Reconstuctution of primary states in the physical space  
     
     % Define inhomogeneous contributions due to disturbances
     if (psize.nw>0)
     wvec(:,1)=double(subs(uinput.w(:),tt));

      if ~isempty(Tw.Rx0) 
     bc{1}=Tw.Rx0*wvec;
     else
     bc{1}=zeros(size(Tw,1));
     end

     if ~isempty(Tw.Ry0) 
     bc{2}=Tw.Ry0*wvec;
     else
     bc{2}=zeros(size(Tw,1));
     end

     if ~isempty(Tw.R20) 
     bc{3}=Tw.R20*wvec;
     else
     bc{3}=zeros(size(Tw,1));
     end
     
     nshift=[0 nx nx+ny];
     for k=1:3

     if (isa(bc{k},'polynomial') & bc{k}.maxdeg>0)     
     coeff=zeros(bc{k}.maxdeg+1,size(bc{k},1));
     if(~isempty(bc{k}.degmat))
         for nterms=1:bc{k}.nterms 
         coeff(bc{k}.degmat(nterms)+1,:)=bc{k}.coefficient(nterms,:);
         end
     else
         coeff(:)=bc{k}.coefficient(:);
     end
     coeff=flipud(coeff);
     for n=1:size(coeff,2);
     bcw_input(:,n+nshift(k))=polyval(coeff(:,n),grid.comp);
     end
     else
    for n=1:length(bc{k})
         if (k<3)
     bcw_input{n+nshift(k)}(1:length(grid.comp),:)=bc{k}(n);
         else
     bcw_input{n+nshift(k)}(1:length(grid.comp),1:length(grid.comp))=bc{k}(n);
         end
    end % n=1:length(bc{k})
     end % isa(bc{k},'polynomial')

     end %k=1:3
     end % psize.nw>0

     % Add disturbances to regulated and observed outputs
     if(psize.nro>0)
     solution.timedep.regulated(:,ntime)=solution.timedep.regulated(:,ntime)+D11op*wvec(:,1);
     end
     if(psize.noo>0)
     solution.timedep.observed(:,ntime)=solution.timedep.observed(:,ntime)+D21op*wvec(:,1);
     end

     % Define inhomogeneous contributions due to controlled inputs
     if (psize.nu>0)
     uvec(:,1)=double(subs(uinput.u(:),tt));
     bc=Tu.Q2*uvec;
     if (isa(bc,'polynomial') & bc{k}.maxdeg>0)  
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
     if(psize.nro>0)
     solution.timedep.regulated(:,ntime)=solution.timedep.regulated(:,ntime)+D12op*uvec(:,1);
     end
     if(psize.noo>0)
     solution.timedep.observed(:,ntime)=solution.timedep.observed(:,ntime)+D22op*uvec(:,1);
     end
     end % if (psize.nu>0)
     

     % Reconstruct solution using inverse Chebyshev transform 

     % Reconstruct 1D states

     for n=1:nx
         n1=no+1+(n-1)*(N+1);
         n2=no+n*(N+1);
     acheb_p_local=acheb_p(n1:n2);
     solution.timedep.pde{1}(:,n,ntime) = ifcht(acheb_p_local);
      % Add contribution to solution due to inhomogeneous boundary inputs
      % and disturbances 
     if(psize.nw>0 && ~isempty(bcw_input) && (~isdouble(bcw_input) || ~all(all(bcw_input==0))))
         solution.timedep.pde{1}(:,n,ntime)=solution.timedep.pde{1}(:,n,ntime)+bcw_input{n};
     end % endif
     if(psize.nu>0 && ~isempty(bcu_input) && (~isdouble(bcu_input) || ~all(all(bcu_input==0))))
     solution.timedep.pde{1}(:,n,ntime)=solution.timedep.pde{1}(:,n,ntime)+bcu_input(:,n);
     end % endif
     end % ns

     for n=1:ny
         ng=nx+n;
         n1=no+nx*(N+1)+1+(n-1)*(N+1);
         n2=no+nx*(N+1)+n*(N+1);
     acheb_p_local=acheb_p(n1:n2);
     solution.timedep.pde{1}(:,ng,ntime) = ifcht(acheb_p_local);
      % Add contribution to solution due to inhomogeneous boundary inputs
      % and disturbances 
     if(psize.nw>0 && ~isempty(bcw_input) && (~isdouble(bcw_input) || ~all(all(bcw_input==0)))) 
         solution.timedep.pde{1}(:,ng,ntime)=solution.timedep.pde{1}(:,ng,ntime)+bcw_input{ng};
     end % endif
     if(psize.nu>0 && ~isempty(bcu_input) && (~isdouble(bcu_input) || ~all(all(bcu_input==0))))
     solution.timedep.pde{1}(:,ng,ntime)=solution.timedep.pde{1}(:,ng,ntime)+bcu_input(:,ng);
     end % endif
     end % ns

     % Reconstruct 2D states
     for n=1:ns
     ng=nx+ny+n;
         n1=no+(nx+ny)*(N+1)+1+(n-1)*(N+1)^2;
         n2=no+(nx+ny)*(N+1)+n*(N+1)^2;
     acheb_p_local=acheb_p(n1:n2);
     acheb_local_2D=reshape(acheb_p_local,N+1,N+1);
     solution.timedep.pde{2}(:,:,n,ntime)=fcgltran2d(acheb_local_2D,0);

%     acheb_p_local=acheb_p(no+1+(n-1)*(N+1):no+n*(N+1));

      % Add contribution to solution due to inhomogeneous boundary inputs
      % and disturbances 
     if(psize.nw>0 && ~isempty(bcw_input) && (~isdouble(bcw_input) || ~all(all(bcw_input==0))))
         solution.timedep.pde{2}(:,:,n,ntime)=solution.timedep.pde{2}(:,:,n,ntime)+bcw_input{ng};
     end
     if(psize.nu>0 && ~isempty(bcu_input) && (~isdouble(bcu_input) || ~all(all(bcu_input==0))))
     solution.timedep.pde{2}(:,:,n,ntime)=solution.timedep.pde{2}(:,:,n,ntime)+bcu_input(:,ng);
     end
     end

    
 end
     
     
end
%  
