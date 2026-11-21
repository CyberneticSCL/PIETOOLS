%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PIESIM_transform_to_solution.m     PIETOOLS 2025
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function solution=PIESIM_transform_to_solution(psize, PIE, Dop, uinput, solcoeff, opts);
% This routine transforms solution from the Chebyshev coefficient space of the fundamental states to
% the physical solution of the primary states
% It first performs a transform of Chebshev coefficients of fundamental
% states to Chebyshev coefficients of primary states using Tcheb_2PDEstate
% operator
% It then transforms Chebysehv coefficients of the primary states to the
% physical solution of the preimary states
%   coefficients of the primary states: Tmap*a + Twmap*w_c*w(t) +
%   Tumap*u_c*u(t); z = C1*a + D11*w_c*w(t) + D12*u_c*u(t); y = C2*a +
%   D21*w_c*w(t) + D22*u_c*u(t), with a the coefficients of the fundamental
%   states and w_c, u_c the spatial coefficients of the inputs.
%
% Inputs: 
% 1) psize - size of the PIE problem: N, n0, n, nw, nu, nr0, nrx, no0, nox
% 2) PIE structure - not used
% 3) Dop- discretized PIE operators
%    ('struct' with fields Tcheb_2PDEstate, Tucheb_2PDEstate,
%    Twcheb_2PDEstate, optionally Tchebmap_2PDEstate, Tuchebmap_2PDEstate,
%    Twchebmap_2PDEstate (used instead when present), C1cheb, C2cheb,
%    D11cheb, D12cheb, D21cheb, D22cheb)
% 4) uinput -  user-defined boundary inputs, forcing functions and initial conditions
%    ('struct' with fields w, u, cells of function handles of t)
% 6) solcoeff - Chebyshev coefficients of the time-dependent PDE solution
%    ('struct' with fields final, tf, w, u and, for BDF, timedep.dtime,
%    timedep.coeff)
% 7) opts - options for simulation parameters
%    ('struct' with fields tf, intScheme)
% 
% Outputs: 
% 1) solution 
% solution is a strucutre with the following fields
% --- solution.tf - scalar - actual final time of the solution
% --- solution.final.primary{1,2} - cell array containing final value of
% solution states
% --- solution.final.primary{1} - array of size n0 - ode (finite-dimensional) solutions at a final time 
% --- solution.final.primary{2} - array of size (N+1) x nx - pde (distributed state) solutions at a final time
% --- solution.final.observed{1,2} - cell array containing final value of observed outputs 
% --- solution.final.observed[1} - array of size no0 containing final
%      value of finite-dimensional observed outputs
% --- solution.final.observed{2} - array of size (N+1) x nox 
%      containing final value of infinite-dimensional observed outputs
% --- solution.final.regulated{1,2} - cell array containing final value of regulated outputs 
% --- solution.final.regulated[1} - array of size nr0 containing
%     final value of finite-dimensional regulated outputs
% --- solution.final.regulated{2} - array of size (N+1) x nrx 
%      containing final value of infinite-dimensional regulated ooutputs
%
% IF OPTS.INTSCHEME=1 (BDF) OPTION, there are additional outputs
% --- solution.timedep.dtime - array of size 1 x Nsteps - array of temporal stamps (discrete time values) of the time-dependent solution
% --- solution.timedep.primary{1,2} - cell array containing final value of
% solution states
% --- solution.timedep.primary{1} - array of size n0 x Nsteps - time-dependent solutions of n0 ODE (finite-dimensional) states
% --- solution.timedep.primary{2} - array of size (N+1) x nx x Nsteps - time-dependent
%     solution of ns PDE (distributed) states of the primary PDE system
% --- solution.timedep.observed{1,2} - cell array containing time-dependent observed outputs 
% --- solution.timedep.observed[1} - array of size no0 x Nsteps -
%     time-dependent value of finite-dimensional observed outputs
% --- solution.timedep.observed{2} - array of size (N+1) x nox x Nsteps
%      containing infinite-dimensional observed outputs
% --- solution.timedep.regulated{1,2} - cell array containing time-dependent regulated outputs 
% --- solution.timedep.regulated[1} - array of size nr0 x Nsteps -
%     time-dependent value of finite-dimensional regulated outputs
% --- solution.timedep.regulated{2} - array of size (N+1) x nrx x Nsteps
% containing infinite-dimensional regulated ooutputs

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PIETOOLS - PIESIM_transform_to_solution
%
% Copyright (C) 2026 PIETOOLS Team
%
% This program is free software; you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation; either version 2 of the License, or
% (at your option) any later version.
%
% This program is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
%
% You should have received a copy of the GNU General Public License
% along with this program; if not, write to the Free Software
% Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% If you modify this code, document all changes carefully and include date
% authorship, and a brief description of modifications
%
% For support, contact M. Peet, Arizona State University at mpeet@asu.edu,
% or Y. Peet at ypeet@asu.edu
% Initial coding YP  - 9_19_2021
% YP, 6/16/2022 - Added reconstruction of observed and regulated outputs.
% Added a support of case when LHS PDE operator is not the same as the
% state map operator. note: input arguments changed.
% VJ, YP, MMP, 11/21/2026: Changes since PIETOOLS 2025: input grid removed;
%   solution fields ode, pde replaced by primary{1,2}, and observed,
%   regulated are cells {1,2} (finite-dimensional outputs; outputs of s in
%   physical space); the input terms enter through Dop.Tucheb_2PDEstate,
%   Twcheb_2PDEstate (or their map versions) with uinput.w, uinput.u cells
%   of function handles of t; Dop.Tcheb_2PDEstate (or Tchebmap_2PDEstate)
%   replaces Mcheb_nonsquare; psize fields n0, nr0, no0, nrx, nox; for a
%   diverged (NaN) solution timedep.primary{2} is empty (was timepde.pde).

%----------------------------------------   
% We first transform the final solution, since it is available for all
% integration schemes
%----------------------------------------

% Define local variables

     N=psize.N;
     n0=psize.n0;
     nx=sum(psize.n);

% Determine number of observed and regulated outputs

    if ~isfield(psize,'nrx')
        psize.nrx=0;
    end 
    if ~isfield(psize,'nox')
       psize.nox=0;
    end

     % Maps to the primary-state coefficients; Tmap, Tumap, Twmap replace
     % T, Tu, Tw when present
     Tcheb_2PDEstate=Dop.Tcheb_2PDEstate; 
     Tucheb_2PDEstate=Dop.Tucheb_2PDEstate;
     Twcheb_2PDEstate=Dop.Twcheb_2PDEstate;

      if isfield(Dop,'Tchebmap_2PDEstate') 
          Tcheb_2PDEstate=Dop.Tchebmap_2PDEstate;
      end
      if isfield(Dop,'Tuchebmap_2PDEstate') 
          Tucheb_2PDEstate=Dop.Tuchebmap_2PDEstate;
      end
       if isfield(Dop,'Twchebmap_2PDEstate') 
          Twcheb_2PDEstate=Dop.Twchebmap_2PDEstate;
      end

     % Output operators
     C1op=Dop.C1cheb;
     C2op=Dop.C2cheb;
     D11op=Dop.D11cheb;
     D12op=Dop.D12cheb;
     D21op=Dop.D21cheb;
     D22op=Dop.D22cheb;

     % Final time of the computed solution
     if (solcoeff.tf~=opts.tf) 
     disp('Warning: solution final time does n0t match user input final time');
     disp('Defaulting to solution final time');
     end
     tf=solcoeff.tf;
     solution.tf=tf;
%      

 % Reconstruction of the coefficients of the primary states
     acheb_p = Tcheb_2PDEstate*solcoeff.final;
          
 % Reconstuctution of primary states in the physical space    
     
     % Define inhomogeneous contributions due to boundary disturbances
     if (psize.nw>0)
        for k = 1:numel(uinput.w)
        wvec(k,1)=uinput.w{k}(tf); 
        end
        acheb_p=acheb_p+Twcheb_2PDEstate*solcoeff.w*wvec;
     end % psize.nw>0
     
     
     % Define inhomogeneous contributions due to boundary inputs
     if (psize.nu>0)
        for k = 1:numel(uinput.u)
        uvec(k,1)=uinput.u{k}(tf); 
        end
     acheb_p=acheb_p+Tucheb_2PDEstate*solcoeff.u*uvec;
     end %psize.nu>0
     
     % Reconstruct solution using inverse Chebyshev transform (ifcht
     % function)
     if(~isnan(acheb_p))
     for n=1:nx
     acheb_p_local=acheb_p(n0+1+(n-1)*(N+1):n0+n*(N+1));
     solution.final.primary{2}(:,n) = ifcht(acheb_p_local);
     end
     if (n0>0)
     solution.final.primary{1}=solcoeff.final(1:n0);
     else
     solution.final.primary{1}=[]; 
     end
     else
         solution.final.primary{1}=[];
%        if (ns>0)                                                          % VJ, YP, MMP, 11/21/2026 (was)
         if (nx>0)  % nx, the number of PDE states (ns unset)               % VJ, YP, MMP, 11/21/2026
             solution.final.primary{2}=[];
         end
     end

  %   Reconstruction of observed and regulated outputs

 solution.final.regulated=[];
 solution.final.observed=[];

        % Don't execute if solution blew up
        if(~isnan(solcoeff.final)) 

if (psize.nr0+psize.nrx>0)
%   Reconstruction of regulated outputs

%    Chebyshev coefficients of state-to-output portion
        coeff_final_regulated=C1op*solcoeff.final;

% Add contribution to coefficients from disturbances and control inputs
        if (psize.nw>0)
         coeff_final_regulated=coeff_final_regulated+D11op*solcoeff.w*wvec(:,1);
        end
         if (psize.nu>0)
        coeff_final_regulated=coeff_final_regulated+D12op*solcoeff.u*uvec(:,1);
         end

% Transform to physical space

% Finite-dimensional outputs
   if(psize.nr0>0)
    solution.final.regulated{1}=coeff_final_regulated(1:psize.nr0);
   end 
% Infinite-dimensional outputs
    for n=1:psize.nrx
     coeff_final_regulated_local=coeff_final_regulated(psize.nr0+1+(n-1)*(N+1):psize.nr0+n*(N+1));
     solution.final.regulated{2}(:,n) = ifcht(coeff_final_regulated_local);
    end

end % if (psize.nr+psize.nrx>0)

     %   Reconstruction of observed outputs

     if (psize.no0+psize.nox>0)

%    Chebyshev coefficients of state-to-output portion
        coeff_final_observed=C2op*solcoeff.final;

% Add contribution to coefficients from disturbances and control inputs
        if (psize.nw>0)
         coeff_final_observed=coeff_final_observed+D21op*solcoeff.w*wvec(:,1);
        end
         if (psize.nu>0)
        coeff_final_observed=coeff_final_observed+D22op*solcoeff.u*uvec(:,1);
         end

% Transform to physical space

% Finite-dimensional outputs
   if(psize.no0>0)
    solution.final.observed{1}=coeff_final_observed(1:psize.no0);
   end 
% Infinite-dimensional outputs
    for n=1:psize.nox
     coeff_final_observed_local=coeff_final_observed(psize.no0+1+(n-1)*(N+1):psize.no0+n*(N+1));
     solution.final.observed{2}(:,n) = ifcht(coeff_final_observed_local);
    end

end % if (psize.no0+psize.nox>0)
     
        end % if (~isnaan)
     
%----------------------------------------   
% We now transform time-dependent solution, which is available for BDF
% scheme only (opts.IntScheme=1) and at final time not equal to 0
%----------------------------------------
          
if (opts.intScheme==1 & solution.tf~=0)
    
    % Define ODE solution and temporal stamps array
     solution.timedep.dtime=solcoeff.timedep.dtime;
     if (~isnan(solcoeff.timedep.coeff) & n0>0)
     solution.timedep.primary{1}=solcoeff.timedep.coeff(1:n0,:);
    else
    solution.timedep.primary{1}=[];
     end % ~isnan(solcoeff.timedep.coeff)
         
 %---------------------------------------------    
 % Reconstruct PDE and output solutions for every time step
 %--------------------------------------------

 solution.timedep.regulated=[];
 solution.timedep.observed=[];

     for ntime=1:size(solution.timedep.dtime,2);

 % Reconstruction of the coefficients of the primary states
         acheb_p = Tcheb_2PDEstate*solcoeff.timedep.coeff(:,ntime);
         
         tt=solution.timedep.dtime(ntime);
         
         
 % Reconstuctution of primary states in the physical space  
     
     % Define inhomogeneous contributions due to disturbances
     if (psize.nw>0)

       for k = 1:numel(uinput.w)
         wvec(k,1)= uinput.w{k}(tt);
       end
     acheb_p=acheb_p+Twcheb_2PDEstate*solcoeff.w*wvec;
     end % if (psize.nw>0)

     % Define inhomogeneous contributions due to controlled inputs
     if (psize.nu>0)

     for k = 1:numel(uinput.u)
     uvec(k,1)=uinput.u{k}(tt); 
     end
     acheb_p=acheb_p+Tucheb_2PDEstate*solcoeff.u*uvec;
     end % if (psize.nu>0)
  
     % Reconstruct solution using inverse Chebyshev transform (ifcht function)
     if (~isnan(acheb_p))
     for n=1:nx
     acheb_p_local=acheb_p(n0+1+(n-1)*(N+1):n0+n*(N+1));
     solution.timedep.primary{2}(:,n,ntime) = ifcht(acheb_p_local);
      % Add contribution to solution due to inhomogeneous boundary inputs
      % and disturbances 
     end % for n=1:ns
     else
%    solution.timepde.pde=[];                                               % VJ, YP, MMP, 11/21/2026 (was)
     solution.timedep.primary{2}=[];  % as final.primary{2} above           % VJ, YP, MMP, 11/21/2026
     end % ~isnan(acheb_p)

%   Reconstruction of regulated outputs

    if (psize.nr0+psize.nrx>0)

%    Chebyshev coefficients of state-to-output portion
        coeff_timedep_regulated=C1op*solcoeff.timedep.coeff(:,ntime);

% Add contribution to coefficients from disturbances and control inputs
        if (psize.nw>0)
         coeff_timedep_regulated=coeff_timedep_regulated+D11op*solcoeff.w*wvec(:,1);
        end
         if (psize.nu>0)
        coeff_timedep_regulated=coeff_timedep_regulated+D12op*solcoeff.u*uvec(:,1);
         end

% Transform to physical space

% Finite-dimensional outputs
   if(psize.nr0>0)
    solution.timedep.regulated{1}(:,ntime)=coeff_timedep_regulated(1:psize.nr0);
   end 
% Infinite-dimensional outputs
    for n=1:psize.nrx
     coeff_timedep_regulated_local=coeff_timedep_regulated(psize.nr0+1+(n-1)*(N+1):psize.nr0+n*(N+1));
     solution.timedep.regulated{2}(:,n,ntime) = ifcht(coeff_timedep_regulated_local);
    end

    end %   if (psize.nr+psize.nrx>0)

     %   Reconstruction of observed outputs

     if (psize.no0+psize.nox>0)

%    Chebyshev coefficients of state-to-output portion
        coeff_timedep_observed=C2op*solcoeff.timedep.coeff(:,ntime);

% Add contribution to coefficients from disturbances and control inputs
        if (psize.nw>0)
         coeff_timedep_observed=coeff_timedep_observed+D21op*solcoeff.w*wvec(:,1);
        end
         if (psize.nu>0)
        coeff_timedep_observed=coeff_timedep_observed+D22op*solcoeff.u*uvec(:,1);
         end

% Transform to physical space

% Finite-dimensional outputs
   if(psize.no0>0)
    solution.timedep.observed{1}(:,ntime)=coeff_timedep_observed(1:psize.no0);
   end 
% Infinite-dimensional outputs
    for n=1:psize.nox
     coeff_timedep_observed_local=coeff_timedep_observed(psize.no0+1+(n-1)*(N+1):psize.no0+n*(N+1));
     solution.timedep.observed{2}(:,n,ntime) = ifcht(coeff_timedep_observed_local);
    end

     end % if (psize.no0+psize.nox>0)
     
     end % for ntime

end % if (options)
%  
