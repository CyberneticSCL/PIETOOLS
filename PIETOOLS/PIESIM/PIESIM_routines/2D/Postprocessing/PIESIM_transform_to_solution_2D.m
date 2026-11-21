%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PIESIM_transform_to_solution_2D.m     PIETOOLS 2025
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function solution=PIESIM_transform_to_solution_2D(psize, PIE, Dop, uinput, gridall, solcoeff, opts);
% This routine transforms solution from the Chebyshev coefficient space of the fundamental states to
% the physical solution of the primary states in 2D
% It first performs a transform of Chebshev coefficients of fundamental
% states to Chebyshev coefficients of primary states using Tcheb_2PDEstate
% operator
% It then transforms Chebysehv coefficients of the primary states to the
% physical solution of the preimary states
%   coefficients of the primary states: Tmap*a + Twmap*w_c*w(t) +
%   Tumap*u_c*u(t); z = C1*a + D11*w_c*w(t) + D12*u_c*u(t); y = C2*a +
%   D21*w_c*w(t) + D22*u_c*u(t), with a the coefficients of the fundamental
%   states and w_c, u_c the spatial coefficients of the inputs (for a
%   non-separable input recomputed at each time, with value 1).
%
% Inputs: 
% 1) psize - size of the PIE problem: all variables defining the size of the PIE problem
%    ('struct' with fields N, n0, nx, ny, n, nw0, nwx, nwy, nw2, nu0, nux,
%    nuy, nu2, nr0, nrx, nry, nr2, no0, nox, noy, no2, nw, nu)
% 2) PIE structure - not used
% 3) Dop- discretized PIE operators
%    ('struct' with fields Tcheb_2PDEstate, Tucheb_2PDEstate,
%    Twcheb_2PDEstate, optionally Tchebmap_2PDEstate, Tuchebmap_2PDEstate,
%    Twchebmap_2PDEstate (used instead when present), C1cheb, C2cheb,
%    D11cheb, D12cheb, D21cheb, D22cheb)
% 4) uinput -  user-defined boundary inputs, forcing functions and initial conditions
%    ('struct' with fields w, u (function handles of t if separable, sym
%    in st, sx, sy otherwise), wsep, usep)
% 5) gridall - physical and computational gridall for states differentiable up to order zero (corresponding to a primary = PDE state discretization)
%    ('struct' with cell fields x, y)
% 6) solcoeff - Chebyshev coefficients of the time-dependent PDE solution
%    ('struct' with fields final, tf, w, u and, for BDF, timedep.dtime,
%    timedep.coeff)
% 7) opts - options for simulation parameters
%    ('struct' with fields tf, intScheme)
% 
% Outputs: 
% 1) solution 
% solution is a structure with the following fields
% --- solution.tf - scalar - actual final time of the solution
% --- solution.final.primary{1,2,3,4} - cell array containing all primary state solutions (ode and pde) at a final time
% --- solution.final.primary{1} - array of size n0 - ode (finite-dimensional) solutions at a final time 
% --- solution.final.primary{2} - array containing the solution for states that are only the functions of s1 - 
%      array of size (N(1)+1) x nx, nx - number of states depending only on s1
% --- solution.final.primary{3} - array containing the solution for states that are only the functions of s2 - 
%      array of size (N(2)+1) x ny, ny - number of states depending only on s2
% --- solution.final.primary{4} - array containing the solution for states that are the functions of two variables - 
% it is array of size (N(1)+1) x (N(2)+1) x n2, n2 - number of states depending on both s1 and s2
% --- solution.final.observed{1,2,3,4} - cell array containing final value of observed outputs 
% --- solution.final.observed{1} - array of size no0  - final value of finite-dimensional observed outputs
% --- solution.final.observed{2} - array containing final value of infinite-dimnesional 
%      observed outputs that are only the functions of s1 - 
%      array of size (N(1)+1) x nox, nox - number of outputs depending only on s1
% --- solution.final.observed{3} - array containing final value of infinite-dimnesional 
%      observed outputs that are only the functions of s2 - 
%      array of size (N(2)+1) x noy, noy - number of outputs depending
%      only on s2
% --- solution.final.observed{4} - array containing final value of observed outputs that are the functions of two variables - 
%      array of size (N(1)+1) x (N(2)+1) x no2, no2 - number of outputs depending on both s1 and s2
% --- solution.final.regulated{1,2,3,4} - cell array containing final value of regulatedd outputs 
% --- solution.final.regulated{1} - array of size nr0  - final value of finite-dimensional regulated outputs
% --- solution.final.regulated{2} - array containing final value of infinite-dimensional 
%      regulated outputs that are only the functions of s1 - 
%      array of size (N(1)+1) x nrx, nrx - number of outputs depending only on s1  
%      solution.final.regulated{3} - array containing final value of infinite-dimensional 
%      regulated outputs that are only the functions of s2 - 
%      array of size (N(2)+1) x nry, nry - number of outputs depending only on s2 
% --- solution.final.regulated{4} - array containing final value of regulated outputs that are the functions of two variables - 
% It is array of size (N(1)+1) x (N(2)+1) x nr2, nr2 - number of outputs depending on both s1 and s2
%
% IF OPTS.INTSCHEME=1 (BDF) OPTION, there are additional outputs
% --- solution.timedep.dtime - array of size 1 x Nsteps - array of temporal stamps (discrete time values) of the time-dependent solution
% --- solution.timedep.primary{1,2,3,4} - cell array containing all
% time-dependent primary state solutions (ode and pde)
% --- solution.timedep.primary{1} - array of size n0 x Nsteps - time-dependent solution of no ODE (finite-dimensional) states
% --- solution.timedep.primary{2} - array containing the solution for states that are only the functions of s1 - 
%      array of size (N(1)+1) x nx x Nsteps, nx - number of states depending only on s1
% --- solution.timedep.primary{3} - array containing the solution for states that are only the functions of s2 - 
%      array of size (N(2)+1) x ny x Nsteps, ny - number of states depending only on s2
% --- solution.timedep.primary{4} - array containing the solution for states that are the functions of two variables - 
%      array of size (N(1)+1) x (N(2)+1) x n2 x Nsteps, n2 - number of states depending on both s1 and s2
% --- solution.timedep.observed{1,2,3,4} - cell array containing time-dependent 
%      observed outputs 
% --- solution.timedep.observed[1} - array of size no0 x Nsteps -
%     time-dependent value of finite-dimensional observed outputs
% --- solution.timedep.observed{2} - array containing observed outputs that are only the functions of s1 - 
%     array of size (N(1)+1) x nox x Nsteps, nox - number of observed outputs depending only on s1
% --- solution.timedep.observed{3} - array containing observed outputs that are only the functions of s2 - 
%     array of size (N(2)+1) x noy x Nsteps, noy - number of observed outputs depending only on s2
% --- solution.timedep.observed{4} - array containing observed outputs that are the functions of two variables - 
%      array of size (N(1)+1) x (N(2)+1) x no2 x Nsteps, no2 - number of outputs depending on both s1 and s2
% --- solution.timedep.regulated{1,2,3,4} - array containing time-dependent infinite-dimnesional 
%      regulated outputs 
% --- solution.timedep.regulated{1} - array of size nr0 x Nsteps -
%     time-dependent value of finite-dimensional regulated outputs
% --- solution.timedep.regulated{2} - array containing regulated outputs that are only the functions of s1 - 
%     array of size (N(1)+1) x nrx x Nsteps, nrx - number of regulated outputs depending only on s1
% --- solution.timedep.regulated{3} - array containing regulated outputs that are only the functions of s2 -
%     array of size (N(2)+1) x nry x Nsteps, nry - number of regulated outputs depending only on s2
% --- solution.timedep.regulated{4} - array containing regulated outputs that are the functions of two variables -
%     array of size (N(1)+1) x (N(2)+1) x nr2 x Nsteps, nr2 - number of outptus depending on both s1 and s2

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PIETOOLS - PIESIM_transform_to_solution_2D
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
% Initial coding - YP, 4/16/2024
% YP, MMP, 11/21/2026: Changes since PIETOOLS 2025: input grid renamed
%   gridall; solution fields ode, pde replaced by primary{1,2,3,4} (ODE
%   states; states of s1, of s2, of (s1,s2)), and observed, regulated are
%   cells {1,2,3,4} likewise (outputs of s1, s2, (s1,s2) in physical space);
%   uinput.w, uinput.u are cells, a non-separable input has its spatial
%   coefficients recomputed at each time; Dop fields *_2PDEstate (or
%   *map_2PDEstate) replace *_nonsquare; psize.n0 replaces psize.no; N is
%   [N_x, N_y].

%----------------------------------------   
% We first transform the final solution, since it is available for all
% integration schemes
%----------------------------------------

    syms st;

% Define local variables
    N=psize.N;
    % Inputs grouped by kind (R, L2[x], L2[y], L2[x,y]), Np(g) coefficients
    % per input of kind g; offsets of the groups in solcoeff.w, solcoeff.u
    Np=[1, N(1)+1, N(2)+1, prod(N+1)];

    nw_groups=[psize.nw0, psize.nwx, psize.nwy, psize.nw2];
    nw_points=nw_groups.*Np;
    cum_nwg = cumsum(nw_groups);
    cum_nwp = cumsum(nw_points);

    nu_groups=[psize.nu0, psize.nux, psize.nuy, psize.nu2];
    nu_points=nu_groups.*Np;
    cum_nug = cumsum(nu_groups);
    cum_nup = cumsum(nu_points);

     % Numbers of states and outputs
     n0=psize.n0;
     nx=sum(psize.nx,'all');
     ny=sum(psize.ny,'all');
     n2=sum(psize.n,'all');
     nr_sum=psize.nr0+psize.nrx+psize.nry+psize.nr2;
     no_sum=psize.no0+psize.nox+psize.noy+psize.no2;
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
     disp('Warning: solution final time does not match user input final time');
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
             if (uinput.wsep{k})
        wvec(k,1)=uinput.w{k}(tf); 
             else
        uinput.wspace{k}=subs(uinput.w{k},st,tf);
        % Relate the disturbance index to its type
        group = find(k <= cum_nwg, 1);
        index=cum_nwp(group-1) +(k-cum_nwg(group-1)-1)*Np(group)+1;  
             switch group
             case 2
        solcoeff.w(index:index+N(1),k)=PIESIM_NonPoly2Mat_cheb(N(1), uinput.wspace{k}, 0, gridall.x);
             case 3
        solcoeff.w(index:index+N(2),k)=PIESIM_NonPoly2Mat_cheb(N(2), uinput.wspace{k}, 0, gridall.y);
             case 4
        solcoeff.w(index:index+prod(N+1)-1,k)=PIESIM_NonPoly2Mat_cheb_2D(N, uinput.wspace{k}, 0, gridall);
             end
        wvec(k,1)=1;
             end
      end
     acheb_p=acheb_p+Twcheb_2PDEstate*solcoeff.w*wvec;
     end % psize.nw>0


     % Define inhomogeneous contributions due to boundary inputs
     if (psize.nu>0)
        for k = 1:numel(uinput.u)
             if (uinput.usep{k})
        uvec(k,1)=uinput.u{k}(tf); 
             else uinput.uspace{k}=subs(uinput.u{k},st,tf);
        % Relate the control input index to its type
        group = find(k <= cum_nug, 1);
        index=cum_nup(group-1) +(k-cum_nug(group-1)-1)*Np(group)+1;  
             switch group
             case 2
        solcoeff.u(index:index+N(1),k)=PIESIM_NonPoly2Mat_cheb(N(1), uinput.uspace{k}, 0, gridall.x);
             case 3
        solcoeff.u(index:index+N(2),k)=PIESIM_NonPoly2Mat_cheb(N(2), uinput.uspace{k}, 0, gridall.y);
             case 4
        solcoeff.u(index:index+prod(N+1)-1,k)=PIESIM_NonPoly2Mat_cheb_2D(N, uinput.uspace{k}, 0, gridall);
             end
        uvec(k,1)=1;
             end
      end
     acheb_p=acheb_p+Tucheb_2PDEstate*solcoeff.u*uvec;
     end % psize.nu>0
     
     
     % Reconstruct solution using inverse Chebyshev transform 

     % Reconstruct 0D states

     if (n0>0)
     solution.final.primary{1}=solcoeff.final(1:n0);
     else
     solution.final.primary{1}=[];
     end

     % Reconstruct 1D states

     for n=1:nx
         ni=n0+1+(n-1)*(N(1)+1);
         nii=n0+n*(N(1)+1);
     acheb_p_local=acheb_p(ni:nii);
     solution.final.primary{2}(:,n) = ifcht(acheb_p_local);
     end % nx

     for n=1:ny
         nprev=n0+nx*(N(1)+1);
         ni=nprev+1+(n-1)*(N(2)+1);
         nii=nprev+n*(N(2)+1);
     acheb_p_local=acheb_p(ni:nii);
     solution.final.primary{3}(:,n) = ifcht(acheb_p_local);
     end % ny

     % Reconstruct 2D states
     for n=1:n2
         nprev=n0+(N+1)*[nx;ny];
     ni=nprev+1+(n-1)*prod(N+1);
     nii=nprev+n*prod(N+1);
     acheb_p_local_2D=acheb_p(ni:nii);
     acheb_local_2D=reshape(acheb_p_local_2D,N(1)+1,N(2)+1);
     solution.final.primary{4}(:,:,n)=fcgltran2d(acheb_local_2D,0);
     end % ns

%   Reconstruction of observed and regulated outputs

 solution.final.regulated=[];
 solution.final.observed=[];

% Reconstruction of regulated outputs

    if(nr_sum>0)
    coeff_final_regulated=C1op*solcoeff.final;
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

 % Infinite-dimensional 1D outputs

    for n=1:psize.nrx
     coeff_final_regulated_local=coeff_final_regulated(psize.nr0+1+(n-1)*(N(1)+1):psize.nr0+n*(N(1)+1));
     solution.final.regulated{2}(:,n) = ifcht(coeff_final_regulated_local);
    end

     for n=1:psize.nry
         nprev=psize.nr0+psize.nrx*(N(1)+1);
         ni=nprev+1+(n-1)*(N(2)+1);
         nii=nprev+n*(N(2)+1);
     coeff_final_regulated_local=coeff_final_regulated(ni:nii);
     solution.final.regulated{3}(:,n) = ifcht(coeff_final_regulated_local);
     end

 % Infinite-dimensional 2D outputs
     for n=1:psize.nr2
         nprev=psize.nr0+(N+1)*[psize.nrx;psize.nry];
     ni=nprev+1+(n-1)*prod(N+1);
     nii=nprev+n*prod(N+1);
     coeff_final_regulated_local=coeff_final_regulated(ni:nii);
     coeff_final_regulated_local_2D=reshape(coeff_final_regulated_local,N(1)+1,N(2)+1);
     solution.final.regulated{4}(:,:,n)=fcgltran2d(coeff_final_regulated_local_2D,0);
     end % nr2
  
    end % nr_sum>0

    % Reconstruction of observed outputs

    if(no_sum>0)

    coeff_final_observed=C2op*solcoeff.final;
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

   % Infinite-dimensional 1D outputs

    for n=1:psize.nox
     coeff_final_observed_local=coeff_final_observed(psize.no0+1+(n-1)*(N(1)+1):psize.no0+n*(N(1)+1));
     solution.final.observed{2}(:,n) = ifcht(coeff_final_observed_local);
    end

     for n=1:psize.noy
         nprev=psize.no0+psize.nox*(N(1)+1);
         ni=nprev+1+(n-1)*(N(2)+1);
         nii=nprev+n*(N(2)+1);
     coeff_final_observed_local=coeff_final_observed(ni:nii);
     solution.final.observed{3}(:,n) = ifcht(coeff_final_observed_local);
     end

 % Infinite-dimensional 2D outputs
     for n=1:psize.no2
         nprev=psize.no0+(N+1)*[psize.nox;psize.noy];
     ni=nprev+1+(n-1)*prod(N+1);
     nii=nprev+n*prod(N+1);
     coeff_final_observed_local=coeff_final_observed(ni:nii);
     coeff_final_observed_local_2D=reshape(coeff_final_observed_local,N(1)+1,N(2)+1);
     solution.final.observed{4}(:,:,n)=fcgltran2d(coeff_final_observed_local_2D,0);
     end % nr2

    end  % no_sum
     
     
%----------------------------------------   
% We now transform time-dependent solution, which is available for BDF
% scheme only (opts.IntScheme=1)
%----------------------------------------
          
if (opts.intScheme==1&opts.tf~=0)
    
    % Define ODE solution and temporal stamps array
     solution.timedep.dtime=solcoeff.timedep.dtime;
     if (n0>0)
     solution.timedep.primary{1}=solcoeff.timedep.coeff(1:n0,:);
     else
     solution.timedep.primary{1}=[];
     end
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
             if (uinput.wsep{k})
        wvec(k,1)=uinput.w{k}(tt); 
         else
        uinput.wspace{k}=subs(uinput.w{k},st,tt);
        % Relate the disturbance index to its type
        group = find(k <= cum_nwg, 1);
        index=cum_nwp(group-1) +(k-cum_nwg(group-1)-1)*Np(group)+1;  
             switch group
             case 2
        solcoeff.w(index:index+N(1),k)=PIESIM_NonPoly2Mat_cheb(N(1), uinput.wspace{k}, 0, gridall.x);
             case 3
        solcoeff.w(index:index+N(2),k)=PIESIM_NonPoly2Mat_cheb(N(2), uinput.wspace{k}, 0, gridall.y);
             case 4
        solcoeff.w(index:index+prod(N+1)-1,k)=PIESIM_NonPoly2Mat_cheb_2D(N, uinput.wspace{k}, 0, gridall);
             end % swtich group
        wvec(k,1)=1;
             end %  if (uinput.wsep{k})
             end %  for k = 1:numel(uinput.w)
     acheb_p=acheb_p+Twcheb_2PDEstate*solcoeff.w*wvec;
     end % psize.nw>0

     % Define inhomogeneous contributions due to controlled inputs
    if (psize.nu>0)
    for k = 1:numel(uinput.u)
             if (uinput.usep{k})
    uvec(k,1)=uinput.u{k}(tt); 
             else
                  uinput.uspace{k}=subs(uinput.u{k},st,tt);
        % Relate the control input index to its type
        group = find(k <= cum_nug, 1);
        index=cum_nup(group-1) +(k-cum_nug(group-1)-1)*Np(group)+1;  
             switch group
             case 2
        solcoeff.u(index:index+N(1),k)=PIESIM_NonPoly2Mat_cheb(N(1), uinput.uspace{k}, 0, gridall.x);
             case 3
        solcoeff.u(index:index+N(2),k)=PIESIM_NonPoly2Mat_cheb(N(2), uinput.uspace{k}, 0, gridall.y);
             case 4
        solcoeff.u(index:index+prod(N+1)-1,k)=PIESIM_NonPoly2Mat_cheb_2D(N, uinput.uspace{k}, 0, gridall);
             end % swtich group
        uvec(k,1)=1;
             end %  if (uinput.usep{k})
             end %  for k = 1:numel(uinput.u)
     acheb_p=acheb_p+Tucheb_2PDEstate*solcoeff.u*uvec;
 end % psize.nu>0
     

     % Reconstruct solution using inverse Chebyshev transform 

     % Reconstruct 1D states

     for n=1:nx
         ni=n0+1+(n-1)*(N(1)+1);
         nii=n0+n*(N(1)+1);
     acheb_p_local=acheb_p(ni:nii);
     solution.timedep.primary{2}(:,n,ntime) = ifcht(acheb_p_local);
     end % nx

     for n=1:ny
         nprev=n0+nx*(N(1)+1);
         ni=nprev+1+(n-1)*(N(2)+1);
         nii=nprev+n*(N(2)+1);
     acheb_p_local=acheb_p(ni:nii);
     solution.timedep.primary{3}(:,n,ntime) = ifcht(acheb_p_local);
     end % ny

     % Reconstruct 2D states
     for n=1:n2
         nprev=n0+(N+1)*[nx;ny];
         ni=nprev+1+(n-1)*prod(N+1);
         nii=nprev+n*prod(N+1);
     acheb_p_local=acheb_p(ni:nii);
     acheb_local_2D=reshape(acheb_p_local,N(1)+1,N(2)+1);
     solution.timedep.primary{4}(:,:,n,ntime)=fcgltran2d(acheb_local_2D,0);
     end % ns


% Reconstruction of regulated outputs

    if(nr_sum>0)
    coeff_timedep_regulated=C1op*solcoeff.timedep.coeff(:,ntime);
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

 % Infinite-dimensional 1D outputs

    for n=1:psize.nrx
     coeff_timedep_regulated_local=coeff_timedep_regulated(psize.nr0+1+(n-1)*(N(1)+1):psize.nr0+n*(N(1)+1));
     solution.timedep.regulated{2}(:,n,ntime) = ifcht(coeff_timedep_regulated_local);
    end

     for n=1:psize.nry
         nprev=psize.nr0+psize.nrx*(N(1)+1);
         ni=nprev+1+(n-1)*(N(2)+1);
         nii=nprev+n*(N(2)+1);
     coeff_timedep_regulated_local=coeff_timedep_regulated(ni:nii);
     solution.timedep.regulated{3}(:,n,ntime) = ifcht(coeff_timedep_regulated_local);
     end

 % Infinite-dimensional 2D outputs
     for n=1:psize.nr2
         nprev=psize.nr0+(N+1)*[psize.nrx;psize.nry];
     ni=nprev+1+(n-1)*prod(N+1);
     nii=nprev+n*prod(N+1);
     coeff_timedep_regulated_local=coeff_timedep_regulated(ni:nii);
     coeff_timedep_regulated_local_2D=reshape(coeff_timedep_regulated_local,N(1)+1,N(2)+1);
     solution.timedep.regulated{4}(:,:,n,ntime)=fcgltran2d(coeff_timedep_regulated_local_2D,0);
     end % nr2
  
    end % nr_sum>0

    % Reconstruction of observed outputs

    if(no_sum>0)

    coeff_timedep_observed=C2op*solcoeff.timedep.coeff(:,ntime);
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

   % Infinite-dimensional 1D outputs

    for n=1:psize.nox
     coeff_timedep_observed_local=coeff_timedep_observed(psize.no0+1+(n-1)*(N(1)+1):psize.no0+n*(N(1)+1));
     solution.timedep.observed{2}(:,n,ntime) = ifcht(coeff_timedep_observed_local);
    end

     for n=1:psize.noy
         nprev=psize.no0+psize.nox*(N(1)+1);
         ni=nprev+1+(n-1)*(N(2)+1);
         nii=nprev+n*(N(2)+1);
     coeff_timedep_observed_local=coeff_timedep_observed(ni:nii);
     solution.timedep.observed{3}(:,n,ntime) = ifcht(coeff_timedep_observed_local);
     end

 % Infinite-dimensional 2D outputs
     for n=1:psize.no2
         nprev=psize.n0+(N+1)*[psize.nox;psize.noy];
     ni=nprev+1+(n-1)*prod(N+1);
     nii=nprev+n*prod(N+1);
     coeff_timedep_observed_local=coeff_timedep_observed(ni:nii);
     coeff_timedep_observed_local_2D=reshape(coeff_timedep_observed_local,N(1)+1,N(2)+1);
     solution.timedep.observed{4}(:,:,n,ntime)=fcgltran2d(coeff_timedep_observed_local_2D,0);
     end % nr2

    end  % no_sum
     

    
     end % ntime

end %  if (opts.intScheme==1&opts.tf~=0)
     
%  
