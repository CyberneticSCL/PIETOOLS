%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PIESIM_time_integrate.m    PIETOOLS 2021b
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Performs time-advancement of a discretized PIE with the scheme specified
% by "opts.intScheme"=1 or 2
% opts.intScheme=1 - Backward Differentiation Formula
% opts.intScheme=2 - Analytical (symbolic) integration of the matrix
% exponential formula
%
% Inputs:
% psize - size of the PIE problem: nw, nu, nf, nx 
% opts - options for temporal scheme parameters
% uinput   - user-defined boundary inputs
% coeff - Chebyshev coefficients of initial conditions and forcing terms, if any
% Dop - discretized PIE operators
%
% Output:
% solcoeff - contains:
% 1) solcoeff.acheb_f - Chebyshev coefficients of the final
% solution, ODE+PDE states
% 2) solcoeff.tf - final time of the solution
%----------------------------------------
% NOTE: the following three options are only avaialble for BDF scheme (opts.intScheme=1),
% which performs a temporal discretiztion of the PIE equation. The other two
% schemes solve for the temporal integral in a convolution form directly,
% either analytically (opts.intScheme=3) or with approximate Gauss inegration (opts.intScheme=2), 
% thus a temporal history of solution is not available with these two schemes, and only the final answer is. 
%-----------------------------------------
% 3) solcoeff.timedep.ode -
% time-dependent Chebyshev coefficients of the ODE states of PIE system for
% output (for ODE states, Chebyshev coefficients are equal to the solution
% of the states, since the ODE states are not spatially dependent)
% 4) solcoeff.timedep.pde -
% time-dependent Chebyshev coefficients of the PDE states of PIE system for output 
% (for PDE states, an inverse Fourier transform needs to be performed to recover physical solution from its )
% Chebyshev coefficients)
% 5)solcoeff.timedep.dtime - temporal stamps (discrete time values) of the time-dependent solution
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% If you modify this code, document all changes carefully and include date
% authorship, and a brief description of modifications
%
% Initial coding YP  - 9_18_2021

function solcoeff=PIESIM_time_integrate(psize, opts, uinput, coeff, Dop)
N=psize.N;
switch opts.intScheme 
     case 1    
     solcoeff=PIESIM_int_bdf(psize, opts, uinput, coeff, Dop);
     case 2
         [Dop.V,Dop.D] = eig(Dop.Atotal);
         Nsize=size(Dop.Atotal,1);
         if (rank(Dop.V)==Nsize)
             if (N>20)
                 disp('Note: symbolic integration with high N might take a long time');
                 disp('For faster results, use Backward Differentiation Scheme'); 
             end
         solcoeff=PIESIM_int_symbolic(psize, opts, uinput, coeff, Dop);
         else
             disp("Error: matrix Atotal is not diagonalizabe, and analytical integration is not possible"); 
             disp("Resorting to a default scheme of BDF2");
             opts.Norder=2;
             dt=1e-3;
             opts.Nsteps=floor(opts.tf/dt);
             opts.dt=opts.tf/opts.Nsteps;
             solcoeff=PIESIM_int_bdf(psize, opts, uinput, coeff, Dop);
         end
 end