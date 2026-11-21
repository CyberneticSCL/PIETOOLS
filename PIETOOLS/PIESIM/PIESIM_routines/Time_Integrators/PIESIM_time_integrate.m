%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PIESIM_time_integrate.m    PIETOOLS 2025
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Performs time-advancement of a discretized PIE with the scheme specified
% by "opts.intScheme"=1 or 2
% opts.intScheme=1 - Backward Differentiation Formula
% opts.intScheme=2 - Analytical (symbolic) integration of the matrix
% exponential formula
%
% Inputs:
% psize - size of the PIE problem: N, dim, n0, nw, nu (nw0, nwx, nwy, nw2,
% nu0, nux, nuy, nu2 in 2D)
% opts - options for temporal scheme parameters
%    ('struct' with fields intScheme, Norder, dt, Nsteps, tf, N, ploteig)
% uinput   - user-defined boundary inputs
%    ('struct' with fields w, wdot, u, udot, wsep, usep, wspace, uspace)
% coeff - Chebyshev coefficients of initial conditions and forcing terms, if any
%    ('struct' with fields acheb_f0, w, u)
% gridall - physical grids for solution states of different degree of
% differentiability
% Dop - discretized PIE operators
%    ('struct' with fields Acheb, Atotal = inv(T)*A, Tcheb_inv, B1cheb,
%    B2cheb, Twcheb, Tucheb)
%
% Output:
% solcoeff - contains:
% 1) solcoeff.final - Chebyshev coefficients of the final
% solution, ODE+PDE states
% 2) solcoeff.tf - final time of the solution
%----------------------------------------
% NOTE: the following options are only avaialble for BDF scheme (opts.intScheme=1),
% which performs a temporal discretiztion of the PIE equation. The other
% scheme solves for the temporal integral in a convolution form directly,
% analytically (opts.intScheme=2),
% thus a temporal history of solution is not available with this scheme, and only the final answer is.
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
% 6) solcoeff.timedep.coeff - coefficients of all states, one column per
% time step; solcoeff.w, solcoeff.u - coeff.w, coeff.u at the final time
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PIETOOLS - PIESIM_time_integrate
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
% Initial coding YP  - 9_18_2021
% YP 6/16/2022 - added stability check
% YP, MMP, 11/21/2026: Changes since PIETOOLS 2025: new input gridall,
%   passed to the BDF integrator, also when symbolic integration falls back
%   to BDF2.

function solcoeff=PIESIM_time_integrate(psize, opts, uinput, coeff, gridall, Dop)
N=psize.N;

% Stability of the time scheme for the eigenvalues of Atotal
PIESIM_stability_check(opts, Dop.Atotal);
    
switch opts.intScheme 
     case 1    
     % BDF of order opts.Norder
     solcoeff=PIESIM_int_bdf(psize, opts, uinput, coeff, gridall, Dop);
     case 2
         % Symbolic integration through Atotal = V*D/V; if Atotal is not
         % diagonalizable, BDF2 with the time step tf/floor(tf/1e-3)
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
%            solcoeff=PIESIM_int_bdf(psize, opts, uinput, coeff, Dop);      % YP, MMP, 11/21/2026 (was)
             solcoeff=PIESIM_int_bdf(psize, opts, uinput, coeff, gridall, Dop); % YP, MMP, 11/21/2026
         end
 end
