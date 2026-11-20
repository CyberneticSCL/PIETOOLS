%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PIESIM_discretize_icf.m     PIETOOLS 2025
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Perform discretization of initial conditions and forcing functions 
% The initial condition of a PDE state differentiable p times, sampled on
% its N-p+1 Chebyshev-Gauss-Lobatto nodes, gives its N-p+1 Chebyshev
% coefficients; the spatial part of an input on [a,b] gives N+1.
%
% Inputs:
% 1) uinput - user-defined boundary inputs, forcing and initial conditions
%    ('struct' with fields a, b (the domain [a,b]); ic, the initial
%    conditions (ODE states first, constants; then the fundamental PDE
%    states, sym in sx); optionally wspace, uspace, the spatial parts of
%    the inputs (sym in sx), and Bpw_nonpol, the non-polynomial forcing)
% 2) psize - size of the problem
%    ('struct' with fields N, n0, n, nw, nu, nw0, nwx, nu0, nux)
% 3) gridall - cell array containing grids of different degrees of differentiability
%
% Outputs:
% 1) coeff - Chebyshev coefficients for initial conditions and forcing functions
%    ('struct': acheb_f0, column, the ODE initial values, then the
%    coefficients of each PDE state; w (u), (nw0+(N+1)*nwx) x nw
%    ((nu0+(N+1)*nux) x nu), column k a unit entry for a
%    finite-dimensional input, the N+1 coefficients of its spatial part
%    for an input on [a,b]; w = 1 (u = 1) without uinput.wspace
%    (uinput.uspace))
% 2) B1_nonpol - contribution to PIE B1 operator arising from non-polyn0mial forcing
%    ([] without uinput.Bpw_nonpol or if nw = 0)
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PIETOOLS - PIESIM_discretize_icf
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
% Initial coding YP  11_05_2021
% YP, 6/16/2022: Renamed uinput.B21_nonpol to uinput.Bpw_nonpol, updated description of
% outputs
% YP, 6/1/2025: Added functionality to support infinite-dimensional disturbances
% through parser 
% YP, 1/7/2026 - changed the notaiton of ic.PDE to ic.PIE for PIE simulation
% YP, 2/17/2026: Changed all initial conditions to uinput.ic format
% MMP, 11/05/2026: No non-polynomial forcing without disturbances (nw=0).
%   The batch and legacy terms input checks set nw=0 when uinput.w is not
%   given but keep uinput.Bpw_nonpol with its nw columns, so the sum
%   Dop.B1cheb+B1_nonpol in PIESIM_discretize_all failed (incompatible sizes).
% MMP, 11/19/2026: Release comment pass (comment-only).

function [coeff,B1_nonpol]=PIESIM_discretize_icf(uinput,psize,gridall);

% Define local variables

n0=psize.n0;
nw=psize.nw;
nu=psize.nu;
nx=sum(psize.n);
a=uinput.a;
b=uinput.b;
N=psize.N;

% Define degree of smoothness p

p = repelem(0:length(psize.n)-1, psize.n);

% Define initial conditions on states on the physical gridall for states.
% var_f denotes the value of the solution variable of the fundamental
% states.
% Define Chebyshev coefficients of initial conditions in the same loop
% acheb_f den0tes the Chebyshev coefficients of the fundamental states
% Each state vector array coefficients are arranged into a global column
% vector

ic=uinput.ic;

for i=1:nx
     acheb=fcht(double(subs(ic(n0+i),gridall{p(i)+1})));
     acheb_glob{i}=reshape(acheb, [], 1);
     clear('acheb');
end


% Concatenate coefficients of all states into a single column vector

 acheb_f0=cat(1, acheb_glob{:});

% Add initial conditions on ODE states to the front of initial conditions
if (n0>0)
    if (size(ic,1)>1) 
        ic=ic';
    end
acheb_f0=cat(1,double(ic(1:n0))',acheb_f0);
end

coeff.acheb_f0=acheb_f0;


% Discretize matrix operator for non-polyn0mial in space forcing

% Bpw_nonpol multiplies w, so nw=0 means no forcing, also when uinput       % MMP, 11/05/2026
% keeps Bpw_nonpol (nw set to 0 for lack of uinput.w).                      % MMP, 11/05/2026
% if isfield(uinput,'Bpw_nonpol')                                           % MMP, 11/05/2026 (was)
if isfield(uinput,'Bpw_nonpol') && nw>0                                     % MMP, 11/05/2026
B1_nonpol = PIESIM_NonPoly2Mat_cheb(N, uinput.Bpw_nonpol, p, gridall);
else
B1_nonpol=[];   
end

coeff.w=1;
coeff.u=1;

% Discretize spatial contribution of u and w disturbances entered in parser
% format as spatially-varying disturbances

if isfield(uinput,'wspace')
     Nforce=psize.nw0+(N+1)*psize.nwx;
     coeff.w=zeros(Nforce,nw);
         k=0;
         index=1;
         % Finite-dimensional disturbances: one unit entry each
         for kk=1:psize.nw0
         k=k+1;
             coeff.w(index,k)=1;
             index=index+1;
         end
         % Disturbances on [a,b]: N+1 coefficients of the spatial part
         for kk=1:psize.nwx
             k=k+1;
             coeff.w(index:index+N,k)=PIESIM_NonPoly2Mat_cheb(N, uinput.wspace(k), 0, gridall);
             index=index+N+1;
         end 
 end

 % Control inputs, as the disturbances
 if isfield(uinput,'uspace')
     Nforce=psize.nu0+(N+1)*psize.nux;
     coeff.u=zeros(Nforce,nu);
         k=0;
         index=1;
         for kk=1:psize.nu0
         k=k+1;
             coeff.u(index,k)=1;
             index=index+1;
         end
         for kk=1:psize.nux
             k=k+1;
             coeff.u(index:index+N,k)=PIESIM_NonPoly2Mat_cheb(N, uinput.uspace(k), 0, gridall);
             index=index+N+1;
         end 
 end












  



