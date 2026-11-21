%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PIESIM_discretize_icf_2D.m     PIETOOLS 2025
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Perform discretization of initial conditions and forcing functions in 2D
% The initial condition of a state differentiable p_x (p_y) times in x (y),
% sampled on its Chebyshev-Gauss-Lobatto grid, gives its N_x-p_x+1
% (N_y-p_y+1) Chebyshev coefficients per direction; the spatial part of an
% input gives N_x+1, N_y+1 or prod(N+1) coefficients.
%
% Inputs:
% 1) uinput - user-defined boundary inputs, forcing and initial conditions
%    ('struct' with fields dom, ic (ODE states first, constants; then the
%    states of x, of y, of (x,y)), w, wspace, wsep, uspace)
% 2) psize - size of the PIE problem: all variables defining the size of the PIE problem
%    ('struct' with fields N = [N_x, N_y], n0, nx, ny, n, nw, nw0, nwx,
%    nwy, nw2, nu0, nux, nuy, nu2)
% 3) gridall - fields containing the following sub-fields:
%  gridall.x - cell array containing gridalls in x direction of different degrees of differentiability
%  gridall.y - cell array containing gridalls in y direction of different degrees of differentiability
%
% Output:
% coeff - Chebyshev coefficients for initial conditions and forcing functions
%    ('struct': acheb_f0, column, the ODE initial values, then the
%    coefficients of the states of x, of y, of (x,y); w, u, one column
%    per input; w = 1 (u = 1) without uinput.wspace (uinput.uspace))
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PIETOOLS - PIESIM_discretize_icf_2D
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
% Initial coding YP  04_16_2024
% YP, 1/7/2026: Changed the notaiton of ic.PDE to ic.PIE for PIE simulation
% YP, 2/17/2026: Changed all initial conditions to uinput.ic format
% MMP, 11/19/2026: Release comment pass (comment-only).
% MMP, 11/20/2026: nu assigned; coeff.u read it unassigned, an error for every
%   2-D simulation with control inputs.

function coeff=PIESIM_discretize_icf_2D(uinput,psize,gridall)

syms sx sy st;

% Define local variables

n0=psize.n0;
nw=psize.nw;
nu=psize.nu;  % number of control inputs, read by coeff.u below             % MMP, 11/20/2026
a=uinput.dom(1,1);
b=uinput.dom(1,2);
c=uinput.dom(2,1);
d=uinput.dom(2,2);
N=psize.N;

        nx=sum(psize.nx,'all');
        ny=sum(psize.ny,'all');
        n2=sum(psize.n,'all');

% Define degree of smoothness p of 2D-1var states

px = repelem(0:length(psize.nx)-1, psize.nx);
py = repelem(0:length(psize.ny)-1, psize.ny);
% Define degree of smoothness p of 2D-2var states

[cols, rows] = meshgrid(0:size(psize.n,2)-1, 0:size(psize.n,1)-1);
p = [repelem(cols(:), psize.n(:))'; repelem(rows(:), psize.n(:))'];

% Define initial conditions on states on the physcal gridall for states.
% var_f denotes the value of the solution variable of the fundamental
% states.
% Define Chebyshev coefficients of initial conditions in the same loop
% acheb_f denotes the Chebyshev coefficients of the fundamental states
% Each state vector array coefficients are arranged into a global column
% vector

ic=uinput.ic;

acheb_glob_x{1}=[];
acheb_glob_y{1}=[];
acheb_glob{1}=[];

% x states only 
for i=1:nx
     acheb=fcht(double(subs(ic(n0+i),gridall.x{px(i)+1})));
     acheb_glob_x{i}=reshape(acheb, [], 1);
     clear('acheb');
end

% y states only 
for i=1:ny
     acheb=fcht(double(subs(ic(n0+nx+i),gridall.y{py(i)+1})));
     acheb_glob_y{i}=reshape(acheb, [], 1);
     clear('acheb');
end

% 2D states (x,y)
for i=1:n2
     acheb=fcgltran2d(double(subs(subs(ic(n0+nx+ny+i),sx,gridall.x{p(1,i)+1}),sy,gridall.y{p(2,i)+1}')),1);
     acheb_glob{i}=reshape(acheb, [], 1);
     clear('acheb');
end


% Concatenate coefficients of all states into a single column vector

 acheb_f0=cat(1, acheb_glob_x{:}, acheb_glob_y{:}, acheb_glob{:});

% Add initial conditions on ODE states to the front of initial conditions
if (n0>0)
    if (size(ic,1)>1) 
        ic=ic';
    end
acheb_f0=cat(1,double(ic(1:n0))',acheb_f0);
end

coeff.acheb_f0=acheb_f0;

% Discretize spatial contribution of u and w disturbnaces

coeff.u=1;
coeff.w=1;

 % Column k of coeff.w: 1 for a finite-dimensional disturbance, else the
 % Chebyshev coefficients of its spatial part (of w at t=0 if w does not
 % separate), stacked as x, y and 2D disturbances
 if isfield(uinput,'wspace')
     Nforce=psize.nw0+(N+1)*[psize.nwx;psize.nwy]+prod(N+1)*psize.nw2;
     coeff.w=zeros(Nforce,nw);
         k=0;
         index=1;
         % Finite-dimensional disturbances
         for kk=1:psize.nw0
         k=k+1;
             coeff.w(index,k)=1;
             index=index+1;
         end
         % Disturbances of x only: N_x+1 coefficients
         for kk=1:psize.nwx
             k=k+1;
             if (~uinput.wsep{k})
                 uinput.wspace{k}=subs(uinput.w{k},st,0);
             end
             coeff.w(index:index+N(1),k)=PIESIM_NonPoly2Mat_cheb(N(1), uinput.wspace{k}, 0, gridall.x);
             index=index+N(1)+1;
         end 
         % Disturbances of y only: N_y+1 coefficients
         for kk=1:psize.nwy
             k=k+1;
              if (~uinput.wsep{k})
                 uinput.wspace{k}=subs(uinput.w{k},st,0);
             end
             coeff.w(index:index+N(2),k)=PIESIM_NonPoly2Mat_cheb(N(2), uinput.wspace{k}, 0, gridall.y);
             index=index+N(2)+1;
         end
         % Disturbances of (x,y): prod(N+1) coefficients
         for kk=1:psize.nw2
             k=k+1;
              if (~uinput.wsep{k})
                 uinput.wspace{k}=subs(uinput.w{k},st,0);
             end
             coeff.w(index:index+prod(N+1)-1,k)=PIESIM_NonPoly2Mat_cheb_2D(N, uinput.wspace{k}, 0, gridall);
             index=index+prod(N+1);
         end
 end

 % Column k of coeff.u: the same for the control inputs
 if isfield(uinput,'uspace')
     Nforce=psize.nu0+(N+1)*[psize.nux;psize.nuy]+prod(N+1)*psize.nu2;
     coeff.u=zeros(Nforce,nu);
         k=0;
         index=1;
         for kk=1:psize.nu0
         k=k+1;
             coeff.u(index,k)=1;
             index=index+1;
         end
         % Control inputs of x only, of y only, of (x,y)
         for kk=1:psize.nux
             k=k+1;
             coeff.u(index:index+N(1),k)=PIESIM_NonPoly2Mat_cheb(N(1), uinput.uspace(k), 0, gridall.x);
             index=index+N(1)+1;
         end 
         for kk=1:psize.nuy
             k=k+1;
             coeff.u(index:index+N(2),k)=PIESIM_NonPoly2Mat_cheb(N(2), uinput.uspace(k), 0, gridall.y);
             index=index+N(2)+1;
         end
         for kk=1:psize.nu2
             k=k+1;
             coeff.u(index:index+prod(N+1)-1,k)=PIESIM_NonPoly2Mat_cheb_2D(N, uinput.uspace(k), 0, gridall);
             index=index+prod(N+1);
         end
 end
