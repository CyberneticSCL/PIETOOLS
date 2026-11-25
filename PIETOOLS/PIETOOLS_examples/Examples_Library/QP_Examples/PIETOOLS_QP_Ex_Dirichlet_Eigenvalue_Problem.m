%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PIETOOLS_QP_Ex_Dirichlet_Eigenvalue_Problem.m     PIETOOLS 2026
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Dirichlet eigenvalue problem solved with qp_copvar. This file is meant to
% be run from the Matlab editor window. It computes the eigenfunctions of
%
%       x = lambda*x_ss,        x(0) = x(b) = 0,        on [0,b],
%
% with the least squares solver qp_copvar, for the known eigenvalues
% lambda_n = -(b/(n*pi))^2, and compares them with the exact eigenfunctions
% x_n = sqrt(2/b)*sin(n*pi*s/b), n = 1, 2, ...
%
% Method.
%  Fundamental state:  v = x_ss in L_2[0,b].
%  PIE of the PDE u_t = u_ss, u(0) = u(b) = 0:  T v_t = A v with A = I (since
%       u_ss = v) and T = PIE.T the Green's operator of v = x_ss with
%       x(0) = x(b) = 0, so x = T v. The problem x = lambda*x_ss reads
%       (T - lambda*I) v = 0.
%  Eigenvalues:  the known eigenvalues lambda_n = -(b/(n*pi))^2 are used; they
%       are not computed here.
%  Eigenfunctions:  v = x_ss is a multiple of sin(n*pi*s/b), so v(0) = 0, and
%       the normalization is v(0) = 0, v'(0) = 1: v = s + s^2*Z with Z a
%       polynomial of degree 14. For each lambda, qp_copvar solves the least
%       squares problem  min_Z ||(T - lambda*I)(s + s^2*Z)||^2 = ||A Z - B||^2
%       with A = (T - lambda*I) S^2 and B = -(T - lambda*I) s (S the
%       multiplication by s, and s the function s as a map from R to L_2);
%       then x = T v = lambda*v is computed with apply_sopvar and normalized to
%       ||x||_{L2} = 1.
%  Only the first modes (about 6 at the polynomial degree 16 used) are
%  resolved: the polynomial cannot represent the oscillations of
%  sin(n*pi*s/b) for larger n.
%
% INPUT
% - none; Step 1 sets b = 1, the number of modes nm = 6 and plotting = true.
%
% OUTPUT
% - none; the workspace receives lams (the eigenvalues used), sg (a grid on [0,b]),
%   Xvals (the eigenfunctions x on the grid, one column per mode, ||x|| = 1)
%   and Vvals (v = x_ss). The errors of the eigenfunctions are printed, and
%   the eigenfunctions plotted.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PIETOOLS - PIETOOLS_QP_Ex_Dirichlet_Eigenvalue_Problem
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
% AT, 11/25/2026: Initial coding. 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

close all; clc; clear;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Step 1: Choose the interval and the number of modes
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Length of the interval [0,b]
b = 1;
% Number of eigenfunctions to compute (1 through 6)
nm = 6;
% Plot the eigenfunctions (true) or not (false)
plotting = true;

if nm < 1 || nm ~= round(nm),   error('The number of modes nm must be a positive integer.');   end
if nm > 6,   error('At most 6 modes are resolved at degree 16; reduce nm.');   end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Step 2: Declare the variables and the PIE
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Spatial variable s1, its dummy variable s1_dum, time t; interval [0,b]
pvar s1 s1_dum t
dom = [0,b];

% PDE variables on which the operators act: y in L2[0,b] (used to declare
% operators on L2[0,b]), x0 in R (used for maps from R to L2) and u2 in
% H2[0,b], the state of the PDE below. Operators are expressed in terms of
% these variables and extracted with `get_op`.
clear stateNameGenerator
x0 = pde_var(1);
y  = pde_var(s1,dom);
u2 = pde_var(s1,dom,2);

% PIE of u_t = u_ss with u(0) = u(b) = 0
PDE = [diff(u2,t) == diff(u2,s1,2);
       subs(u2,s1,0) == 0;   subs(u2,s1,b) == 0];
PIE = convert(PDE,'pie');

% Operators acting on the fundamental state v = x_ss in L2[0,b]:
% T v = x -- the Green's operator of v = x_ss with x(0) = x(b) = 0 (PIE.T); the eigenvalue
% problem x = lambda*x_ss reads (T - lambda*I) v = 0
T = PIE.T;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Step 3: The eigenvalues
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% The known eigenvalues lambda_n = -(b/(n*pi))^2 (not computed)
lams = -(b./(pi*(1:nm)')).^2;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Step 4: Compute the eigenfunctions
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% v = x_ss is a multiple of sin(n*pi*s/b), so v(0) = 0 and the normalization
% is v(0) = 0, v'(0) = 1: v = s + s^2*Z with Z a polynomial of degree 14.
% For each lambda, qp_copvar solves the least squares problem
%   min_Z ||(T - lambda*I)(s + s^2*Z)||^2 = ||A Z - B||^2,
%   A = (T - lambda*I) S^2,   B = -(T - lambda*I) s,
% (S the multiplication by s, s the function s as an R -> L2 map); then x = T v = lambda*v
% is computed with apply_sopvar and normalized to ||x||_{L2} = 1, with x'(0) > 0.

% (S z)(s) = s*z(s) -- multiplication operator by s (L2 -> L2)
Sm = get_op(s1*y);
% Es c = c*s -- the function s, as a map from R to L2 (R -> L2)
Es = get_op(s1*x0,y);
% T*S^2 and T*s, the pieces of (T - lambda*I)(s + s^2*Z) = T*s + T*S^2*Z - lambda*(s + s^2*Z)
TS2 = T*(Sm*Sm);   TEs = T*Es;
qopts = struct('deg',14,'basis','legendre','assembly','apply','solver','backslash');

% Solve for each eigenfunction and evaluate on a grid sg; Xex is the exact eigenfunction for comparison
sg = linspace(0,b,400)';
Xvals = zeros(400,nm);   Vvals = zeros(400,nm);   Xex = zeros(400,nm);
for k = 1:nm
    Zk = qp_copvar(TS2 - lams(k)*(Sm*Sm),lams(k)*Es - TEs,qopts);
    blk = Zk.C{1,1};
    cz = full(blk.params{1});   ez = blk.ZL{1}(:);     % monomial coefficients of Z
    c = zeros(17,1);   c(2) = 1;                       % v = s + s^2*Z
    c(ez+3) = c(ez+3) + cz;
    vpoly = 0*s1;
    for j = 1:numel(c),  vpoly = vpoly + c(j)*s1^(j-1);  end
    xpoly = apply_sopvar(T.C{1,1},vpoly);              % x = T v
    Vvals(:,k) = polyval(flipud(c),sg);
    Xvals(:,k) = arrayfun(@(z) double(subs(xpoly,s1,z)),sg);
    Xex(:,k) = sqrt(2/b)*sin(k*pi*sg/b);               % exact, ||x|| = 1, x'(0) > 0
end

% Normalize x to ||x||_{L2} = 1 and take the sign of the exact eigenfunction
Xvals = Xvals./sqrt(trapz(sg,Xvals.^2));
sg0 = sign(Xvals(40,:) .* Xex(40,:));   Xvals = Xvals.*sg0;   Vvals = Vvals.*sg0;   % sign of the exact one

% Report the error of the eigenfunctions
fprintf('\n  k   max |x - exact|   (x normalized to ||x|| = 1, exact x = sqrt(2/b) sin(k*pi*s/b))\n');
for k = 1:nm
    fprintf(' %2d   %10.2e\n',k,max(abs(Xvals(:,k)-Xex(:,k))));
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Step 5: Plot the eigenfunctions x = T v
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if plotting
    co = [0 0.447 0.741; 0.85 0.325 0.098; 0.929 0.694 0.125; 0.494 0.184 0.556; 0.466 0.674 0.188; 0.301 0.745 0.933; 0.635 0.078 0.184];
    figure; hold on;
    for k = nm:-1:1
        plot(sg,Xvals(:,k),'Color',co(mod(k-1,7)+1,:),'LineWidth',2);
    end
    xlim([0 b]);  xlabel('s');  title('Computed eigenfunctions, x = T v  (||x||_{L2} = 1)');
end
