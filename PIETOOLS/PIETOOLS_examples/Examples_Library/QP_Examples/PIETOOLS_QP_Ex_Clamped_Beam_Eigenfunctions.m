%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PIETOOLS_QP_Ex_Clamped_Beam_Eigenfunctions.m     PIETOOLS 2026
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Clamped beam eigenfunctions computed with qp_opvar. This file is meant to
% be run from the Matlab editor window. For the eigenvalue problem
%
%       u_ss = lambda*u_ssss,       u(0) = u(b) = u_s(0) = u_s(b) = 0,
%
% on [0,b], and given eigenvalues lambda, it computes the eigenfunctions u
% with the least squares solver qp_opvar. The eigenvalues are precomputed
% and defined manually in Step 1 (six values for b = 1, scaled by b^2); they
% are the largest-magnitude eigenvalues of this problem, the largest being
% lambda_1 = -(b/(2*pi))^2 (buckling load of the clamped column; the
% constant of the LPI ||u_s||^2 <= C ||u_ss||^2 is C = -lambda_1, see
% example 15 of examples_FI_library_PIETOOLS). The exact eigenvalues are
% lambda = -(b/beta)^2, with beta = 2*pi*j or beta = 2*x, tan(x) = x.
%
% Method.
%  Fundamental state:  v = u_ssss in L_2[0,b]  (a polynomial of degree 16).
%  PIE of the PDE u_t = u_ss + u_ssss with the same boundary conditions:
%       T v_t = A v, where T = PIE.T gives u = T v (it integrates v four times
%       with the clamped boundary conditions), and A = P + I with P v = u_ss,
%       so P = PIE.A - I. The eigenvalue problem reads (P - lambda*I) v = 0.
%  The least squares problem  min ||(P - lambda*I) X||^2  s.t.  X(0) = 1  is
%       solved by qp_opvar itself, with the constraint built into the unknown,
%       X = 1 + s*Z with Z a polynomial of degree 15:
%
%           min_Z ||A Z - B||^2,   A = (P - lambda*I) S,   B = -(P - lambda*I) 1,
%
%       where S is the multiplication by s and 1 the constant function (an R -> L_2
%       map). qp_opvar returns Z (assembly 'apply', Legendre basis), X = v is the
%       eigenfunction, u = T v is computed with apply_sopvar and normalized to
%       ||u||_{L2} = 1. The signs of modes 1, 3, 4, 5 are flipped.
%  Only the first 6 modes are used (the polynomial of degree 16 resolves about
%  7): the polynomial cannot represent the oscillations of higher modes.
%
% INPUT
% - none; Step 1 sets b = 1, the number of modes nm = 6 and plotting = true.
%
% OUTPUT
% - none; the workspace receives lams (the eigenvalues used), lam_exact (the
%   exact ones), sg (a grid on [0,b]), Uvals (the eigenfunctions u on the grid,
%   one column per mode, ||u|| = 1) and Vvals (v = u_ssss). The eigenvalues and
%   the clamped boundary values of u are printed, and the eigenfunctions plotted.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PIETOOLS - PIETOOLS_QP_Ex_Clamped_Beam_Eigenfunctions
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
% Step 1: Choose the interval, the number of modes and the eigenvalues
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Length of the interval [0,b]
b = 1;
% Number of eigenfunctions to compute (1 through 6)
nm = 6;
% Plot the eigenfunctions (true) or not (false)
plotting = true;
% Degree of the polynomial v = u_ssss
deg = 16;

% The eigenvalues, defined manually (precomputed for b = 1; the eigenvalues of
% u_ss = lambda*u_ssss scale with b^2); the first six, most negative first
lam_b1 = [-0.025330379244080; -0.012382004040601; -0.006332657311635; ...
          -0.004189125353983; -0.002814560658210; -0.002102692882017];
if nm < 1 || nm ~= round(nm),   error('The number of modes nm must be a positive integer.');   end
if nm > numel(lam_b1)
    error('Only %d eigenvalues are defined; reduce nm.',numel(lam_b1));
end
lams = b^2*lam_b1(1:nm);

% The exact eigenvalues, for comparison: beta = 2*pi*j (cosine family) and beta = 2*x, tan(x) = x
xt = arrayfun(@(j) fzero(@(z) tan(z)-z,[j*pi+0.01,(j+0.5)*pi-0.01]),1:nm);
bb = sort([2*pi*(1:nm), 2*xt]);     lam_exact = -b^2./bb(1:nm).^2;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Step 2: Declare the variables and the PIE
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Spatial variable s1, its dummy variable s1_dum, time t; interval [0,b]
pvar s1 s1_dum t
dom = [0,b];

% PDE variables on which the operators act: y in L2[0,b] (used to declare
% operators on L2[0,b]), x0 in R (used for maps from R to L2) and u4 in
% H4[0,b], the state of the PDE below. Operators are expressed in terms of
% these variables and extracted with `get_op`.
clear stateNameGenerator
x0 = pde_var(1);
y  = pde_var(s1,dom);
u4 = pde_var(s1,dom,4);

% PIE of u_t = u_ss + u_ssss with the four clamped boundary conditions
% (u(0) = u(b) = u_s(0) = u_s(b) = 0)
PDE = [diff(u4,t) == diff(u4,s1,2) + diff(u4,s1,4);
       subs(u4,s1,0) == 0;             subs(u4,s1,b) == 0;
       subs(diff(u4,s1,1),s1,0) == 0;  subs(diff(u4,s1,1),s1,b) == 0];
PIE = convert(PDE,'pie');

% Operators acting on the fundamental state v = u_ssss in L2[0,b]:
% identity operator on L2[0,b]
Id = get_op(y);
% P v = u_ss -- P = PIE.A - I; the eigenvalue problem u_ss = lambda*u_ssss reads (P - lambda*I) v = 0
P = PIE.A - Id;
% T v = u    -- T = PIE.T, integrates v four times with the clamped boundary conditions
T = PIE.T;

% Operators used to build the normalization v(0) = 1 into the unknown, X = v = 1 + s*Z:
% (S z)(s) = s*z(s) -- multiplication operator by s (L2 -> L2)
Sm = get_op(s1*y);
% E1 c = c*1 -- the constant function 1, as a map from R to L2 (R -> L2)
E1 = get_op(x0,y);
% P*S and P*1, the two pieces of (P - lambda*I)(1 + s*Z) = P*1 + P*S*Z - lambda*(1 + s*Z)
PS = P*Sm;   PE = P*E1;

% X = 1 + s*Z with Z of degree deg-1: qp_opvar solves the least squares problem in Z
qopts = struct('deg',deg-1,'basis','legendre','assembly','apply','solver','backslash');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Step 3: Compute the eigenfunctions
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Solve for each eigenfunction and evaluate on a grid sg
sg = linspace(0,b,400)';
Uvals = zeros(400,nm);   Vvals = zeros(400,nm);   bcs = zeros(nm,2);
sgn = ones(1,nm);   sgn(intersect([1 3 4 5],1:nm)) = -1;
for k = 1:nm
    % QP by qp_opvar: min_Z ||(P - lambda*I)(1 + s*Z)||^2, i.e. ||A Z - B||^2 with
    % A = (P - lambda*I)*S and B = -(P - lambda*I)*1; then v = X = 1 + s*Z, X(0) = 1
    Zk = qp_opvar(PS - lams(k)*Sm,lams(k)*E1 - PE,qopts);
    blk = Zk.C{1,1};
    cz = full(blk.params{1});   ez = blk.ZL{1}(:);     % monomial coefficients of Z, exponents 0..deg-1
    c = zeros(deg+1,1);   c(1) = 1;                    % v = 1 + s*Z in monomials of s
    c(ez+2) = c(ez+2) + cz;
    vpoly = 0*s1;
    for j = 1:numel(c),  vpoly = vpoly + c(j)*s1^(j-1);  end
    upoly = apply_sopvar(T.C{1,1},vpoly);              % u = T v
    Vvals(:,k) = polyval(flipud(c),sg);
    Uvals(:,k) = arrayfun(@(z) double(subs(upoly,s1,z)),sg);
    dupoly = diff(upoly,s1);                           % end values of u_s (clamped: 0)
    bcs(k,:) = [double(subs(dupoly,s1,0)) double(subs(dupoly,s1,b))];
end

% Normalize u to ||u||_{L2} = 1 and flip the signs of modes 1, 3, 4, 5
nrm = sqrt(trapz(sg,Uvals.^2));
Uvals = Uvals./nrm;   bcs = bcs./nrm(:);              % ||u||_{L2} = 1
Uvals = Uvals.*sgn;   Vvals = Vvals.*sgn;

% Report: the eigenvalues, and how well u satisfies the clamped conditions (all should be 0)
fprintf('\n  k   lambda (used)   lambda (exact)    |u(b)|   |u_s(0)|  |u_s(b)|  (clamped: all 0)\n');
for k = 1:nm
    fprintf(' %2d  %14.8f  %14.8f   %9.1e %9.1e %9.1e\n',k,lams(k),lam_exact(k), ...
        abs(Uvals(end,k)),abs(bcs(k,1)),abs(bcs(k,2)));
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Step 4: Plot the eigenfunctions u = T v
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if plotting
    co = [0 0.447 0.741; 0.85 0.325 0.098; 0.929 0.694 0.125; 0.494 0.184 0.556; 0.466 0.674 0.188; 0.301 0.745 0.933; 0.635 0.078 0.184];
    figure; hold on;
    for k = nm:-1:1
        plot(sg,Uvals(:,k),'Color',co(mod(6-k,7)+1,:),'LineWidth',2);
    end
    xlim([0 b]);  xlabel('s');  title('Computed eigenfunctions, u = T v  (||u||_{L2} = 1)');
end

% The largest eigenvalue is lambda_1 = -(b/(2*pi))^2
fprintf('\n Largest eigenvalue lambda_1 = %.6f (exact -(b/(2*pi))^2 = %.6f).\n',lams(1),-(b/(2*pi))^2);
