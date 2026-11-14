% DEMO14_trace_inequality_2D.m
%
% This document extends DEMO13 to a 2D domain: it illustrates how a trace
% (boundary value) inequality on the unit square can be verified, and its
% constant found, using PIETOOLS. The trace constant is the smallest value
% c such that
%   int_{BC} x^2  <=  c (||x||^2 + ||x_{s1}||^2 + ||x_{s2}||^2)
% for all x in H1([0,1]^2), where ||.|| is the L2 norm on the square, and
%   int_{BC} x^2 = int_0^1 x(0,s2)^2 + x(1,s2)^2 ds2
%                      + int_0^1 x(s1,0)^2 + x(s1,1)^2 ds1.
% There are no boundary conditions on x. We pose this as an optimization
% program
%   min_{c} c,
%   s.t. int_{BC} x^2 − c(<x,x> + <x_{s1},x_{s1}> + <x_{s2},x_{s2}>) <= 0
%   x in H1.
% To solve, we take x in the space of functions with x_{s1}, x_{s2}, and
% x_{s1 s2} in L2 (the 2D version of W^1), and solve the boundary value
% problem with only trivial conditions along each direction,
%               x(0,s2) = x(0,s2),     x(s1,0) = x(s1,0).
% PIETOOLS then returns the fundamental state xf, free of boundary
% conditions (four components: a scalar, functions of s1, of s2, and of
% (s1,s2); see Dmap), and the operators such that
%               x = H0op*xf,    x_{s1} = H1op*xf,    x_{s2} = H2op*xf.
% Composing H0op with evaluation at each edge gives the edge values,
%               x(0,s2) = E1op*xf,   x(1,s2) = E2op*xf,
%               x(s1,0) = E3op*xf,   x(s1,1) = E4op*xf.
% Alternatively (option 'inputs' below), the values of x on two edges can be
% declared as inputs, x(s1,0) = v1(s1), x(0,s2) = v2(s2), as in DEMO13. The
% fundamental state is then [x(0,0); v1_{s1}; v2_{s2}; x_{s1 s2}], and
%               x = H0op*[w; xf],  with w the inputs.
% Both options parameterize all of H1 and give the same constant.
% We can then pose the trace inequality as an LPI
%   min_{gam}   gam,
%      s.t.     gam*(H0op'*H0op + H1op'*H1op + H2op'*H2op) - sum_k Ekop'*Ekop ≥ 0.
%
% The 1D version of this inequality (DEMO13) has the constant 2.164 on
% [0,1] for x(0)^2 + x(1)^2. The optimal constant on the unit square is
% about 4.165 (finite elements). The same construction applies to more
% variables, at a higher cost: the certificate grows quickly with the
% number of variables and with the degrees dD, dw.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PIETOOLS - DEMO14
%
% Copyright (C)2026  PIETOOLS Team
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
% Initial coding, AT, 11/10/2026: 2D extension of DEMO13_trace_inequality.
% DJ, 11/11/2026: Use `convert` to solve the BVP;

clc; clear; clear stateNameGenerator;
echo on

% =============================================
% === Declare the operators of interest

% % Choose how the state is parameterized:
% %   'mean'   - trivial boundary conditions; the state is
% %              [int int x; int x_{s1} ds2; int x_{s2} ds1; x_{s1 s2}]
% %              SeDuMi solver compute the upper bound 
% %                                             without numerical problems
% %              
% %   'inputs' - boundary values as inputs, x(s1,c) = v1, x(a,s2) = v2; the
% %              state is [x(a,c); v1_{s1}; v2_{s2}; x_{s1 s2}]
% %              SeDuMi solver computes the upper bound, 
% %              but runs into numerical problems
option = 'mean';        % 'mean' or 'inputs'

% % Declare the BVP
% %     x in W^{[1,1]}([a,b]x[c,d])
a = 0;  b = 1;  c = 0;  d = 1;
pvar s1 s2
x = pde_var([s1;s2],[a,b;c,d],[1;1]);
switch option
    case 'mean'
        % Trivial conditions along each direction
        BVP = [subs(x,s1,a)==subs(x,s1,a);
               subs(x,s2,c)==subs(x,s2,c)];
    case 'inputs'
        % Edge values declared as inputs
        v1 = pde_var('in',1,s1,[a,b]);      % v1(s1) = x(s1,c)
        v2 = pde_var('in',1,s2,[c,d]);      % v2(s2) = x(a,s2)
        BVP = [subs(x,s2,c)==v1;
               subs(x,s1,a)==v2];
    otherwise
        error("option must be 'mean' or 'inputs'.");
end

% % Parameterize all x in H1:  x = H0op*xf, and its first derivatives
% xf = [int_0^1 int_0^1 x; int_0^1 x_{s1} ds2; int_0^1 x_{s2} ds1; x_{s1 s2}]
BVP_sol = convert(BVP);
H0op = BVP_sol.T;               % (H0op*xf) = x
H1op = diff(H0op,s1);           % (H1op*xf) = x_{s1}
H2op = diff(H0op,s2);           % (H2op*xf) = x_{s2}

if strcmp(option,'inputs')
    % Prepend the input components w to the state

    % [w, x] = [x(0,0); x_s1(s1, 0);  x_s2(0, s2); x_s1s2(s1,s2)]
    H0op = [BVP_sol.Tw, H0op];     % (H0op*[w;xf]) = x
    H1op = diff(H0op,s1);          % (H1op*[w;xf]) = x_{s1}
    H2op = diff(H0op,s2);          % (H2op*[w;xf]) = x_{s2}
end

% % Edge values of x, as operators from xf to L2 of the other variable
E1op = subs(H0op,s1,a);             % x(a,s2)
E2op = subs(H0op,s1,b);             % x(b,s2)
E3op = subs(H0op,s2,c);             % x(s1,c)
E4op = subs(H0op,s2,d);             % x(s1,d)


% =============================================
% === Declare the LPI

% % Initialize LPI program
prob = lpiprogram(H0op.vars,H0op.dom);

% % Declare decision variables:
% %   gam \in \R
[prob,gam] = lpidecvar(prob,'gam');     % scalar decision variable

% % Set inequality constraints:
% %   gam*(H0op'*H0op + H1op'*H1op + H2op'*H2op) - sum_k Ekop'*Ekop >= 0
% % (the 'inputs' option needs higher certificate degrees, dD = 1, to reach
% % the same value, which costs about 94000 decision variables)
P = gam*(H0op'*H0op + H1op'*H1op + H2op'*H2op) ...
    - E1op'*E1op - E2op'*E2op - E3op'*E3op - E4op'*E4op;
if strcmp(option,'mean')
    prob = lpi_ineq_sop(prob,P,struct('dD',0,'dw',1));
else
    prob = lpi_ineq_sop(prob,P,struct('dD',1,'dw',1));
end

% % Set objective function:
% %   min gam
prob = lpisetobj(prob, gam);

% % Solve LPI and retrieve solution
prob = lpisolve(prob);
trace_constant = double(lpigetsol_sop(prob,gam));


echo off

fprintf(['\n If successful, ',num2str(trace_constant,4),' is an upper bound on c in int_{bdry} x^2 <= c(||x||^2+||x_s1||^2+||x_s2||^2) on [0,1]^2.\n'])
fprintf([' The optimal value is 4.165 (finite elements).\n']);
