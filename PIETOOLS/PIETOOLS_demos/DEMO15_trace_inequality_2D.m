% DEMO15_trace_inequality_2D.m
%
% This document extends DEMO14 to a 2D domain: it illustrates how a trace
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
% x_{s1 s2} in L2 (the 2D version of W^1). No condition is imposed on x:
% instead, states of lower dimension are declared that x determines, and
% the relations between them form a boundary value problem (BVP). Option
% 'int' (below) uses the means
%       x1(s1) = int_0^1 x(s1,s2) ds2,   x2(s2) = int_0^1 x(s1,s2) ds1,
%       x0 = int_0^1 x1(s1) ds1 = int_0^1 x2(s2) ds2,
% and option 'subs' the edge values
%       x1(s1) = x(s1,0),   x2(s2) = x(0,s2),   x0 = x1(0) = x2(0).
% PIETOOLS then returns the fundamental state xf, free of boundary
% conditions, with four components: x0, x1_{s1}, x2_{s2} and x_{s1 s2}. 
% The operator T has a row for each state x, x0, x1, x2.
% That of x gives the operators such that
%               x = H0op*xf,    x_{s1} = H1op*xf,    x_{s2} = H2op*xf.
% Composing H0op with evaluation at each edge gives the edge values,
%               x(0,s2) = E1op*xf,   x(1,s2) = E2op*xf,
%               x(s1,0) = E3op*xf,   x(s1,1) = E4op*xf.
% Both options parameterize this space and give the same constant.
% We can then pose the trace inequality as an LPI
%   min_{gam}   gam,
%      s.t.     gam*(H0op'*H0op + H1op'*H1op + H2op'*H2op) - sum_k Ekop'*Ekop ≥ 0.
%
% The 1D version of this inequality (DEMO14) has the constant 2.164 on
% [0,1] for x(0)^2 + x(1)^2. The optimal constant on the unit square is
%   c = sqrt(2)/tanh(1/(2*sqrt(2))) = 4.16529,
% attained by x = cosh((s1-1/2)/sqrt(2))*cosh((s2-1/2)/sqrt(2)), which
% solves -x_{s1s1}-x_{s2s2}+x = 0 with x_{n} = x/c on the boundary (n the
% outward normal). The same construction applies to more
% variables, at a higher cost: the certificate grows quickly with the
% number of variables and with the degrees dD, dw.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PIETOOLS - DEMO15
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
% DJ, 11/25/2026: Use ODE and PDE states to declare BVP;

clc; clear; clear stateNameGenerator;
echo on

% =============================================
% === Declare the operators of interest

% % Choose the states of lower dimension declared with x, 'int' or 'subs':
% %   'int'  - means of x; the fundamental state is
% %            [int int x; int x_{s1} ds2; int x_{s2} ds1; x_{s1 s2}];
% %            certificate degree dD = 0 suffices (about 18600 decision
% %            variables, bound 4.1654)
% %   'subs' - edge values of x; the fundamental state is
% %            [x(a,c); x_{s1}(s1,c); x_{s2}(a,s2); x_{s1 s2}];
% %            needs dD = 1 (about 18600 decision variables, bound 4.1654;
% %            dD = 0 gives 4.449)
option = 'int';

% % Declare the BVP
% %     x in W^{[1,1]}([a,b]x[c,d])
a = 0;  b = 1;  c = 0;  d = 1;
pvar s1 s2
x = pde_var([s1;s2],[a,b;c,d],[1;1]);
switch option
    case 'subs'
        % Edge values declared as states
        x0 = pde_var();                     % x(a,c)
        x1 = pde_var(s1,[a,b],1);           % x(s1,c)
        x2 = pde_var(s2,[c,d],1);           % x(a,s2)
        BVP = [subs(x,s1,a)==x2;
               subs(x,s2,c)==x1;
               subs(x1,s1,a)==x0;
               subs(x2,s2,c)==x0];
    case 'int'
        % Mean values declared as states
        x0 = pde_var();                     % int_a^b int_c^d x(s1,s2) ds2 ds1
        x1 = pde_var(s1,[a,b],1);           % int_c^d x(s1,s2) ds2
        x2 = pde_var(s2,[c,d],1);           % int_a^b x(s1,s2) ds1
        BVP = [int(x,s1,[a,b])==x2;
               int(x,s2,[c,d])==x1
               int(x1,s1,[a,b])==x0;
               int(x2,s2,[c,d])==x0];
    otherwise
        error("option must be 'int' or 'subs'.");
end

% % Parameterize all x in H1:  x = H0op*xf, and its first derivatives
% xf = [x0; x1_{s1}; x2_{s2}; x_{s1 s2}]. T has a row for each state
% (x, x0, x1, x2), but only that of x is used
BVP_sol = convert(BVP);
H0op = BVP_sol.T(1,:);          % (H0op*xf) = x
H1op = diff(H0op,s1);           % (H1op*xf) = x_{s1}
H2op = diff(H0op,s2);           % (H2op*xf) = x_{s2}

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
% % (option 'subs' needs the certificate degree dD = 1 to reach the same
% % value)
P = gam*(H0op'*H0op + H1op'*H1op + H2op'*H2op) ...
    - E1op'*E1op - E2op'*E2op - E3op'*E3op - E4op'*E4op;
if strcmp(option,'int')
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

fprintf(['\n If successful, ',num2str(trace_constant,6),' is an upper bound on c in int_{bdry} x^2 <= c(||x||^2+||x_s1||^2+||x_s2||^2) on [0,1]^2.\n']) % DJ (CC), 11/25/2026
fprintf([' The optimal value is sqrt(2)/tanh(1/(2*sqrt(2))) = 4.16529.\n']);
