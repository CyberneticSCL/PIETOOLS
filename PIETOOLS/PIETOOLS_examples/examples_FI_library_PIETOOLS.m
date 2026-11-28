function varargout = examples_FI_library_PIETOOLS(varargin)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% examples_FI_library_PIETOOLS.m     PIETOOLS 2026
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This library file contains the definition of some common functional
% inequalities drawn from the literature, and computes their constants using
% LPIs. To use one of the examples, call this function with as FIRST input
% the index of the example you wish to use. Please check the table of
% contents below to see which index corresponds to which inequality.
%
% Inputs:
%     [C,info] = examples_FI_library_PIETOOLS(example,option1,option2,...)
%       - example - integer (1-15) corresponding to the selected inequality;
%                   with no input, the function prompts for the index and
%                   the options (separated by spaces or commas)
%       - option# - character array 'name=value' that overrides a default
%                   parameter. Allowed names:
%                     b    - length of the interval [0,b] (default 1),
%                     dD   - degree of the diagonal part of the certificate
%                            (default depends on the example),
%                     dw   - degree of the off-diagonal part of the
%                            certificate (default depends on the example),
%                     t    - interior point of example 8, 0 < t < b
%                            (default b/2),
%                     eps  - regularization of the weight in example 10
%                            (default 1e-2),
%                     tmin - the point t ranges over [tmin, b-tmin] in
%                            example 9, 0 <= tmin < b/2 (default 0.05),
%%
% Table of contents:
%  Poincare-type:   1 Wirtinger (Dirichlet)     
%                   2 Wirtinger (periodic)
%                   3 Poincare (mean zero)      
%                   4 Poincare (x(0) = 0)
%                   5 Generalized Poincare
%  Trace:           6 two endpoints             
%                   7 Sobolev trace (one endpoint)
%                   8 Agmon / 1D Sobolev (interior point t)
%                   9 Agmon / 1D Sobolev, sup over t (t a spatial variable)
%  Hardy-type:      10 Hardy (weight s^2)
%  Operators norms: 11 Volterra operator norm
%                   12 Green operator (Dirichlet)
%                   13 Green operator (mixed)
%                   14 Green operator (hinged beam)
%  Eigenvalues:     15 Clamped beam u_ss = lambda*u_ssss
%
%   The Heisenberg-type and Ehrling (Landau-Kolmogorov) inequalities, which
%   come from products of norms, are not in this library. The eigenvalues and
%   eigenfunctions computed with qp_opvar (Dirichlet problem x = lambda*x_ss,
%   clamped beam eigenfunctions) are in the folder
%   PIETOOLS_examples/Examples_Library/QP_Examples.
%
% Outputs:
%       - C    - double; the upper bound on the constant that the LPI returns,
%       - info - struct with fields
%                   name      - description of the inequality,
%                   constant  - the same bound C,
%                   sharp     - the known sharp constant (NaN if not known),
%                   paper     - constant of Barraza Martinez et al. [1], where
%                               applicable (NaN otherwise),
%                   extra     - struct with example-specific data (empty in
%                               all current examples).
%
% In each example the inequality is written as an operator inequality in
% terms of the free state y = x_s (and, where needed, the value x(0)), so
% that each constant is the solution of an LPI
%       min C   s.t.   (quadratic form of operators with C) >= 0.
%
% If you wish to include a new example in our library, please send it to us
% and we will include it in the next release. Please also include the
% citation information, if available.
%
% For support, contact M. Peet, Arizona State University at mpeet@asu.edu

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PIETOOLS - examples_FI_library_PIETOOLS
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
% AT, 11/16/2026: Initial coding;
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Process the user inputs
[index,opts] = process_inputs(varargin,nargin);

% Print the index; spatial variables s1, s2, dummy s1_dum and time t;
% interval [0,b]; outputs default to NaN (sharp, paper) and an empty struct
% (extra)
fprintf([' --- Extracting functional inequality ', num2str(index),' ---\n']);
pvar s1 s2 s1_dum t
b = opts.b;
dom = [0,b];
sharp = NaN;    paper = NaN;    extra = struct();

% PDE variables on which the operators act: x0 in R, y in L2[0,b], y01 in
% L2[0,1] and y12 = (y1,y2) in L2^2[0,1]; u4 in H4[0,b]
% is the state of the PDE of example 15. Operators are expressed in
% terms of integrals of these variables, and extracted with `get_op`.
clear stateNameGenerator
x0  = pde_var(1);
y   = pde_var(s1,dom);
y01 = pde_var(s1,[0,1]);
y12 = pde_var(2,s1,[0,1]);
u4  = pde_var(s1,dom,4);
switch index
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                           Table of Contents
%
%__________________________________________________________________________
%   Inequality                                         | Parameters     Notes / References
%==========================================================================
%       Poincare-type inequalities
%==========================================================================
    case 1
%   min C
%   s.t. ||x||^2 <= C ||x_s||^2,                             | b = 1        Wirtinger
%        x(0) = x(b) = 0                                   |                Sharp C = (b/pi)^2
%   References: the optimal Poincare constant (b/pi)^2 is computed with SOS programs (SOSTOOLS) as a
%   function of the degree of the polynomial h(x) in Valmorbida et al. [2] (Sec. V-A); see also [3].
%   Free variable y in L_2[0, 1]
%   step 1. Introduce v such that: x = int_0^1 v (implies that x_s = v).
%   step 2. Introduce y such that: v = y - int_0^1 y 
%                           (implies that int v =int x_s = x(1) - x(0) = 0)%  
%  x =  T (I - Pi)y , where 
% (Ty)(s) = int_0^s y(theta) dtheta -- Voltera Type operator
% (Piy)(s) = int_0^1 y(theta) dtheta -- Averaging operator
%   LPI:  min C  s.t.  C*I - M'*M >= 0
        name = 'Wirtinger: int x^2 <= C int x_s^2, x(0)=x(b)=0';
        dg = degrees(opts,2,2);

        % Operators acting on the free variable y in L2[0,b]:
        % (Ty)(s) = int_0^s y(theta) dtheta -- Volterra-type operator (integrates y from 0 to s)
        T  = get_op(int(y(s1_dum),s1_dum,[0,s1]));
        % (Pi y)(s) = (1/b) int_0^b y(theta) dtheta -- averaging operator (maps y to its mean, a constant function)
        Pi = get_op((1/b)*int(y,s1,[0,s1,b]));
        % identity operator on L2[0,b]
        I1 = get_op(y);

        % M = T*(I - Pi) maps y to x:
        % (My)(s) = int_0^s y(theta) dtheta - (s/b) int_0^b y(theta) dtheta,
        % so x(0) = x(b) = 0 and x_s = y - mean(y). M kills constants (M*Pi = 0), so the worst case
        % of ||My||^2 <= C ||y||^2 has y with zero mean, where ||y|| = ||x_s||.
        M =  T*(I1 - Pi);

        % optimization problem min C s.t. C*I - M'*M >= 0
        Csol = solve_min(dom,@(C) C - M'*M,dg);

        % sharp value
        sharp = (b/pi)^2;
%--------------------------------------------------------------------------
    case 2
%   OPTIMIZATION
%   min C
%   s.t. ||x||^2 <= C ||x_s||^2,                       |        Wirtinger (periodic)
%        x(0) = x(1) (periodic BC),                    |        Sharp C = (b/(2*pi))^2
%        int_0^1 x = 0
%   References: Valmorbida and Papachristodoulou [3] (Sec. II-A) treat this inequality with INTSOSTOOLS:
%   their Lemma 1 shows that the worst case is x(0) = x(1) = 0 (replace x by x - k), which is why the
%   common value k = x(0) = x(b) is free here, and v = int_0^s x then gives the second-order form
%   int v_s^2 <= C int v_ss^2, v = v_s = 0 at both ends (clamped beam, example 15); see also [2].
%   Free variable y in L_2[0, 1]
%   step 1. Introduce v such that: x = v - int_0^1 v (implies that int_0^1 x = 0).
%   step 2. Introduce w such that: v = int_0^s w (implies that v_s=x_s=w)
%   step 3. Introduce y such that: w = y - int_0^1 y 
%                           (implies that int w =int x_s = x(1) - x(0) = 0)
%  
%  x = (I - Pi) T (I - Pi)y , where 
% (Ty)(s) = int_0^s y(theta) dtheta -- Voltera Type operator
% (Piy)(s) = int_0^1 y(theta) dtheta -- Averaging operator
%   LPI:  min C  s.t.  C*I - M'*M >= 0
        name = 'Wirtinger (periodic): int x^2 <= C int x_s^2, x(0)=x(b), int x = 0';
        dg = degrees(opts,2,2);

        % Operators acting on the free variable y in L2[0,b]:
        % (Ty)(s) = int_0^s y(theta) dtheta -- Volterra-type operator (integrates y from 0 to s)
        T  = get_op(int(y(s1_dum),s1_dum,[0,s1]));
        % (Pi y)(s) = (1/b) int_0^b y(theta) dtheta -- averaging operator (maps y to its mean, a constant function)
        Pi = get_op((1/b)*int(y,s1,[0,s1,b]));
        % identity operator on L2[0,b]
        I1 = get_op(y);

        % M = (I - Pi)*T*(I - Pi) maps y to x in three steps:
        %   w = (I - Pi) y   w has zero mean (w = x_s),
        %   v = T w          v(0) = 0 and v(b) = int w = 0,
        %   x = (I - Pi) v   x has zero mean and x(0) = x(b) = -mean(v) (periodic), with x_s = w.
        M  = (I1-Pi)*T*(I1-Pi);

        % optimization problem min C s.t. C*I - M'*M >= 0
        Csol = solve_min(dom,@(C) C - M'*M,dg);

        % sharp value
        sharp = (b/(2*pi))^2;
%--------------------------------------------------------------------------
    case 3
%   min C
%   s.t. ||x||^2 <= C ||x_s||^2,                             | b = 1        Poincare (Poincare-Wirtinger)
%        int_0^b x = 0                                     |                Sharp C = (b/pi)^2
%   Free variable y in L_2[0, b]
%   step 1. Introduce v such that: x = v - (1/b) int_0^b v
%                                           (implies that int_0^b x = 0).
%   step 2. Introduce y such that: v = int_0^s y   (v(0) = int_0^b v= 0)
%                                           (implies that x_s = v_s = y).
%   x = (I - Pi) T y, where
%   (Ty)(s) = int_0^s y(theta) dtheta -- Volterra-type operator
%   (Pi y)(s) = (1/b) int_0^b y(theta) dtheta -- Averaging operator
%   LPI:  min C  s.t.  C*I - M'*M >= 0,   M = (I - Pi) T
        name = 'Poincare: int x^2 <= C int x_s^2, int x = 0';
        dg = degrees(opts,2,2);

        % Operators acting on the free variable y in L2[0,b]:
        % (Ty)(s) = int_0^s y(theta) dtheta -- Volterra-type operator (integrates y from 0 to s)
        T  = get_op(int(y(s1_dum),s1_dum,[0,s1]));
        % (Pi y)(s) = (1/b) int_0^b y(theta) dtheta -- averaging operator (maps y to its mean, a constant function)
        Pi = get_op((1/b)*int(y,s1,[0,s1,b]));
        % identity operator on L2[0,b]
        I1 = get_op(y);

        % M = (I - Pi)*T maps y to x: first v = T y (v(0) = 0, v_s = y), then x = v - mean(v):
        % (My)(s) = int_0^s y(theta) dtheta - (1/b) int_0^b (b-theta) y(theta) dtheta,
        % so int_0^b My = 0 and (My)_s = y for any y in L2 (no constraint on y)
        M  = (I1-Pi)*T;

        % optimization problem min C s.t. C*I - M'*M >= 0
        Csol = solve_min(dom,@(C) C - M'*M,dg);

        % sharp value
        sharp = (b/pi)^2;
%--------------------------------------------------------------------------
    case 4
%   min C
%   s.t. ||x||^2 <= C ||x_s||^2,                             | b = 1        Poincare, one zero end
%        x(0) = 0                                            |                Sharp C = (2*b/pi)^2
%   Free variable y in L_2[0, b]
%   Introduce y such that: x = int_0^s y
%                           (implies that x(0) = 0 and x_s = y).
%   x = T y, where
%   (Ty)(s) = int_0^s y(theta) dtheta -- Volterra-type operator
%   LPI:  min C  s.t.  C*I - T'*T >= 0
        name = 'Poincare: int x^2 <= C int x_s^2, x(0)=0';
        dg = degrees(opts,2,2);

        % Operator acting on the free variable y in L2[0,b]:
        % (Ty)(s) = int_0^s y(theta) dtheta -- Volterra-type operator; x = T y has x(0) = 0 and x_s = y
        T = get_op(int(y(s1_dum),s1_dum,[0,s1]));

        % optimization problem min C s.t. C*I - T'*T >= 0, i.e. ||x||^2 = ||Ty||^2 <= C ||y||^2 = C ||x_s||^2
        Csol = solve_min(dom,@(C) C - T'*T,dg);

        % sharp value
        sharp = (2*b/pi)^2;
%--------------------------------------------------------------------------
    case 5
%   min C
%   s.t. ||x||^2 <= C (||x_s||^2 + x(0)^2)                   | b = 1        Barraza Martinez et al. [1] (3.2)
%        x in H1[0,b] (no boundary conditions)             |                Paper: C = 2*max(b,b^2)
%   Free variable (x0, y) in R x L_2[0, b]
%   Introduce x such that: x = x0 + int_0^s y
%                           (implies that x(0) = x0 and x_s = y).
%   Define operators T, E0, B such that
%    x = T(x0,y), x(0) = E0(x0,y), x_s = B(x0,y), where
%   T(x0,y)(s) = x0 + int_0^s y(theta) dtheta,  E0(x0,y) = x0,  B(x0,y) = y
%   LPI:  min C  s.t.  C*(B'*B + E0'*E0) - T'*T >= 0
%   For b = 1 the optimal value is 1.351 (finite elements).
        name = 'Generalized Poincare: int x^2 <= C (int x_s^2 + x(0)^2)';
        dg = degrees(opts,2,2);

        % Operators acting on the state (x0,y) in R x L2[0,b], with x0 = x(0) and y = x_s:
        % T(x0,y)(s) := x0 + int_0^s y(theta) dtheta = x(s) -- maps the state to the function x (L2 output)
        T  = get_op(x0 + int(y(s1_dum),s1_dum,[0,s1]));
        % E0(x0,y) := x0 = x(0) -- maps the state to the point value x(0) (an element of R)
        E0 = get_op(x0 + int(0*y,s1,[0,b]));
        % B(x0,y) := (0,y) = (0,x_s) -- the projection (x0,y) -> y in R x L2, so <B.,B.> = ||x_s||^2 and B'*B = B
        % (the code below uses B itself as the quadratic form B'*B)
        B  = get_op([0*x0; y]);

        % optimization problem min C s.t. C*(B'*B + E0'*E0) - T'*T >= 0, i.e. ||x||^2 <= C (||x_s||^2 + x(0)^2)
        Csol = solve_min(dom,@(C) C*(B + E0'*E0) - T'*T,dg);

        % reference value of the paper
        paper = 2*max(b,b^2);
%==========================================================================
%       Trace inequalities
%==========================================================================
    case 6
%   min C
%   s.t. x(0)^2 + x(b)^2 <= C (||x_s||^2 + ||x||^2),                | b = 1   Barraza Martinez et al. [1] (2.5)
%        x in H1[0,b] (no boundary conditions)                      |           Paper: C = max(4*b,4/b)
%   Free variable (x0, y) in R x L_2[0, b]
%   Introduce x0,y such that: x = x0 + int_0^s y
%             (x(0) = x0, x(b) = x0 + int_0^b y and x_s = y).
%   Define operators T, E0, E1, B such that
%   x = T(x0,y), x(0) = E0(x0,y), x(b) = E1(x0,y), x_s = B (x0,y), where
%   T(x0,y)(s) = x0 + int_0^s y(theta) dtheta,  
%   E0(x0,y) = x0,
%   E1(x0,y) = x0 + int_0^b y(theta) dtheta,
%   B(x0,y) = y
%   LPI:  min C  s.t.  C*(X'*X + B'*B) - E0'*E0 - E1'*E1 >= 0
%   For b = 1 the optimal value is 2.164 (finite elements).
        name = 'Trace: x(0)^2 + x(b)^2 <= C (int x_s^2 + int x^2)';
        dg = degrees(opts,2,2);

        % Operators acting on the state (x0,y) in R x L2[0,b], with x0 = x(0) and y = x_s:
        % T(x0,y)(s) := x0 + int_0^s y(theta) dtheta = x(s) -- maps the state to the function x (L2 output)
        T  = get_op(x0 + int(y(s1_dum),s1_dum,[0,s1]));
        % E0(x0,y) := x0 = x(0) -- maps the state to the point value x(0) (an element of R)
        E0 = get_op(x0 + int(0*y,s1,[0,b]));
        % E1(x0,y) := x0 + int_0^b y(theta) dtheta = x(b) -- point value at the right end (an element of R)
        E1 = get_op(x0 + int(y,s1,[0,b]));
        % B(x0,y) := (0,y) = (0,x_s) -- the projection (x0,y) -> y in R x L2, so <B.,B.> = ||x_s||^2 and B'*B = B
        % (the code below uses B itself as the quadratic form B'*B)
        B  = get_op([0*x0; y]);

        % optimization problem min C s.t. C*(T'*T + B'*B) - E0'*E0 - E1'*E1 >= 0,
        % i.e. x(0)^2 + x(b)^2 <= C (||x||^2 + ||x_s||^2)
        Csol = solve_min(dom,@(C) C*(T'*T + B) - E0'*E0 - E1'*E1,dg);

        % reference value of the paper
        paper = max(4*b,4/b);
%--------------------------------------------------------------------------
    case 7
%   min C
%   s.t. x(b)^2 <= C (||x||^2 + ||x_s||^2),                | b = 1        Sobolev-type (single endpoint)
%        x in H1[0,b] (no boundary conditions)             |                Sharp C = coth(b)  (weight of ||x_s||^2 equal to C gives k = 1, C tanh(b) = 1)
%   Free variable (x0, y) in R x L_2[0, b]
%   Introduce x0,y such that: x = x0 + int_0^s y
%                           (implies that x(b) = x0 + int_0^b y and x_s = y).
%   Define operators T, E0, B such that
%   x = T(x0,y), x(b) = E(x0,y), x_s = B(x0,y), where
%   T(x0,y)(s) = x0 + int_0^s y(theta) dtheta,
%   E(x0,y) = x0 + int_0^b y(theta) dtheta  
%   B(x0,y) = y
%   LPI:  min C  s.t.  C*(T'*T + B'*B) - E'*E >= 0
        name = 'Sobolev trace: x(b)^2 <= C (||x||^2 + ||x_s||^2)';
        dg = degrees(opts,[],[]);

        % Operators acting on the state (x0,y) in R x L2[0,b], with x0 = x(0) and y = x_s:
        % X(x0,y)(s) := x0 + int_0^s y(theta) dtheta = x(s) -- maps the state to the function x (L2 output)
        X = get_op(x0 + int(y(s1_dum),s1_dum,[0,s1]));
        % E(x0,y) := x0 + int_0^b y(theta) dtheta = x(b) -- maps the state to the point value x(b) (an element of R)
        E = get_op(x0 + int(y,s1,[0,b]));
        % B(x0,y) := (0,y) = (0,x_s) -- the projection (x0,y) -> y in R x L2, so <B.,B.> = ||x_s||^2 and B'*B = B
        % (the code below uses B itself as the quadratic form B'*B)
        B = get_op([0*x0; y]);

        % optimization problem min C s.t. C*(X'*X + B'*B) - E'*E >= 0,
        % i.e. x(b)^2 <= C (||x||^2 + ||x_s||^2)
        Csol = solve_min(dom,@(C) C*(X'*X + B) - E'*E,dg);

        % sharp value (Euler-Lagrange: weight = C gives k = 1, C tanh(b) = 1)
        sharp = 1/tanh(b);
%--------------------------------------------------------------------------
    case 8
%   min C
%   s.t. x(t)^2 <= C (||x||^2 + ||x_s||^2),                  | b = 1        Agmon / 1D Sobolev (interior point t)
%        x in H1[0,b] (no boundary conditions)             | t = b/2
%                                                            |                Sharp C = 1/(tanh(t)+tanh(b-t))  (weight = C: k = 1)
%   Free variable (x*, y1, y2) in R x L_2[0, 1] x L_2[0, 1]
%   Split [0,b] at t and rescale both pieces to [0,1]:
%   u1(tau) = x(t*tau),   u2(tau) = x(t + (b-t)*tau).
%   Introduce x*, y1, y2 such that: u1(tau) = x* - int_tau^1 y1,  
%                                   u2(tau)= x* + int_0^tau y2
%             (x* = x(t) is shared by both pieces, so x is continuous at t;
%              y1 = t*x_s on [0,t] and y2 = (b-t)*x_s on [t,b]).
%   Define operators U1, U2, B1, B2, E such that
%   u1 = U1 z, u2 = U2 z, y1 = B1 z, y2 = B2 z, x(t) = E z, for z = (x*,y1,y2), where
%   (U1 z)(tau) = x* - int_tau^1 y1(theta) dtheta,
%   (U2 z)(tau) = x* + int_0^tau y2(theta) dtheta,
%   B1 z = y1,  B2 z = y2,  E z = x*
%   Then ||x||^2 = t ||u1||^2 + (b-t) ||u2||^2,
%        ||x_s||^2 = ||y1||^2/t + ||y2||^2/(b-t).
%   LPI:  min C  
%   s.t.  C*(t*U1'*U1 + (b-t)*U2'*U2 + (1/t)*B1'*B1 + (1/(b-t))*B2'*B2) - E'*E >= 0
%   The sup over t of x(t)^2 is attained at t = 0 or b, i.e. by example 7:
%   the constant here is smaller for interior t.
        dg = degrees(opts,[],[]);
        t = opts.t;
        if isempty(t),  t = b/2;  end
        if t <= 0 || t >= b
            error('Option t must satisfy 0 < t < b (t = 0 or b is example 7).');
        end
        name = sprintf('Agmon / 1D Sobolev: x(t)^2 <= C (||x||^2 + ||x_s||^2), t = %g in (0,%g)',t,b);

        % After the split at t both pieces are rescaled to [0,1], so the operators below
        % act on the state z = (x*, y1, y2) = (x0, y12) in R x L2^2[0,1] and live on dom1 = [0,1].
        % The row [a c] multiplies (y1,y2) by (a,c).
        dom1 = [0,1];

        % Operators from the state z to the two pieces of x (L2[0,1] outputs):
        % (U1 z)(tau) := x* - int_tau^1 y1(theta) dtheta  -- left piece u1(tau) = x(t*tau)
        U1 = get_op(x0 - int([1 0]*y12(s1_dum),s1_dum,[s1,1]));
        % (U2 z)(tau) := x* + int_0^tau y2(theta) dtheta  -- right piece u2(tau) = x(t + (b-t)*tau)
        U2 = get_op(x0 + int([0 1]*y12(s1_dum),s1_dum,[0,s1]));

        % Projections of the state onto its derivative components (L2[0,1] outputs):
        % B1 z := y1 = t*x_s on [0,t],   B2 z := y2 = (b-t)*x_s on [t,b]
        B1 = get_op(0*x0 + [1 0]*y12);
        B2 = get_op(0*x0 + [0 1]*y12);

        % Point evaluation at the split point (an element of R):
        % E z := x* = x(t)
        E  = get_op(x0 + int([0 0]*y12,s1,[0,1]));

        % The norms of x and x_s as quadratic forms of z (with the Jacobian factors of the split):
        % <z,Qx z> = ||x||^2 = t ||u1||^2 + (b-t) ||u2||^2
        Qx = t*(U1'*U1) + (b-t)*(U2'*U2);
        % <z,Qs z> = ||x_s||^2 = ||y1||^2/t + ||y2||^2/(b-t)
        Qs = (1/t)*(B1'*B1) + (1/(b-t))*(B2'*B2);

        % optimization problem min C s.t. C*(Qx + Qs) - E'*E >= 0, i.e. x(t)^2 <= C (||x||^2 + ||x_s||^2)
        Csol = solve_min(dom1,@(C) C*(Qx + Qs) - E'*E,dg);

        % sharp value (Euler-Lagrange with the jump condition at t)
        sharp = 1/(tanh(t) + tanh(b-t));
        fprintf('\n t = %g (interval [0,%g]): sharp C(t) = %.6f; the endpoint value (example 7) is %.6f.\n', ...
            t,b,sharp,1/tanh(b));
%--------------------------------------------------------------------------
    case 9
%   Agmon / 1D Sobolev, supremum over t (t is a spatial variable)
%   min C
%   s.t. sup_{t in [tmin, b-tmin] x(t)^2 <= C (||x||^2 + ||x_s||^2)   
%                                                          | b = 1  
%      x in H1[0,b] (no boundary conditions)               | Sharp C = C(tmin) = 1/(tanh(tmin)+tanh(b-tmin)) 
%                                                                               
%   Lift t to a second spatial variable and split [0,b] at t as in example 8:
%   Free variable (x*, y1, y2) in L_2[t] x L_2[tau, t] x L_2[tau, t], with t in [tmin, b-tmin], tau in [0,1]
%   Introduce x*, y1, y2 such that: u1(tau,t) = x*(t) - int_tau^1 y1(theta,t) dtheta,
%                                   u2(tau,t) = x*(t) + int_0^tau y2(theta,t) dtheta
%             (u1(tau,t) = x(t*tau), u2(tau,t) = x(t + (b-t)*tau), x*(t) = x(t),
%              y1 = t*x_s on [0,t] and y2 = (b-t)*x_s on [t,b], for every t).
%   Define operators U1, U2, B1, B2, E (as in example 8, with t a variable) such that
%   u1 = U1 z, u2 = U2 z, y1 = B1 z, y2 = B2 z, x(t) = E z, for z = (x*,y1,y2).
%   Then, for each t, ||x||^2 = t ||u1||^2 + (b-t) ||u2||^2,
%        ||x_s||^2 = ||y1||^2/t + ||y2||^2/(b-t).
%   Multiplying the operator inequality of example 8 by t*(b-t) > 0 makes every
%   coefficient a polynomial in t. The operator acts in tau and only MULTIPLIES in t, so
%   it is a multiplier in the spatial variable t and positivity on L_2[tau,t] is
%   equivalent to positivity for every t in [tmin, b-tmin].
%   LPI:  min C  s.t.  C*( Mw1(t)*U1'*U1 + Mw2(t)*U2'*U2
%                        + Mc(t)*B1'*B1 + Md(t)*B2'*B2 ) - Me(t)*E'*E >= 0,
%         Mw1 = t^2 (b-t), Mw2 = t (b-t)^2, Mc = b-t, Md = t, Me = t (b-t)
%         (multiplication operators in t).
%   The supremum over t of C(t) is attained at the end of the t-interval, so the
%   result is C(tmin). With tmin = 0 the constraint is void at t = 0 and t = b
%   (weight t*(b-t)) and the LPI is degenerate: use example 7 for that case.
        dg = degrees(opts,2,2);
        tmin = opts.tmin;
        if tmin < 0 || tmin >= b/2
            error('Option tmin must satisfy 0 <= tmin < b/2.');
        end
        name = sprintf('Agmon / 1D Sobolev, sup over t in [%g,%g]: x(t)^2 <= C (||x||^2 + ||x_s||^2)',tmin,b-tmin);
        if tmin == 0
            disp('Warning: tmin = 0 makes the lifted LPI degenerate at t = 0 and t = b (use example 7 instead).');
        end

        % The point t is now a second spatial variable, called s2 below, with domain [tmin,b-tmin];
        % tau (called s1 below) lives on [0,1].

        % Lifted state z = (x*, y1, y2): x*(s2) in L2[s2], and y1(tau,s2), y2(tau,s2) in L2[tau,s2]
        % (declared here, since their domain depends on tmin).
        xs = pde_var(s2,[tmin,b-tmin]);
        y1 = pde_var([s1;s2],[0,1; tmin,b-tmin]);
        y2 = pde_var([s1;s2],[0,1; tmin,b-tmin]);

        % Multiplication operators by a function w(s2), acting in s2 only:
        % (mulY(w) z)(tau,s2) = w(s2) z(tau,s2) on L2[tau,s2],   (mulX(w) z)(s2) = w(s2) z(s2) on L2[s2]
        mulY = @(w) get_op(w*y1);
        mulX = @(w) get_op(w*xs);

        % Operators from the state z = (x*,y1,y2) to the two pieces of x, as in example 8 but with t = s2 a variable
        % (each is a 1 x 3 block row acting on (x*,y1,y2); the zero terms only fix the spaces):
        % (U1 z)(tau,s2) := x*(s2) - int_tau^1 y1(theta,s2) dtheta   -- left piece  u1(tau,s2) = x(s2*tau)
        U1 = get_op(xs - int(y1(s1_dum,s2),s1_dum,[s1,1]) + 0*y2);
        % (U2 z)(tau,s2) := x*(s2) + int_0^tau y2(theta,s2) dtheta   -- right piece u2(tau,s2) = x(s2 + (b-s2)*tau)
        U2 = get_op(xs + 0*y1 + int(y2(s1_dum,s2),s1_dum,[0,s1]));

        % Projections onto the derivative components (L2[tau,s2] outputs):
        % B1 z := y1 = s2*x_s on [0,s2],   B2 z := y2 = (b-s2)*x_s on [s2,b]
        B1 = get_op(0*xs + y1 + 0*y2);
        B2 = get_op(0*xs + 0*y1 + y2);

        % Point evaluation x(s2) = x*(s2), an L2[s2] output: a multiplier in s2, so no trace on the diagonal is needed
        % E z := x*(s2) = x(s2)
        E  = get_op(xs + int(0*y1,s1,[0,1]) + int(0*y2,s1,[0,1]));

        % The norms as quadratic forms of z, multiplied by s2*(b-s2) > 0 so that all coefficients are polynomials in s2:
        % <z,Qx z> = s2(b-s2) ||x||^2 = s2(b-s2) ( s2 ||u1||^2 + (b-s2) ||u2||^2 )
        Qx = U1'*mulY(s2^2*(b-s2))*U1 + U2'*mulY(s2*(b-s2)^2)*U2;
        % <z,Qs z> = s2(b-s2) ||x_s||^2 = s2(b-s2) ( ||y1||^2/s2 + ||y2||^2/(b-s2) )
        Qs = B1'*mulY(b-s2)*B1 + B2'*mulY(s2)*B2;
        % <z,QE z> = s2(b-s2) x(s2)^2
        QE = E'*mulX(s2*(b-s2))*E;

        % optimization problem min C s.t. C*(Qx + Qs) - QE >= 0, for all s2 in [tmin,b-tmin] (variables s1 and s2)
        prog = lpiprogram([s1;s2],[],[0 1; tmin b-tmin]);
        [prog,C] = lpidecvar(prog,'C');
        if isempty(dg)
            prog = lpi_ineq(prog,C*(Qx + Qs) - QE);
        else
            prog = lpi_ineq(prog,C*(Qx + Qs) - QE,dg);
        end
        prog = lpisetobj(prog,C);
        prog = lpisolve(prog);
        Csol = double(lpigetsol(prog,C));

        % sharp value: C(t) is largest at the end of the t-interval, so the result is C(tmin)
        sharp = 1/(tanh(tmin) + tanh(b-tmin));
        fprintf('\n tmin = %g, t in [%g,%g]: sharp C(tmin) = %.6f; the endpoint value (example 7) is %.6f.\n', ...
            tmin,tmin,b-tmin,sharp,1/tanh(b));
%==========================================================================
%       Hardy-type and weighted inequalities
%==========================================================================
    case 10
%   min C
%   s.t. int_0^1 x^2 <= C int_0^1 (s^2+eps) x_{s}^2,         | b = 1      Hardy, alpha = 2
%        x(1) = 0                                            | eps = 1e-2   C(eps) -> 4 as eps -> 0 (sharp, not attained)
%   Free variable y in L_2[0, 1]
%   Introduce y such that: x = -int_s^1 y
%                           (implies that x(1) = 0 and x_s = y).
%   x = T y, where
%   (Ty)(s) = -int_s^1 y(theta) dtheta -- Volterra-type operator
%   (Sy)(s) = s*y(s) -- multiplication operator, S'*S is the weight s^2
%   LPI:  min C  s.t.  C*(S'*S + eps*I) - T'*T >= 0
%   The value is compared with 1/lambda_min of -((s^2+eps)x')' = lambda x,
%   x'(0) = 0, x(1) = 0 (finite differences). For eps = 0 the SDP has no
%   strictly feasible point and the solver returns values below 4, so eps > 0.
        name = 'Hardy: int x^2 <= C int (s^2+eps) x_s^2, x(1)=0';
        dom = [0,1];    b = 1;
        dg = degrees(opts,3,3);
        eps_h = opts.eps;

        % Operators acting on the free variable y in L2[0,1]:
        % (Ty)(s) := -int_s^1 y(theta) dtheta -- Volterra-type operator (integrates y from s to 1).
        % T(1) = 0 and (Ty)_s = y, so x = Ty for any y = x_s with x(1) = 0
        T = get_op(-int(y01(s1_dum),s1_dum,[s1,1]));
        % (Sy)(s) := s*y(s) -- multiplication operator by s; S'*S is the multiplication by the weight s^2
        S = get_op(s1*y01);
        % identity operator on L2[0,1]
        I1 = get_op(y01);

        % optimization problem min C s.t. C*(S'*S + eps*I) - T'*T >= 0,
        % i.e. ||x||^2 = ||Ty||^2 <= C int (s^2+eps) y^2 = C int (s^2+eps) x_s^2
        Csol = solve_min(dom,@(C) C*(S'*S + eps_h*I1) - T'*T,dg);

        % finite-difference reference (cell centered, x'(0) = 0, x(1) = 0)
        N = 4000; h = 1/N; sf = (0:N)'*h;
        w = sf.^2 + eps_h;
        main = (w(1:N)+w(2:N+1))/h^2;  main(1) = w(2)/h^2;  main(N) = (w(N)+2*w(N+1))/h^2;
        K = spdiags([-w(2:N)/h^2; 0],-1,N,N) + spdiags([0; -w(2:N)/h^2],1,N,N) + spdiags(main,0,N,N);
        sharp = 1/eigs(K,1,'smallestabs');
%==========================================================================
%       Operator norms
%==========================================================================
    case 11
%   min C
%   s.t. int_0^b (Tx)^2 <= C int_0^b x^2,                | b = 1        Volterra operator norm
%        (Tx)(s) = int_0^s x(r) dr, x in L2[0,b]           |                Sharp ||T|| = 2*b/pi, so C = (2*b/pi)^2
%   Free variable x in L_2[0, b]
%   T x = int_0^s x(theta) dtheta, where
%   (Tx)(s) = int_0^s x(theta) dtheta -- Volterra-type operator
%   LPI:  min C  s.t.  C*I - T'*T >= 0   (the norm is sqrt(C))
%   It is the same operator inequality as example 4 (Poincare with
%   x(0) = 0, x = T*x_s), read as a bound on the norm of T.
        name = 'Volterra operator norm: ||T x||^2 <= C ||x||^2, (T x)(s) = int_0^s x';
        dg = degrees(opts,[],[]);

        % Operator on L2[0,b]:
        % (Tx)(s) := int_0^s x(r) dr -- Volterra-type operator (integrates x from 0 to s)
        T = get_op(int(y(s1_dum),s1_dum,[0,s1]));

        % optimization problem min C s.t. C*I - T'*T >= 0, i.e. ||Tx||^2 <= C ||x||^2 for all x in L2
        Csol = solve_min(dom,@(C) C - T'*T,dg);

        % sharp value (||T|| = 2*b/pi)
        sharp = (2*b/pi)^2;
        fprintf('\n Operator norm ||T|| <= %.6f (exact 2*b/pi = %.6f).\n',sqrt(Csol),2*b/pi);
%--------------------------------------------------------------------------
    case 12
%   min C
%   s.t. int_0^b (Gx)^2 <= C int_0^b x^2,                | b = 1        Green's operator, Dirichlet (-x'' = f, x(0) = x(b) = 0)
%        (Gx)(s) = int_0^b G(s,theta) x(theta) dtheta,     |                Sharp ||G|| = (b/pi)^2, so C = (b/pi)^4
%        x in L2[0,b]
%   Free variable x in L_2[0, b]
%   G(s,theta) = theta*(b-s)/b (theta < s),  s*(b-theta)/b (theta > s), where
%   (Gx)(s) = int_0^s theta*(b-s)/b x(theta) dtheta + int_s^b s*(b-theta)/b x(theta) dtheta
%   is the solution of -u'' = x, u(0) = u(b) = 0, and ||G|| = 1/lambda_1 with lambda_1 = (pi/b)^2.
%   LPI:  min C  s.t.  C*I - G'*G >= 0   (the norm is sqrt(C))
        name = 'Green operator (Dirichlet): ||G x||^2 <= C ||x||^2';
        dg = degrees(opts,3,3);

        % Operator on L2[0,b], with s1 = s and s1_dum = theta:
        % (Gx)(s) := int_0^s theta*(b-s)/b x(theta) dtheta + int_s^b s*(b-theta)/b x(theta) dtheta
        % -- Green's operator of -u'' = x with u(0) = u(b) = 0, i.e. (Gx)'' = -x
        G = get_op(int(s1_dum*(b-s1)/b*y(s1_dum),s1_dum,[0,s1]) + int(s1*(b-s1_dum)/b*y(s1_dum),s1_dum,[s1,b]));

        % optimization problem min C s.t. C*I - G'*G >= 0, i.e. ||Gx||^2 <= C ||x||^2 for all x in L2
        Csol = solve_min(dom,@(C) C - G'*G,dg);

        % sharp value (||G|| = 1/lambda_1 = (b/pi)^2)
        sharp = (b/pi)^4;
        fprintf('\n Operator norm ||G|| <= %.6f (exact (b/pi)^2 = %.6f).\n',sqrt(Csol),(b/pi)^2);
%--------------------------------------------------------------------------
    case 13
%   min C
%   s.t. int_0^b (Gx)^2 <= C int_0^b x^2,                | b = 1        Green's operator, mixed (-x'' = f, x(0) = 0, x_s(b) = 0)
%        (Gx)(s) = int_0^b min(s,theta) x(theta) dtheta,   |                Sharp ||G|| = 4*b^2/pi^2, so C = (4*b^2/pi^2)^2
%        x in L2[0,b]
%   Free variable x in L_2[0, b]
%   G(s,theta) = min(s,theta), where
%   (Gx)(s) = int_0^s theta x(theta) dtheta + int_s^b s x(theta) dtheta
%   is the solution of -u'' = x, u(0) = 0, u_s(b) = 0, and ||G|| = 1/lambda_1 with lambda_1 = (pi/(2*b))^2.
%   LPI:  min C  s.t.  C*I - G'*G >= 0   (the norm is sqrt(C))
        name = 'Green operator (mixed): ||G x||^2 <= C ||x||^2';
        dg = degrees(opts,3,3);

        % Operator on L2[0,b] with kernel min(s,theta), with s1 = s and s1_dum = theta:
        % (Gx)(s) := int_0^s theta x(theta) dtheta + int_s^b s x(theta) dtheta
        % -- Green's operator of -u'' = x with u(0) = 0 and u_s(b) = 0, i.e. (Gx)'' = -x
        G = get_op(int(s1_dum*y(s1_dum),s1_dum,[0,s1]) + int(s1*y(s1_dum),s1_dum,[s1,b]));

        % optimization problem min C s.t. C*I - G'*G >= 0, i.e. ||Gx||^2 <= C ||x||^2 for all x in L2
        Csol = solve_min(dom,@(C) C - G'*G,dg);

        % sharp value (||G|| = 1/lambda_1 = 4*b^2/pi^2)
        sharp = (4*b^2/pi^2)^2;
        fprintf('\n Operator norm ||G|| <= %.6f (exact 4*b^2/pi^2 = %.6f).\n',sqrt(Csol),4*b^2/pi^2);
%--------------------------------------------------------------------------
    case 14
%   min C
%   s.t. int_0^b (Bx)^2 <= C int_0^b x^2,                | b = 1        Green's operator, hinged beam (u'''' = f, u = u'' = 0 at 0 and b)
%        B = G*G, G = Green's operator of example 12        |                Sharp ||B|| = (b/pi)^4, so C = (b/pi)^8
%        x in L2[0,b]
%   Free variable x in L_2[0, b]
%   The solution of u'''' = x, u = u'' = 0 at s = 0 and s = b, is obtained by solving
%   -v'' = x, v(0) = v(b) = 0 and then -u'' = v, u(0) = u(b) = 0, so B = G*G
%   with G the Dirichlet Green's operator (a kernel of degree 3), and ||B|| = 1/lambda_1^2 with lambda_1 = (pi/b)^2.
%   LPI:  min C  s.t.  C*I - B'*B >= 0   (the norm is sqrt(C))
        name = 'Green operator (hinged beam): ||B x||^2 <= C ||x||^2';
        dg = degrees(opts,4,4);

        % Dirichlet Green's operator as in example 12 (s1 = s and s1_dum = theta):
        % (Gx)(s) := int_0^s theta*(b-s)/b x(theta) dtheta + int_s^b s*(b-theta)/b x(theta) dtheta,  (Gx)'' = -x
        G = get_op(int(s1_dum*(b-s1)/b*y(s1_dum),s1_dum,[0,s1]) + int(s1*(b-s1_dum)/b*y(s1_dum),s1_dum,[s1,b]));

        % B := G*G (composition of two Dirichlet solves: x -> v = Gx with -v'' = x, then u = Gv with -u'' = v),
        % so Bx = u solves u'''' = x with u = u'' = 0 at both ends (hinged beam)
        B = G*G;

        % optimization problem min C s.t. C*I - B'*B >= 0, i.e. ||Bx||^2 <= C ||x||^2 for all x in L2
        Csol = solve_min(dom,@(C) C - B'*B,dg);

        % sharp value (||B|| = 1/lambda_1^2 = (b/pi)^4)
        sharp = (b/pi)^8;
        fprintf('\n Operator norm ||B|| <= %.6g (exact (b/pi)^4 = %.6g).\n',sqrt(Csol),(b/pi)^4);
%==========================================================================
%       Eigenvalue problems
%==========================================================================
    case 15
%   min C
%   s.t. int_0^b u_{s}^2 <= C int_0^b u_{ss}^2,          | b = 1        Clamped-clamped problem  u_{ss} = lambda*u_{ssss}
%        u(0) = u(b) = u_s(0) = u_s(b) = 0               |                Sharp C = (b/(2*pi))^2, so lambda_max = -C
%   For u_{ss} = lambda*u_{ssss} with these boundary conditions, multiplying by u and
%   integrating by parts gives -||u_s||^2 = lambda*||u_ss||^2, so lambda = -||u_s||^2/||u_ss||^2.
%   The eigenvalue of largest magnitude is lambda = -C, with C the sharp constant above
%   (the buckling load of the clamped column, P = 1/C = (2*pi/b)^2).
%   Free variable v = u_{ssss} in L_2[0, b]
%   Introduce u such that: u_{ssss} = v, u(0) = u(b) = u_s(0) = u_s(b) = 0
%             (the four boundary conditions determine u from v).
%   Define operators B and P such that, via the PIE of the PDE
%   u_t = u_ss + u_ssss with the same boundary conditions, where
%   B v = u (the map PIE.T) and the PIE operator A = P + I, since
%   (u_ss + u_ssss) = P v + v with v = u_ssss; so P v = u_ss and
%   P = PIE.A - I. Also u_s = (d/ds B) v.
%   LPI:  min C  s.t.  C*P'*P - (d/ds B)'*(d/ds B) >= 0
        name = 'Clamped beam eigenvalue problem u_ss = lambda u_ssss: ||u_s||^2 <= C ||u_ss||^2';
        dg = degrees(opts,3,3);

        % PDE u_t = u_ss + u_ssss with the four clamped boundary conditions
        % (u(0) = u(b) = u_s(0) = u_s(b) = 0), for u = u4 in H4[0,b]; its PIE gives the operators below.
        PDE = [diff(u4,t) == diff(u4,s1,2) + diff(u4,s1,4);
               subs(u4,s1,0) == 0;             subs(u4,s1,b) == 0;
               subs(diff(u4,s1,1),s1,0) == 0;  subs(diff(u4,s1,1),s1,b) == 0];

        % PIE: T v_t = A v for the fundamental state v = u_ssss in L2[0,b], with u = T v
        PIE = convert(PDE,'pie');

        % Operators acting on the fundamental state v = u_ssss in L2[0,b]:
        % B v = u    -- B = PIE.T, the operator that integrates v four times with the clamped boundary conditions
        B = PIE.T;
        % identity operator on L2[0,b]
        Id = get_op(y);
        % P v = u_ss -- since u_ss + u_ssss = A v and u_ssss = v, we have A = P + I, so P = PIE.A - I
        P = PIE.A - Id;
        % H1 v = u_s -- the derivative of B in s, (H1*u_ssss) = u_s
        H1 = diff(B,s1);

        % optimization problem min C s.t. C*P'*P - H1'*H1 >= 0, i.e. ||u_s||^2 <= C ||u_ss||^2
        prog = lpiprogram(B.vars,B.dom);
        [prog,C] = lpidecvar(prog,'C');
        if isempty(dg)
            prog = lpi_ineq(prog,C*(P'*P) - H1'*H1);
        else
            prog = lpi_ineq(prog,C*(P'*P) - H1'*H1,dg);
        end
        prog = lpisetobj(prog,C);
        prog = lpisolve(prog);
        Csol = double(lpigetsol(prog,C));

        % sharp value and the eigenvalue: lambda = -||u_s||^2/||u_ss||^2, so the largest one is lambda = -C
        sharp = (b/(2*pi))^2;
        fprintf('\n Largest eigenvalue lambda = -C = %.6f (exact -(b/(2*pi))^2 = %.6f).\n',-Csol,-sharp);
%__________________________________________________________________________
%==========================================================================
end

% Report the result
fprintf('\n %s\n   LPI bound C = %.6f',name,Csol);
if ~isnan(sharp)
    fprintf(';  sharp / reference value %.6f',sharp);
end
if ~isnan(paper)
    fprintf(';  paper [1] %.6f',paper);
end
fprintf('.\n');

% Specify the outputs of the function
info = struct('name',name,'constant',Csol,'sharp',sharp,'paper',paper,'extra',extra);
varargout{1} = Csol;
if nargout>=2
    varargout{2} = info;
end

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% References %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% % [1] -
% @article{barrazamartinez2021trace,
%   title={On trace theorems and Poincare inequality for one-dimensional Sobolev spaces},
%   author={Barraza Martinez, Bernardo and Gonzalez Ospino, Jorge and Hernandez Monzon, Jairo},
%   journal={arXiv preprint arXiv:2108.13926},
%   year={2021}
% }
%
% % [2] -
% @inproceedings{valmorbida2014sdp,
%   title={Semi-definite programming and functional inequalities for distributed parameter systems},
%   author={Valmorbida, Giorgio and Ahmadi, Mohamadreza and Papachristodoulou, Antonis},
%   booktitle={53rd IEEE Conference on Decision and Control},
%   address={Los Angeles, CA},
%   pages={4304--4309},
%   year={2014}
% }
%
% % [3] -
% @inproceedings{valmorbida2015intsostools,
%   title={Introducing {INTSOSTOOLS}: A {SOSTOOLS} plug-in for integral inequalities},
%   author={Valmorbida, Giorgio and Papachristodoulou, Antonis},
%   booktitle={European Control Conference},
%   address={Linz, Austria},
%   pages={1231--1236},
%   year={2015}
% }

function Csol = solve_min(dom,ineq,dg)
% CSOL = SOLVE_MIN(DOM,INEQ,DG) solves  min C  s.t.  ineq(C) >= 0  on the
% domain dom, with certificate degrees dg, and returns the optimal C.
% INPUT
% - dom:  1 x 2 'double', the interval of the variable s1;
% - ineq: function handle, C -> the operator ineq(C), affine in the
%         scalar decision variable C;
% - dg:   'struct' with fields dD, dw (the degrees passed to
%         lpi_ineq), or [] for the default degrees.
% OUTPUT
% - Csol: 'double', the solved value of C.
pvar s1 s1_dum
% Program in s1 (dummy s1_dum) with the decision variable C; minimize C
% subject to ineq(C) >= 0, solve and extract C
prog = lpiprogram([s1,s1_dum],dom);
[prog,C] = lpidecvar(prog,'C');
if isempty(dg)
    prog = lpi_ineq(prog,ineq(C));
else
    prog = lpi_ineq(prog,ineq(C),dg);
end
prog = lpisetobj(prog,C);
prog = lpisolve(prog);
Csol = double(lpigetsol(prog,C));
end

function dg = degrees(opts,dD,dw)
% DG = DEGREES(OPTS,DD,DW) returns the certificate degrees: the user's
% values if specified, else the defaults dD, dw of the example (empty: the
% default degrees of lpi_ineq).
% INPUT
% - opts: 'struct' of options, fields dD and dw ([] if not specified);
% - dD, dw: 'double' scalars or [], the defaults of the example.
% OUTPUT
% - dg:   'struct' with fields dD, dw, or [] if either degree is empty.
if ~isempty(opts.dD),   dD = opts.dD;   end
if ~isempty(opts.dw),   dw = opts.dw;   end
if isempty(dD) || isempty(dw)
    dg = [];
else
    dg = struct('dD',dD,'dw',dw);
end
end

function [index,opts] = process_inputs(varargin0,nargin1)
% [INDEX,OPTS] = PROCESS_INPUTS(VARARGIN0,NARGIN1) reads the example index
% and the 'name=value' options; with no input, it prompts for them.
% INPUT
% - varargin0: 'cell', the inputs of the library function;
% - nargin1:   'double', their number.
% OUTPUT
% - index: 'double', the example index, 1 to 15;
% - opts:  'struct' with fields b, dD, dw, eps, t, tmin: the
%          defaults, overridden by the recognized options; an option that
%          is not understood is ignored with a warning.

% Number of examples and default options
n_examples = 15;
opts = struct('b',1,'dD',[],'dw',[],'eps',1e-2,'t',[], ...
              'tmin',0.05);

% Collect the inputs
if nargin1==0 %<-- If there is no input, pause the script and let the user specify an example
    userinp = input(sprintf('\n Select an example (1 through %d) \n ---> ',n_examples),'s');
    varargin0 = split(userinp,[" ",","]).';
    varargin0{1} = str2double(varargin0{1});
end
if ~isnumeric(varargin0{1}) || ~isscalar(varargin0{1}) || ~any(varargin0{1}==1:n_examples)
    error('First argument must be an integer between 1 and %d, corresponding to an example listed in this file',n_examples)
end
index = varargin0{1};
% Parse each 'name=value' option with a numeric value and a known name
for j=2:numel(varargin0)
    tok = regexp(strtrim(varargin0{j}),'^(\w+)\s*=\s*(.+?);?$','tokens','once');
    if ischar(varargin0{j}) && ~isempty(tok) && isfield(opts,tok{1}) && ~isnan(str2double(tok{2}))
        opts.(tok{1}) = str2double(tok{2});
    else
        disp(['Warning: Argument ',num2str(j),' is not understood, and is therefore ignored']);
    end
end

end
