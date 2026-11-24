% This document illustrates how PIETOOLS can be used to declare PI
% operators in 1D and 2D, acting on coupled finite- and
% infinite-dimensional vectors.
% The script declares, as 'copvar' objects, the 3-PI operators A and B and
% the 4-PI operator C on intervals, the 9-PI operators D and E on
% rectangles, and the 2D PI operator F, whose blocks map between R, L2[x]
% and L2[x,y]. Each operator is extracted (get_op) from the PDE terms that
% it produces when applied to PDE variables. It prints nothing.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PIETOOLS - Code Illustrations
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
% If you modify this code, make sure to change the code in the manual as
% well, and vice versa. Document all changes carefully and include date
% authorship, and a brief description of modifications
%
% For support, contact M. Peet, Arizona State University at mpeet@asu.edu
% DJ, 12/30/2024: Initial coding;
% DJ, 11/23/2026: Declare the operators from PDE terms with get_op;

% Clear the workspace
clear

%% 3.1.1 Declaring 3-PI Operators

% Declare a copvar object representing the 3-PI operator
% A: L2^2[-1,1] --> L2^2[-1,1] defined by
%   (A*v)(s) = [int_{-1}^{s} v1(th) + 2*v2(th) dth]
%              [int_{-1}^{s} 3*v1(th) + 4*v2(th) dth]
% for v=(v1,v2) in L2^2[-1,1]
% The operator is applied to a PDE variable v in L2^2[-1,1], and extracted
% from the resulting terms with "get_op"; v(s_dum) is v evaluated at the
% dummy variable of integration s_dum
pvar s s_dum
v = pde_var(2,s,[-1,1]);                                                    % DJ, 11/23/2026
A = get_op(int([1,2; 3,4]*v(s_dum),s_dum,[-1,s]));

% Declare a copvar object representing the 3-PI operator
% B: L2[0,1] --> L2^2[0,1] defined by
%   (B*v)(s) = [v(s)    ] + int_{0}^{s} [2s v(th)    ] dth
%              [s^2 v(s)]               [s(s-th)v(th)]
%                               + int_{s}^{1} [3th v(th)       ] dth
%                                             [3/4 (s^2-s)v(th)]
% for v in L2[0,1].
pvar s s_dum
R0 = [1; s^2];
R1 = [2*s; s*(s-s_dum)];
R2 = [3*s_dum; (3/4)*(s^2-s)];
v = pde_var(s,[0,1]);                                                       % DJ, 11/23/2026
B = get_op(R0*v + int(R1*v(s_dum),s_dum,[0,s]) ...
                + int(R2*v(s_dum),s_dum,[s,1]));


%% 3.1.4 Declaring 4-PI Operators

% Declare a copvar object representing the 4-PI operator
% C: R^2 x L2[0,3] --> R x L2^2[0,3], defined by
%   (C*q)(s) = [-x1 +2*x2 +int_{0}^{3}(3-s^2)*v(s)ds                              ]
%              [-s*x2 +v(s) +int_{0}^{s}(s-th)*v(th)dth +int_{s}^{3}s*v(th)dth    ]
%              [s*x1 +s^3*v(s) +int_{0}^{s}th*v(th)dth +int_{s}^{3}(th-s)*v(th)dth]
% for q=(x1,x2,v) in R^2 x L2[0,3].
P = [-1,2];
Q1 = (3-s^2);
Q2 = [0,-s; s,0];
R0 = [1; s^3]; R1 = [s-s_dum; s_dum]; R2 = [s; s_dum-s];   % dummy variable s_dum
% q = (x,v), with x in R^2 and v in L2[0,3]; the first row of terms is
% the output in R, the second that in L2^2[0,3]
x = pde_var(2);     v = pde_var(s,[0,3]);                                   % DJ, 11/23/2026
C = get_op([P*x + int(Q1*v,s,[0,3]);
            Q2*x + R0*v + int(R1*v(s_dum),s_dum,[0,s]) + int(R2*v(s_dum),s_dum,[s,3])]);


%% 3.2.1 Declaring 9-PI Operators

% Declare a copvar object representing the 9-PI operator
% D: L2^2[[0,1]x[1,2]] --> L2^2[[0,1]x[1,2]] defined by
%   D*v(x,y) = int_{0}^{x} int_{y}^{2} [x^2*v1(th,nu) +x*y*v2(th,nu)] dnu dth 
%                                      [x*y*v1(th,nu) +y^2*v2(th,nu)]
% for v = (v1,v2) in L2^2[[0,1]x[1,2]]
pvar s1 s2 s1_dum s2_dum
N12 = [s1^2, s1*s2; s1*s2, s2^2];
v = pde_var(2,[s1;s2],[0,1; 1,2]);
D = get_op(int(N12*v(s1_dum,s2_dum),[s1_dum;s2_dum],[0,s1; s2,2]));         % DJ, 11/23/2026

% Declare a copvar object representing the 9-PI operator
% E: L2[[0,1]x[-1,1]] --> L2[[0,1]x[-1,1]] defined by
%   E*v(x,y) = x^2*y^2*v(x,y) +int_{-1}^{y}x*(y-nu)*v(x,nu) dnu
%               +int_{x}^{1}(x-th)*y*v(th,y)dth
%                   +int_{x}^{1}int_{-1}^{y}(x-th)*(y-nu)*v(th,nu) dnu dth
% for v in L2[[0,1]x[-1,1]]
pvar s1 s2 s1_dum s2_dum
N00 = s1^2 * s2^2;      N01 = s1*(s2-s2_dum);
N20 = (s1-s1_dum)*s2;   N21 = (s1-s1_dum)*(s2-s2_dum);
v = pde_var([s1;s2],[0,1; -1,1]);                                           % DJ, 11/23/2026
E = get_op(N00*v + int(N01*v(s1,s2_dum),s2_dum,[-1,s2]) ...
         + int(N20*v(s1_dum,s2),s1_dum,[s1,1]) ...
         + int(N21*v(s1_dum,s2_dum),[s1_dum;s2_dum],[s1,1; -1,s2]));


%% 3.2.2 Declaring General 2D PI Operators
% Declare a copvar object representing the 2D PI operator
% F: R x L2[0,2] x L2[[0,2]x[2,3]] --> L2^2[0,2] x L2[[0,2]x[2,3]] defined by
%   F*v(x,y) 
%   = [v0 +x*v1(x) +int_{x}^{2}v1(th) dth +int_{0}^{x}int_{2}^{3}y*v(th,y) dy dth] 
%     [x*v0 +x^2*v1(x) +int_{x}^{2}(th-x)*v1(th)dth +int_{0}^{x}int_{2}^{3}y^2(x-th)*v(th,y) dy dth] 
%     [y^2*v1(x) +int_{0}^{x}y*v1(th)dth +int_{0}^{x}int_{2}^{y}th*nu*v2(th,nu) dth dnu]
% for v=(v0,v1,v2) in R x L2[0,2] x L2[[0,2]x[2,3]]
pvar x y x_dum y_dum
Rx0 = [1; x];
Rxx_0 = [x; x^2]; Rxx_2 = [1; x_dum-x];
Rx2_1 = [y; y^2 * (x-x_dum)];
R2x_0 = y^2; R2x_1 = y;
R22_11 = x_dum*y_dum;
% v = (v0,v1,v2); the first row of terms is the output in L2^2[0,2], the
% second that in L2[[0,2]x[2,3]]; the kernel Rx2_1 is integrated over y
v0 = pde_var(1);  v1 = pde_var(x,[0,2]);  v2 = pde_var([x;y],[0,2; 2,3]);  % DJ, 11/23/2026
F = get_op([Rx0*v0 + Rxx_0*v1 + int(Rxx_2*v1(x_dum),x_dum,[x,2]) ...
                + int(Rx2_1*v2(x_dum,y),[x_dum;y],[0,x; 2,3]);
            R2x_0*v1 + int(R2x_1*v1(x_dum),x_dum,[0,x]) ...
                + int(R22_11*v2(x_dum,y_dum),[x_dum;y_dum],[0,x; 2,y])]);