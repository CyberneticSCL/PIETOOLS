% This document illustrates how PIETOOLS can be used to declare PI
% operators in 1D and 2D, acting on coupled finite- and
% infinite-dimensional vectors.
% We refer to Chapter 3 of the manual for more context on the codes.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PIETOOLS - Code Illustrations
%
% Copyright (C)2024  PIETOOLS Team
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
% DJ, 12/30/2024: Initial coding;
clear

%% 3.1.1 Declaring 3-PI Operators

% Declare an opvar object representing the 3-PI operator
% A: L2^2[-1,1] --> L2^2[-1,1] defined by
%   (A*v)(s) = [int_{-1}^{s} v1(th) + 2*v2(th) dth]
%              [int_{-1}^{s} 3*v1(th) + 4*v2(th) dth]
% for v=(v1,v2) in L2^2[-1,1]
opvar A;
A.I = [-1,1];
A.R.R1 = [1,2; 3,4];

% Declare an opvar object representing the 3-PI operator
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
opvar B;
B.I = [0,1];
B.var1 = s;      B.var2 = s_dum;
B.R.R0 = R0;     B.R.R1 = R1;     B.R.R2 = R2;


%% 3.1.4 Declaring 4-PI Operators

% Declare an opvar object representing the 4-PI operator
% C: R^2 x L2[0,3] --> R x L2^2[0,3], defined by
%   (C*q)(s) = [-x1 +2*x2 +int_{0}^{3}(3-s^2)*v(s)ds                              ]
%              [-s*x2 +v(s) +int_{0}^{s}(s-th)*v(th)dth +int_{s}^{3}s*v(th)dth    ]
%              [s*x1 +s^3*v(s) +int_{0}^{s}th*v(th)dth +int_{s}^{3}(th-s)*v(th)dth]
% for q=(x,v1,v2) in R x L2^2[0,3].
pvar s tt
P = [-1,2];
Q1 = (3-s^2);
Q2 = [0,-s; s,0];
R0 = [1; s^3]; R1 = [s-tt; tt]; R2 = [s; tt-s];
opvar C;
C.I = [0,3];
C.var1 = s; C.var2 = tt;
C.P = P;
C.Q1 = Q1;
C.Q2 = Q2;
C.R.R0 = R0; C.R.R1 = R1; C.R.R2 = R2;


%% 3.2.1 Declaring 9-PI Operators

% Declare an opvar2d object representing the 9-PI operator
% D: L2^2[[0,1]x[1,2]] --> L2^2[[0,1]x[1,2]] defined by
%   D*v(x,y) = int_{0}^{x} int_{y}^{2} [x^2*v1(th,nu) +x*y*v2(th,nu)] dnu dth 
%                                      [x*y*v1(th,nu) +y^2*v2(th,nu)]
% for v = (v1,v2) in L2^2[[0,1]x[1,2]]
pvar s1 s2
N12 = [s1^2, s1*s2; s1*s2, s2^2];
opvar2d D;
D.var1 = [s1;s2];
D.I = [0,1; 1,2];
D.R22{2,3} = N12;

% Declare an opvar2d object representing the 9-PI operator
% E: L2[[0,1]x[-1,1]] --> L2[[0,1]x[-1,1]] defined by
%   E*v(x,y) = x^2*y^2*v(x,y) +int_{-1}^{y}x*(y-nu)*v(x,nu) dnu
%               +int_{x}^{1}(x-th)*y*v(th,y)dth
%                   +int_{x}^{1}int_{-1}^{y}(x-th)*(y-nu)*v(th,nu) dnu dth
% for v in L2[[0,1]x[-1,1]]
pvar s1 s2 th1 th2
N00 = s1^2 * s2^2;      N01 = s1*(s2-th2);
N20 = (s1-th1)*s2;      N21 = (s1-th1)*(s2-th2);
opvar2d E;
E.var1 = [s1;s2]; E.var2 = [th1; th2];
E.I = [0,1; -1,1];
E.R22{1,1} = N00; E.R22{1,2} = N01;
E.R22{3,1} = N20; E.R22{3,2} = N21;


%% 3.2.2 Declaring General 2D PI Operators
% Declare an opvar2d object representing the 2D PI operator
% F: R x L2[0,2] x L2[[0,2]x[2,3]] --> L2^2[0,2] x L2[[0,2]x[2,3]] defined by
%   F*v(x,y) 
%   = [v0 +x*v1(x) +int_{x}^{2}v1(th) dth +int_{0}^{x}int_{2}^{3}y*v(th,y) dy dth] 
%     [x*v0 +x^2*v1(x) +int_{x}^{2}(th-x)*v1(th)dth +int_{0}^{x}int_{2}^{3}y^2(x-th)*v(th,y) dy dth] 
%     [y^2*v1(x) +int_{0}^{x}y*v1(th)dth +int_{0}^{x}int_{2}^{y}th*nu*v2(th,nu) dth dnu]
% for v=(v0,v1,v2) in R x L2[0,2] x L2[[0,2]x[2,3]]

pvar x y theta nu
Rx0 = [1; x];
Rxx_0 = [x; x^2]; Rxx_2 = [1; theta-x];
Rx2_1 = [y; y^2 * (x-theta)];
R2x_0 = y^2; R2x_1 = y;
R22_11 = theta*nu;
opvar2d F;
F.var1 = [x; y]; F.var2 = [theta; nu];
F.I = [0,2; 2,3];
F.Rx0 = Rx0;
F.Rxx{1} = Rxx_0; F.Rxx{3} = Rxx_2;
F.Rx2{2} = Rx2_1;
F.R2x{1} = R2x_0; F.R2x{2} = R2x_1;
F.R22{2,2} = R22_11;