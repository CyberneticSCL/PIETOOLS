% This document illustrates how PIETOOLS can be used to declare PI
% operators in 1D and 2D, acting on coupled finite- and
% infinite-dimensional vectors.
% We refer to Chapter 3 of the manual for more context on the codes.

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
% DJ, 12/30/2024: Initial coding;
% DJ, 11/08/2026: Declare the 2D operators D, E and F (3.2.1, 3.2.2) as
%                   copvar objects, with int2copvar and int2sopvar;
clear

%% 3.1.1 Declaring 3-PI Operators

% Declare an opvar object representing the 3-PI operator
% A: L2^2[-1,1] --> L2^2[-1,1] defined by
%   (A*v)(s) = [int_{-1}^{s} v1(th) + 2*v2(th) dth]
%              [int_{-1}^{s} 3*v1(th) + 4*v2(th) dth]
% for v=(v1,v2) in L2^2[-1,1]
% opvar A;                                                                  % DB, 11/06/2026 (was)
% A.I = [-1,1];                                                             % DB, 11/06/2026 (was)
% A.R.R1 = [1,2; 3,4];                                                      % DB, 11/06/2026 (was)
% copvar objects are declared with "pi4_copvar_sop(P,Q1,Q2,R,dims,s,s_dum,dom)",% DB, 11/06/2026
% where dims = [q0,p0; q1,p1] are the dimensions of the R and L2 parts of the% DB, 11/06/2026
% output (q) and input (p), and R = {R0,R1,R2} (with [] for a zero parameter)% DB, 11/06/2026
pvar s s_dum                                                                % DB, 11/06/2026
A = pi4_copvar_sop([],[],[],{[],[1,2; 3,4],[]},[0 0; 2 2],s,s_dum,[-1,1]);  % DB, 11/06/2026

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
% opvar B;                                                                  % DB, 11/06/2026 (was)
% B.I = [0,1];                                                              % DB, 11/06/2026 (was)
% B.var1 = s;      B.var2 = s_dum;                                          % DB, 11/06/2026 (was)
% B.R.R0 = R0;     B.R.R1 = R1;     B.R.R2 = R2;                            % DB, 11/06/2026 (was)
B = pi4_copvar_sop([],[],[],{R0,R1,R2},[0 0; 2 1],s,s_dum,[0,1]);           % DB, 11/06/2026


%% 3.1.4 Declaring 4-PI Operators

% Declare an opvar object representing the 4-PI operator
% C: R^2 x L2[0,3] --> R x L2^2[0,3], defined by
%   (C*q)(s) = [-x1 +2*x2 +int_{0}^{3}(3-s^2)*v(s)ds                              ]
%              [-s*x2 +v(s) +int_{0}^{s}(s-th)*v(th)dth +int_{s}^{3}s*v(th)dth    ]
%              [s*x1 +s^3*v(s) +int_{0}^{s}th*v(th)dth +int_{s}^{3}(th-s)*v(th)dth]
% for q=(x,v1,v2) in R x L2^2[0,3].
% pvar s tt                                                                 % DB, 11/06/2026 (was)
% P = [-1,2];                                                               % DB, 11/06/2026 (was)
% Q1 = (3-s^2);                                                             % DB, 11/06/2026 (was)
% Q2 = [0,-s; s,0];                                                         % DB, 11/06/2026 (was)
% R0 = [1; s^3]; R1 = [s-tt; tt]; R2 = [s; tt-s];                           % DB, 11/06/2026 (was)
% opvar C;                                                                  % DB, 11/06/2026 (was)
% C.I = [0,3];                                                              % DB, 11/06/2026 (was)
% C.var1 = s; C.var2 = tt;                                                  % DB, 11/06/2026 (was)
% C.P = P;                                                                  % DB, 11/06/2026 (was)
% C.Q1 = Q1;                                                                % DB, 11/06/2026 (was)
% C.Q2 = Q2;                                                                % DB, 11/06/2026 (was)
% C.R.R0 = R0; C.R.R1 = R1; C.R.R2 = R2;                                    % DB, 11/06/2026 (was)
P = [-1,2];                                                                 % DB, 11/06/2026
Q1 = (3-s^2);                                                               % DB, 11/06/2026
Q2 = [0,-s; s,0];                                                           % DB, 11/06/2026
R0 = [1; s^3]; R1 = [s-s_dum; s_dum]; R2 = [s; s_dum-s];   % dummy variable s_dum% DB, 11/06/2026
C = pi4_copvar_sop(P,Q1,Q2,{R0,R1,R2},[1 2; 2 1],s,s_dum,[0,3]);            % DB, 11/06/2026


%% 3.2.1 Declaring 9-PI Operators

% Declare a copvar object representing the 9-PI operator
% D: L2^2[[0,1]x[1,2]] --> L2^2[[0,1]x[1,2]] defined by
%   D*v(x,y) = int_{0}^{x} int_{y}^{2} [x^2*v1(th,nu) +x*y*v2(th,nu)] dnu dth 
%                                      [x*y*v1(th,nu) +y^2*v2(th,nu)]
% for v = (v1,v2) in L2^2[[0,1]x[1,2]]
pvar s1 s2
N12 = [s1^2, s1*s2; s1*s2, s2^2];
D = int2copvar(N12,{'s1','s2'},[0,1; 1,2],[0,s1; s2,2]);

% Declare a copvar object representing the 9-PI operator
% E: L2[[0,1]x[-1,1]] --> L2[[0,1]x[-1,1]] defined by
%   E*v(x,y) = x^2*y^2*v(x,y) +int_{-1}^{y}x*(y-nu)*v(x,nu) dnu
%               +int_{x}^{1}(x-th)*y*v(th,y)dth
%                   +int_{x}^{1}int_{-1}^{y}(x-th)*(y-nu)*v(th,nu) dnu dth
% for v in L2[[0,1]x[-1,1]]
% pvar s1 s2 th1 th2                                                        % DB, 11/06/2026 (was)
% N00 = s1^2 * s2^2;      N01 = s1*(s2-th2);                                % DB, 11/06/2026 (was)
% N20 = (s1-th1)*s2;      N21 = (s1-th1)*(s2-th2);                          % DB, 11/06/2026 (was)
% (for a copvar, the dummy variables must be named s1_dum, s2_dum)          % DB, 11/06/2026
pvar s1 s2 s1_dum s2_dum                                                    % DB, 11/06/2026
N00 = s1^2 * s2^2;      N01 = s1*(s2-s2_dum);                               % DB, 11/06/2026
N20 = (s1-s1_dum)*s2;   N21 = (s1-s1_dum)*(s2-s2_dum);                      % DB, 11/06/2026
vE = {'s1','s2'};       domE = [0,1; -1,1];
E = int2copvar(N00,vE,domE,[s1,s1; s2,s2]) + int2copvar(N01,vE,domE,[s1,s1; -1,s2]) ...
    + int2copvar(N20,vE,domE,[s1,1; s2,s2]) + int2copvar(N21,vE,domE,[s1,1; -1,s2]);


%% 3.2.2 Declaring General 2D PI Operators
% Declare a copvar object representing the 2D PI operator
% F: R x L2[0,2] x L2[[0,2]x[2,3]] --> L2^2[0,2] x L2[[0,2]x[2,3]] defined by
%   F*v(x,y) 
%   = [v0 +x*v1(x) +int_{x}^{2}v1(th) dth +int_{0}^{x}int_{2}^{3}y*v(th,y) dy dth] 
%     [x*v0 +x^2*v1(x) +int_{x}^{2}(th-x)*v1(th)dth +int_{0}^{x}int_{2}^{3}y^2(x-th)*v(th,y) dy dth] 
%     [y^2*v1(x) +int_{0}^{x}y*v1(th)dth +int_{0}^{x}int_{2}^{y}th*nu*v2(th,nu) dth dnu]
% for v=(v0,v1,v2) in R x L2[0,2] x L2[[0,2]x[2,3]]

% pvar x y theta nu                                                         % DB, 11/06/2026 (was)
% Rx0 = [1; x];                                                             % DB, 11/06/2026 (was)
% Rxx_0 = [x; x^2]; Rxx_2 = [1; theta-x];                                   % DB, 11/06/2026 (was)
% Rx2_1 = [y; y^2 * (x-theta)];                                             % DB, 11/06/2026 (was)
% R2x_0 = y^2; R2x_1 = y;                                                   % DB, 11/06/2026 (was)
% R22_11 = theta*nu;                                                        % DB, 11/06/2026 (was)
% (for a copvar, the dummy variables must be named x_dum, y_dum)            % DB, 11/06/2026
pvar x y x_dum y_dum                                                        % DB, 11/06/2026
Rx0 = [1; x];                                                               % DB, 11/06/2026
Rxx_0 = [x; x^2]; Rxx_2 = [1; x_dum-x];                                     % DB, 11/06/2026
Rx2_1 = [y; y^2 * (x-x_dum)];                                               % DB, 11/06/2026
R2x_0 = y^2; R2x_1 = y;                                                     % DB, 11/06/2026
R22_11 = x_dum*y_dum;                                                       % DB, 11/06/2026
dx = [0,2];     dy = [2,3];
sp = @(vin,vout) struct('in',{vin},'out',{vout});       % variables of the input and output space
dm = @(din,dout) struct('in',din,'out',dout);           % their domains
F11 = int2sopvar(Rx0,sp({},{'x'}),dm(zeros(0,2),dx));                                   % R --> L2^2[x]
F12 = int2sopvar(Rxx_0,sp({'x'},{'x'}),dm(dx,dx),[x,x]) ...
    + int2sopvar(Rxx_2,sp({'x'},{'x'}),dm(dx,dx),[x,2]);                                % L2[x] --> L2^2[x]
F13 = int2sopvar(Rx2_1,sp({'x','y'},{'x'}),dm([dx;dy],dx),[0,x; 2,3]);                  % L2[x,y] --> L2^2[x]
F22 = int2sopvar(R2x_0,sp({'x'},{'x','y'}),dm(dx,[dx;dy]),[x,x]) ...
    + int2sopvar(R2x_1,sp({'x'},{'x','y'}),dm(dx,[dx;dy]),[0,x]);                       % L2[x] --> L2[x,y]
F23 = int2sopvar(R22_11,sp({'x','y'},{'x','y'}),dm([dx;dy],[dx;dy]),[0,x; 2,y]);        % L2[x,y] --> L2[x,y]
F = copvar({F11, F12, F13; [], F22, F23});