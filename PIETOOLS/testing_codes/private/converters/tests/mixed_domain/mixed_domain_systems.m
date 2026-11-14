function S = mixed_domain_systems(name,p)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% S = MIXED_DOMAIN_SYSTEMS(NAME,P) returns one test system that couples
% states of different spatial dimension (R, L2[s1], L2[s2], L2[s1,s2],
% L2[s1,s2,s3]): its PDE, and its semantics written by hand, independently
% of the parser. S = MIXED_DOMAIN_SYSTEMS('list') returns the catalogue.
%
% INPUT
% - name:   char, the name of a system (see the 'list' output), or 'list';
% - p:      (optional) 'struct' of parameter values overriding the defaults
%           of the system (fields as named in each case below).
%
% OUTPUT
% - S:      for a system, 'struct' with fields
%             name     char, NAME;
%             PDE      'pde_struct', the PDE (not initialized);
%             svars    1 x nx cell, svars{k} the cellstr of the variables of
%                      state k ('s1','s2','s3');
%             wvars    1 x nw cell, the same for the inputs;
%             xs       function handle @(c,W): 1 x nx cell of 'polynomial'
%                      states satisfying the BCs for the inputs W (1 x nw
%                      cell of 'polynomial'), c a vector of 40 coefficients;
%             rhs      @(X,W): 1 x nx cell, the right-hand side of each
%                      state equation at the states X and inputs W;
%             bc       @(X,W): cell of 'polynomial' residuals that vanish
%                      iff X satisfies the BCs;
%             z        @(X,W): cell, the regulated outputs ({} if none);
%             expect   char, '' for a system the converter must accept,
%                      else a keyword of the error it must raise;
%             par      'struct', P;
%           for 'list', a struct array with fields name, spaces, what and
%           baseline, one element per system.
%
% For support, contact M. Peet, Arizona State University at mpeet@asu.edu

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PIETOOLS - mixed_domain_systems
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
% Initial coding MMP, 11/12/2026
% MMP, 11/14/2026: S10, S11, S13 expect 'constrain a state of lower
%   dimension', the refusal of the converter since dbf458a2 (Declan), which
%   no longer names the compatibility conditions; S19 converts since
%   dbf458a2 (D2 fixed), its row updated.

if nargin<2, p = struct(); end
if strcmp(name,'list')
    S = catalogue();
    return
end
pvar s1 s2 s3 t
% Derivative of order k, and polynomial families with prescribed BCs on [0,1]
D = @(q,v,k) dk(q,v,k);
nd  = @(s,c) c(1)+c(2)*(2*s^3-3*s^2)+s^2*(1-s)^2*(c(3)+c(4)*s);        % q'(0)=q'(1)=0
nD  = @(s,c) c(1)*(1-s)^2*(1+2*s)+s^2*(1-s)*(c(2)+c(3)*s);            % q'(0)=0, q(1)=0
dd  = @(s,c) s*(1-s)*(c(1)+c(2)*s+c(3)*s^2);                          % q(0)=q(1)=0
S = struct('name',name,'PDE',[],'svars',{{}},'wvars',{{}},'xs',[],'rhs',[],'bc',[], ...
           'z',[],'expect','','par',p);
S.z = @(X,W) {};
clear stateNameGenerator
switch name
    case 'S1_ode_heat1d_bc'
        % ODE -> 1-D heat through a Dirichlet BC; X' = -a X, reaction r
        a = getp(p,'a',1);  r = getp(p,'r',0);
        X = pde_var(1);  x1 = pde_var(1,s1,[0,1]);
        S.PDE = [diff(X,t)==-a*X; diff(x1,t)==diff(x1,s1,2)+r*x1;
                 subs(x1,s1,0)==X; subs(x1,s1,1)==0];
        S.svars = {{},{'s1'}};
        S.xs  = @(c,W) {c(1), c(1)*(1-s1)+dd(s1,c(2:4))};
        S.rhs = @(Xp,W) {-a*Xp{1}, D(Xp{2},s1,2)+r*Xp{2}};
        S.bc  = @(Xp,W) {subs(Xp{2},s1,0)-Xp{1}, subs(Xp{2},s1,1)};
    case 'S2_ode_heat2d_edge'
        % ODE -> 2-D heat through the edge s2=0, Neumann in s1
        a = getp(p,'a',1);
        X = pde_var(1);  x2 = pde_var(1,[s1;s2],[0,1;0,1]);
        S.PDE = [diff(X,t)==-a*X; diff(x2,t)==diff(x2,s1,2)+diff(x2,s2,2);
                 subs(diff(x2,s1),s1,0)==0; subs(diff(x2,s1),s1,1)==0;
                 subs(x2,s2,0)==X; subs(x2,s2,1)==0];
        S.svars = {{},{'s1','s2'}};
        S.xs  = @(c,W) {c(1), c(1)*(1-s2)+s2*(1-s2)*(nd(s1,c(2:5))*(c(6)+c(7)*s2)+nd(s1,c(8:11))*s2^2)};
        S.rhs = @(Xp,W) {-a*Xp{1}, D(Xp{2},s1,2)+D(Xp{2},s2,2)};
        S.bc  = @(Xp,W) {subs(D(Xp{2},s1,1),s1,0), subs(D(Xp{2},s1,1),s1,1), subs(Xp{2},s2,0)-Xp{1}, subs(Xp{2},s2,1)};
    case 'S3_edge1d_heat2d'
        % 1-D heat on the edge s2=0 of a 2-D heat; reaction r in x1
        r = getp(p,'r',0);
        x1 = pde_var(1,s1,[0,1]);  x2 = pde_var(1,[s1;s2],[0,1;0,1]);
        S.PDE = [diff(x1,t)==diff(x1,s1,2)+r*x1; diff(x2,t)==diff(x2,s1,2)+diff(x2,s2,2);
                 subs(x1,s1,0)==0; subs(x1,s1,1)==0;
                 subs(x2,s1,0)==0; subs(x2,s1,1)==0;
                 subs(x2,s2,0)==x1; subs(x2,s2,1)==0];
        S.svars = {{'s1'},{'s1','s2'}};
        S.xs  = @(c,W) {dd(s1,c(1:3)), (1-s2)*dd(s1,c(1:3))+s1*(1-s1)*s2*(1-s2)*(c(4)+c(5)*s1+c(6)*s2+c(7)*s1*s2^2)};
        S.rhs = @(Xp,W) {D(Xp{1},s1,2)+r*Xp{1}, D(Xp{2},s1,2)+D(Xp{2},s2,2)};
        S.bc  = @(Xp,W) {subs(Xp{1},s1,0), subs(Xp{1},s1,1), subs(Xp{2},s1,0), subs(Xp{2},s1,1), subs(Xp{2},s2,0)-Xp{1}, subs(Xp{2},s2,1)};
    case 'S4_wentzell'
        % 2-D heat with a dynamic (Wentzell) condition on s2=0: x1 = x2(s1,0),
        % x1_t = x1_s1s1 + c x2_s2(s1,0) + r x1,  x2_t = lap x2 + r x2
        c0 = getp(p,'c',1);  r = getp(p,'r',0);
        x1 = pde_var(1,s1,[0,1]);  x2 = pde_var(1,[s1;s2],[0,1;0,1]);
        S.PDE = [diff(x1,t)==diff(x1,s1,2)+c0*subs(diff(x2,s2),s2,0)+r*x1;
                 diff(x2,t)==diff(x2,s1,2)+diff(x2,s2,2)+r*x2;
                 subs(x1,s1,0)==0; subs(x1,s1,1)==0;
                 subs(x2,s1,0)==0; subs(x2,s1,1)==0;
                 subs(x2,s2,0)==x1; subs(x2,s2,1)==0];
        S.svars = {{'s1'},{'s1','s2'}};
        S.xs  = @(c,W) {dd(s1,c(1:3)), (1-s2)*dd(s1,c(1:3))+s1*(1-s1)*s2*(1-s2)*(c(4)+c(5)*s1+c(6)*s2+c(7)*s1*s2^2)};
        S.rhs = @(Xp,W) {D(Xp{1},s1,2)+c0*subs(D(Xp{2},s2,1),s2,0)+r*Xp{1}, D(Xp{2},s1,2)+D(Xp{2},s2,2)+r*Xp{2}};
        S.bc  = @(Xp,W) {subs(Xp{1},s1,0), subs(Xp{1},s1,1), subs(Xp{2},s1,0), subs(Xp{2},s1,1), subs(Xp{2},s2,0)-Xp{1}, subs(Xp{2},s2,1)};
    case 'S5_chain_ode_1d_2d'
        % ODE -> 1-D (Neumann data) -> 2-D (edge s2=0), and the ODE in the
        % s1=0 condition of the 2-D state
        a = getp(p,'a',1);
        X = pde_var(1);  x1 = pde_var(1,s1,[0,1]);  x2 = pde_var(1,[s1;s2],[0,1;0,1]);
        S.PDE = [diff(X,t)==-a*X; diff(x1,t)==diff(x1,s1,2);
                 diff(x2,t)==diff(x2,s1,2)+diff(x2,s2,2);
                 subs(diff(x1,s1),s1,0)==X; subs(x1,s1,1)==0;
                 subs(diff(x2,s1),s1,0)==(1-s2)*X; subs(x2,s1,1)==0;
                 subs(x2,s2,0)==x1; subs(x2,s2,1)==0];
        S.svars = {{},{'s1'},{'s1','s2'}};
        S.xs  = @(c,W) states_S5(c,s1,s2,nD);
        S.rhs = @(Xp,W) {-a*Xp{1}, D(Xp{2},s1,2), D(Xp{3},s1,2)+D(Xp{3},s2,2)};
        S.bc  = @(Xp,W) {subs(D(Xp{2},s1,1),s1,0)-Xp{1}, subs(Xp{2},s1,1), ...
                       subs(D(Xp{3},s1,1),s1,0)-(1-s2)*Xp{1}, subs(Xp{3},s1,1), ...
                       subs(Xp{3},s2,0)-Xp{2}, subs(Xp{3},s2,1)};
    case 'S6_indomain_1d_2d'
        % 1-D <-> 2-D in-domain coupling (integral over s2, multiplier),
        % Dirichlet in s1, Neumann in s2
        g = getp(p,'g',1);
        x1 = pde_var(1,s1,[0,1]);  x2 = pde_var(1,[s1;s2],[0,1;0,1]);
        S.PDE = [diff(x1,t)==diff(x1,s1,2)+g*int(x2,s2,[0,1]);
                 diff(x2,t)==diff(x2,s1,2)+diff(x2,s2,2)+g*x1;
                 subs(x1,s1,0)==0; subs(x1,s1,1)==0;
                 subs(x2,s1,0)==0; subs(x2,s1,1)==0;
                 subs(diff(x2,s2),s2,0)==0; subs(diff(x2,s2),s2,1)==0];
        S.svars = {{'s1'},{'s1','s2'}};
        S.xs  = @(c,W) {dd(s1,c(1:3)), dd(s1,c(4:6))*nd(s2,c(7:10))+dd(s1,c(11:13))*nd(s2,c(14:17))};
        S.rhs = @(Xp,W) {D(Xp{1},s1,2)+g*int(Xp{2},s2,0,1), D(Xp{2},s1,2)+D(Xp{2},s2,2)+g*Xp{1}};
        S.bc  = @(Xp,W) {subs(Xp{1},s1,0), subs(Xp{1},s1,1), subs(Xp{2},s1,0), subs(Xp{2},s1,1), ...
                       subs(D(Xp{2},s2,1),s2,0), subs(D(Xp{2},s2,1),s2,1)};
    case 'S7_ode_2d_indomain'
        % ODE <-> 2-D in-domain coupling, Neumann on every edge
        a = getp(p,'a',1);  b = getp(p,'b',1);  g = getp(p,'g',0.5);
        X = pde_var(1);  x2 = pde_var(1,[s1;s2],[0,1;0,1]);
        S.PDE = [diff(X,t)==-a*X+g*int(x2,[s1;s2],[0,1;0,1]);
                 diff(x2,t)==diff(x2,s1,2)+diff(x2,s2,2)-b*x2+g*X;
                 subs(diff(x2,s1),s1,0)==0; subs(diff(x2,s1),s1,1)==0;
                 subs(diff(x2,s2),s2,0)==0; subs(diff(x2,s2),s2,1)==0];
        S.svars = {{},{'s1','s2'}};
        S.xs  = @(c,W) {c(1), nd(s1,c(2:5))*nd(s2,c(6:9))+nd(s1,c(10:13))*nd(s2,c(14:17))};
        S.rhs = @(Xp,W) {-a*Xp{1}+g*int(int(Xp{2},s1,0,1),s2,0,1), D(Xp{2},s1,2)+D(Xp{2},s2,2)-b*Xp{2}+g*Xp{1}};
        S.bc  = @(Xp,W) {subs(D(Xp{2},s1,1),s1,0), subs(D(Xp{2},s1,1),s1,1), subs(D(Xp{2},s2,1),s2,0), subs(D(Xp{2},s2,1),s2,1)};
    case 'S8_two_edges'
        % 1-D states on the edges s2=0 (in s1) and s1=1 (in s2) of a 2-D
        % state; reaction r3 in the s2 edge state
        r3 = getp(p,'r3',0);
        x1 = pde_var(1,s1,[0,1]);  x3 = pde_var(1,s2,[0,1]);  x2 = pde_var(1,[s1;s2],[0,1;0,1]);
        S.PDE = [diff(x1,t)==diff(x1,s1,2); diff(x3,t)==diff(x3,s2,2)+r3*x3;
                 diff(x2,t)==diff(x2,s1,2)+diff(x2,s2,2);
                 subs(x1,s1,0)==0; subs(x1,s1,1)==0;
                 subs(x3,s2,0)==0; subs(x3,s2,1)==0;
                 subs(x2,s1,0)==0; subs(x2,s1,1)==x3;
                 subs(x2,s2,0)==x1; subs(x2,s2,1)==0];
        S.svars = {{'s1'},{'s2'},{'s1','s2'}};
        S.xs  = @(c,W) {dd(s1,c(1:3)), dd(s2,c(4:6)), (1-s2)*dd(s1,c(1:3))+s1*dd(s2,c(4:6))+s1*(1-s1)*s2*(1-s2)*(c(7)+c(8)*s1*s2)};
        S.rhs = @(Xp,W) {D(Xp{1},s1,2), D(Xp{2},s2,2)+r3*Xp{2}, D(Xp{3},s1,2)+D(Xp{3},s2,2)};
        S.bc  = @(Xp,W) {subs(Xp{1},s1,0), subs(Xp{1},s1,1), subs(Xp{2},s2,0), subs(Xp{2},s2,1), ...
                       subs(Xp{3},s1,0), subs(Xp{3},s1,1)-Xp{2}, subs(Xp{3},s2,0)-Xp{1}, subs(Xp{3},s2,1)};
    case 'S9_ode_2d_flux'
        % ODE <-> 2-D: the ODE fixes the edge value x2(s1,0) = X, the edge
        % flux drives the ODE; w into the ODE, z = X
        a = getp(p,'a',0.5);  c0 = getp(p,'c',1);
        X = pde_var(1);  x2 = pde_var(1,[s1;s2],[0,1;0,1]);
        w = pde_var('in',1);  z = pde_var('out',1);
        S.PDE = [diff(X,t)==a*X+c0*int(subs(diff(x2,s2),s2,0),s1,[0,1])+w;
                 diff(x2,t)==diff(x2,s1,2)+diff(x2,s2,2);
                 z==X;
                 subs(diff(x2,s1),s1,0)==0; subs(diff(x2,s1),s1,1)==0;
                 subs(x2,s2,0)==X; subs(x2,s2,1)==0];
        S.svars = {{},{'s1','s2'}};  S.wvars = {{}};
        S.xs  = @(c,W) {c(1), c(1)*(1-s2)+s2*(1-s2)*(nd(s1,c(2:5))*(c(6)+c(7)*s2)+nd(s1,c(8:11))*s2^2)};
        S.rhs = @(Xp,W) {a*Xp{1}+c0*int(subs(D(Xp{2},s2,1),s2,0),s1,0,1)+W{1}, D(Xp{2},s1,2)+D(Xp{2},s2,2)};
        S.bc  = @(Xp,W) {subs(D(Xp{2},s1,1),s1,0), subs(D(Xp{2},s1,1),s1,1), subs(Xp{2},s2,0)-Xp{1}, subs(Xp{2},s2,1)};
        S.z   = @(Xp,W) {Xp{1}};
    case 'S10_incompatible_corner'
        % the ODE as Dirichlet data on the edge s2=0 of a Dirichlet square:
        % the corner values clash unless X = 0
        X = pde_var(1);  x2 = pde_var(1,[s1;s2],[0,1;0,1]);
        S.PDE = [diff(X,t)==-X; diff(x2,t)==diff(x2,s1,2)+diff(x2,s2,2);
                 subs(x2,s1,0)==0; subs(x2,s1,1)==0;
                 subs(x2,s2,0)==X; subs(x2,s2,1)==0];
        S.svars = {{},{'s1','s2'}};
%       S.expect = 'compatibility';                                         % MMP, 11/14/2026 (was)
        S.expect = 'constrain a state of lower dimension';                  % MMP, 11/14/2026
    case 'S11_incompatible_neumann_corner'
        % x2_s1(0,s2)=0 and x2(s1,0)=x1(s1) require x1'(0)=0, which the
        % Dirichlet conditions of x1 leave free
        x1 = pde_var(1,s1,[0,1]);  x2 = pde_var(1,[s1;s2],[0,1;0,1]);
        S.PDE = [diff(x1,t)==diff(x1,s1,2); diff(x2,t)==diff(x2,s1,2)+diff(x2,s2,2);
                 subs(x1,s1,0)==0; subs(x1,s1,1)==0;
                 subs(diff(x2,s1),s1,0)==0; subs(x2,s1,1)==0;
                 subs(x2,s2,0)==x1; subs(x2,s2,1)==0];
        S.svars = {{'s1'},{'s1','s2'}};
%       S.expect = 'compatibility';                                         % MMP, 11/14/2026 (was)
        S.expect = 'constrain a state of lower dimension';                  % MMP, 11/14/2026
    case 'S12_domain_mismatch'
        % one variable s1 on two intervals: [0,1] for x1, [0,2] for x2
        x1 = pde_var(1,s1,[0,1]);  x2 = pde_var(1,[s1;s2],[0,2;0,1]);
        S.PDE = [diff(x1,t)==diff(x1,s1,2); diff(x2,t)==diff(x2,s1,2)+diff(x2,s2,2);
                 subs(x1,s1,0)==0; subs(x1,s1,1)==0;
                 subs(x2,s1,0)==0; subs(x2,s1,2)==0;
                 subs(x2,s2,0)==x1; subs(x2,s2,1)==0];
        S.svars = {{'s1'},{'s1','s2'}};
        S.expect = 'domain';
    case 'S13_chain_1d_2d_3d_incompatible'
        % 1-D edge -> 2-D face -> 3-D heat with x3 = 0 on the face s2=0:
        % on the edge s2=s3=0, x3 = 0 and x3 = x2(s1,0) = x1(s1) clash
        x1 = pde_var(1,s1,[0,1]);  x2 = pde_var(1,[s1;s2],[0,1;0,1]);  x3 = pde_var(1,[s1;s2;s3],[0,1;0,1;0,1]);
        S.PDE = [diff(x1,t)==diff(x1,s1,2); diff(x2,t)==diff(x2,s1,2)+diff(x2,s2,2);
                 diff(x3,t)==diff(x3,s1,2)+diff(x3,s2,2)+diff(x3,s3,2);
                 subs(x1,s1,0)==0; subs(x1,s1,1)==0;
                 subs(x2,s1,0)==0; subs(x2,s1,1)==0; subs(x2,s2,0)==x1; subs(x2,s2,1)==0;
                 subs(x3,s1,0)==0; subs(x3,s1,1)==0; subs(x3,s2,0)==0; subs(x3,s2,1)==0;
                 subs(x3,s3,0)==x2; subs(x3,s3,1)==0];
        S.svars = {{'s1'},{'s1','s2'},{'s1','s2','s3'}};
%       S.expect = 'compatibility';                                         % MMP, 11/14/2026 (was)
        S.expect = 'constrain a state of lower dimension';                  % MMP, 11/14/2026
    case 'S13_chain_1d_2d_3d'
        % 1-D edge -> 2-D face -> 3-D heat; x1 also enters the 3-D face
        % s2=0 as (1-s3) x1(s1), so the edge s2=s3=0 is compatible
        x1 = pde_var(1,s1,[0,1]);  x2 = pde_var(1,[s1;s2],[0,1;0,1]);  x3 = pde_var(1,[s1;s2;s3],[0,1;0,1;0,1]);
        S.PDE = [diff(x1,t)==diff(x1,s1,2); diff(x2,t)==diff(x2,s1,2)+diff(x2,s2,2);
                 diff(x3,t)==diff(x3,s1,2)+diff(x3,s2,2)+diff(x3,s3,2);
                 subs(x1,s1,0)==0; subs(x1,s1,1)==0;
                 subs(x2,s1,0)==0; subs(x2,s1,1)==0; subs(x2,s2,0)==x1; subs(x2,s2,1)==0;
                 subs(x3,s1,0)==0; subs(x3,s1,1)==0; subs(x3,s2,0)==(1-s3)*x1; subs(x3,s2,1)==0;
                 subs(x3,s3,0)==x2; subs(x3,s3,1)==0];
        S.svars = {{'s1'},{'s1','s2'},{'s1','s2','s3'}};
        S.xs  = @(c,W) states_S13(c,s1,s2,s3,dd);
        S.rhs = @(Xp,W) {D(Xp{1},s1,2), D(Xp{2},s1,2)+D(Xp{2},s2,2), D(Xp{3},s1,2)+D(Xp{3},s2,2)+D(Xp{3},s3,2)};
        S.bc  = @(Xp,W) {subs(Xp{1},s1,0), subs(Xp{1},s1,1), subs(Xp{2},s1,0), subs(Xp{2},s1,1), ...
                       subs(Xp{2},s2,0)-Xp{1}, subs(Xp{2},s2,1), subs(Xp{3},s1,0), subs(Xp{3},s1,1), ...
                       subs(Xp{3},s2,0)-(1-s3)*Xp{1}, subs(Xp{3},s2,1), subs(Xp{3},s3,0)-Xp{2}, subs(Xp{3},s3,1)};
    case 'S14_input_in_source_group_2d'
        % an input in the s1=1 condition of a 2-D state whose edge s2=0 is a
        % 1-D state (unsupported by the converter, as documented)
        x1 = pde_var(1,s1,[0,1]);  x2 = pde_var(1,[s1;s2],[0,1;0,1]);  w = pde_var('in',1);
        S.PDE = [diff(x1,t)==diff(x1,s1,2); diff(x2,t)==diff(x2,s1,2)+diff(x2,s2,2);
                 subs(x1,s1,0)==0; subs(diff(x1,s1),s1,1)==w;
                 subs(x2,s1,0)==0; subs(diff(x2,s1),s1,1)==(1-s2)*w;
                 subs(x2,s2,0)==x1; subs(x2,s2,1)==0];
        S.svars = {{'s1'},{'s1','s2'}};
        S.expect = 'input';
    case 'S15_ode_input_same_bc'
        % 1-D: x1(0) = X + w, the input lifted into the ODE source
        X = pde_var(1);  x1 = pde_var(1,s1,[0,1]);  w = pde_var('in',1);  z = pde_var('out',1);
        S.PDE = [diff(X,t)==-X; diff(x1,t)==diff(x1,s1,2);
                 z==int(x1,s1,[0,1]);
                 subs(x1,s1,0)==X+w; subs(x1,s1,1)==0];
        S.svars = {{},{'s1'}};  S.wvars = {{}};
        S.xs  = @(c,W) {c(1), (c(1)+W{1})*(1-s1)+dd(s1,c(2:4))};
        S.rhs = @(Xp,W) {-Xp{1}, D(Xp{2},s1,2)};
        S.bc  = @(Xp,W) {subs(Xp{2},s1,0)-Xp{1}-W{1}, subs(Xp{2},s1,1)};
        S.z   = @(Xp,W) {int(Xp{2},s1,0,1)};
    case 'S16_io_indomain'
        % S6 with a disturbance w(s1,s2) in the 2-D equation, a scalar output
        % and an L2[s1] output (the trace of x2 on s2=1)
        g = getp(p,'g',1);
        x1 = pde_var(1,s1,[0,1]);  x2 = pde_var(1,[s1;s2],[0,1;0,1]);
        w = pde_var('in',1,[s1;s2],[0,1;0,1]);  z = pde_var('out',1);  y = pde_var('out',1,s1,[0,1]);
        S.PDE = [diff(x1,t)==diff(x1,s1,2)+g*int(x2,s2,[0,1]);
                 diff(x2,t)==diff(x2,s1,2)+diff(x2,s2,2)+g*x1+w;
                 z==int(x1,s1,[0,1]);  y==subs(x2,s2,1);
                 subs(x1,s1,0)==0; subs(x1,s1,1)==0;
                 subs(x2,s1,0)==0; subs(x2,s1,1)==0;
                 subs(diff(x2,s2),s2,0)==0; subs(diff(x2,s2),s2,1)==0];
        S.svars = {{'s1'},{'s1','s2'}};  S.wvars = {{'s1','s2'}};
        S.xs  = @(c,W) {dd(s1,c(1:3)), dd(s1,c(4:6))*nd(s2,c(7:10))+dd(s1,c(11:13))*nd(s2,c(14:17))};
        S.rhs = @(Xp,W) {D(Xp{1},s1,2)+g*int(Xp{2},s2,0,1), D(Xp{2},s1,2)+D(Xp{2},s2,2)+g*Xp{1}+W{1}};
        S.bc  = @(Xp,W) {subs(Xp{1},s1,0), subs(Xp{1},s1,1), subs(Xp{2},s1,0), subs(Xp{2},s1,1), ...
                       subs(D(Xp{2},s2,1),s2,0), subs(D(Xp{2},s2,1),s2,1)};
        S.z   = @(Xp,W) {int(Xp{1},s1,0,1), subs(Xp{2},s2,1)};
    case 'S17_transport_edge_diffusion'
        % 1-D transport on the edge s2=0 of a state diffusing in s2 only
        % (order 0 in s1); reaction r2 in x2
        r2 = getp(p,'r2',0);
        x1 = pde_var(1,s1,[0,1]);  x2 = pde_var(1,[s1;s2],[0,1;0,1],[0;2]);
        S.PDE = [diff(x1,t)==-diff(x1,s1); diff(x2,t)==diff(x2,s2,2)+r2*x2;
                 subs(x1,s1,0)==0;
                 subs(x2,s2,0)==x1; subs(x2,s2,1)==0];
        S.svars = {{'s1'},{'s1','s2'}};
        S.xs  = @(c,W) {s1*(c(1)+c(2)*s1+c(3)*s1^2), (1-s2)*s1*(c(1)+c(2)*s1+c(3)*s1^2)+s2*(1-s2)*(c(4)+c(5)*s1+c(6)*s1*s2)};
        S.rhs = @(Xp,W) {-D(Xp{1},s1,1), D(Xp{2},s2,2)+r2*Xp{2}};
        S.bc  = @(Xp,W) {subs(Xp{1},s1,0), subs(Xp{2},s2,0)-Xp{1}, subs(Xp{2},s2,1)};
    case 'S18_lib43_interior_delay'
        % library PDE 43: heat with an interior delay, the delay as a 2-D
        % transport state; sgn=+1 is the library file as coded (+a1 term),
        % sgn=-1 the system of its header and of the cited reference
        tau = getp(p,'tau',1);  a0 = getp(p,'a0',1.9);  a1 = getp(p,'sgn',1)*getp(p,'a1',1);  c0 = 1;
        a = 0;  b = pi;
        x1 = pde_var(s1,[a,b]);
        x2 = pde_var([s1;s2],[a,b;0,tau],[2;1]);
        S.PDE = [diff(x1,t)==c0*diff(x1,s1,2)+a0*x1+a1*subs(x2,s2,tau);
                 diff(x2,t)==-diff(x2,s2);
                 subs([x1;x2],s1,a)==0;  subs([x1;x2],s1,b)==0;
                 subs(x2,s2,0)==x1];
        S.svars = {{'s1'},{'s1','s2'}};
        sb = @(c) s1*(pi-s1)*(c(1)+c(2)*s1+c(3)*s1^2);
        S.xs  = @(c,W) {sb(c(1:3)), sb(c(1:3))+s2*s1*(pi-s1)*(c(4)+c(5)*s1+c(6)*s2+c(7)*s1*s2)};
        S.rhs = @(Xp,W) {c0*D(Xp{1},s1,2)+a0*Xp{1}+a1*subs(Xp{2},s2,tau), -D(Xp{2},s2,1)};
        S.bc  = @(Xp,W) {subs(Xp{1},s1,a), subs(Xp{1},s1,b), subs(Xp{2},s1,a), subs(Xp{2},s1,b), subs(Xp{2},s2,0)-Xp{1}};
    case 'S19_2d_integral_in_1d_bc'
        % a HIGHER-dimensional state in the BC of a 1-D state:
        % x1(0) = int int x2, x2 driven in-domain by x1
        x1 = pde_var(1,s1,[0,1]);  x2 = pde_var(1,[s1;s2],[0,1;0,1]);
        S.PDE = [diff(x1,t)==diff(x1,s1,2); diff(x2,t)==diff(x2,s1,2)+diff(x2,s2,2)+x1;
                 subs(x1,s1,0)==int(x2,[s1;s2],[0,1;0,1]); subs(x1,s1,1)==0;
                 subs(x2,s1,0)==0; subs(x2,s1,1)==0; subs(x2,s2,0)==0; subs(x2,s2,1)==0];
        S.svars = {{'s1'},{'s1','s2'}};
        S.xs  = @(c,W) states_S19(c,s1,s2,dd);
        S.rhs = @(Xp,W) {D(Xp{1},s1,2), D(Xp{2},s1,2)+D(Xp{2},s2,2)+Xp{1}};
        S.bc  = @(Xp,W) {subs(Xp{1},s1,0)-int(int(Xp{2},s1,0,1),s2,0,1), subs(Xp{1},s1,1), ...
                       subs(Xp{2},s1,0), subs(Xp{2},s1,1), subs(Xp{2},s2,0), subs(Xp{2},s2,1)};
    case 'S20_reversed_var_order'
        % S4 with the 2-D state declared on [s2;s1]
        x1 = pde_var(1,s1,[0,1]);  x2 = pde_var(1,[s2;s1],[0,1;0,1]);
        S.PDE = [diff(x1,t)==diff(x1,s1,2)+subs(diff(x2,s2),s2,0);
                 diff(x2,t)==diff(x2,s1,2)+diff(x2,s2,2);
                 subs(x1,s1,0)==0; subs(x1,s1,1)==0;
                 subs(x2,s1,0)==0; subs(x2,s1,1)==0;
                 subs(x2,s2,0)==x1; subs(x2,s2,1)==0];
        S.svars = {{'s1'},{'s1','s2'}};
        S.xs  = @(c,W) {dd(s1,c(1:3)), (1-s2)*dd(s1,c(1:3))+s1*(1-s1)*s2*(1-s2)*(c(4)+c(5)*s1+c(6)*s2+c(7)*s1*s2^2)};
        S.rhs = @(Xp,W) {D(Xp{1},s1,2)+subs(D(Xp{2},s2,1),s2,0), D(Xp{2},s1,2)+D(Xp{2},s2,2)};
        S.bc  = @(Xp,W) {subs(Xp{1},s1,0), subs(Xp{1},s1,1), subs(Xp{2},s1,0), subs(Xp{2},s1,1), subs(Xp{2},s2,0)-Xp{1}, subs(Xp{2},s2,1)};
    case 'S21_nonunit_domains'
        % S4 on [-1,2] x [0.5,3], with an s1-dependent coefficient on the
        % edge source, x2(s1,0.5) = (2+s1) x1(s1)
        a = -1;  b = 2;  c1 = 0.5;  d1 = 3;
        x1 = pde_var(1,s1,[a,b]);  x2 = pde_var(1,[s1;s2],[a,b;c1,d1]);
        S.PDE = [diff(x1,t)==diff(x1,s1,2)+subs(diff(x2,s2),s2,c1);
                 diff(x2,t)==diff(x2,s1,2)+0.5*diff(x2,s2,2);
                 subs(x1,s1,a)==0; subs(x1,s1,b)==0;
                 subs(x2,s1,a)==0; subs(x2,s1,b)==0;
                 subs(x2,s2,c1)==(2+s1)*x1; subs(x2,s2,d1)==0];
        S.svars = {{'s1'},{'s1','s2'}};
        q1 = @(c) (s1-a)*(b-s1)*(c(1)+c(2)*s1+c(3)*s1^2);
        S.xs  = @(c,W) {q1(c(1:3)), (d1-s2)/(d1-c1)*(2+s1)*q1(c(1:3))+(s1-a)*(b-s1)*(s2-c1)*(d1-s2)*(c(4)+c(5)*s1+c(6)*s2+c(7)*s1*s2^2)};
        S.rhs = @(Xp,W) {D(Xp{1},s1,2)+subs(D(Xp{2},s2,1),s2,c1), D(Xp{2},s1,2)+0.5*D(Xp{2},s2,2)};
        S.bc  = @(Xp,W) {subs(Xp{1},s1,a), subs(Xp{1},s1,b), subs(Xp{2},s1,a), subs(Xp{2},s1,b), subs(Xp{2},s2,c1)-(2+s1)*Xp{1}, subs(Xp{2},s2,d1)};
    case 'S22_vector_states'
        % a 2-vector 1-D state on the edge of a scalar 2-D state, coupled by
        % matrices: x1_t = x1_ss + Am x1 + [1;0] x2_s2(s1,0), x2(s1,0) = Cm x1
        x1 = pde_var(2,s1,[0,1]);  x2 = pde_var(1,[s1;s2],[0,1;0,1]);
        Am = [0 1;-1 0];  Cm = [1 2];
        S.PDE = [diff(x1,t)==diff(x1,s1,2)+Am*x1+[1;0]*subs(diff(x2,s2),s2,0);
                 diff(x2,t)==diff(x2,s1,2)+diff(x2,s2,2);
                 subs(x1,s1,0)==[0;0]; subs(x1,s1,1)==[0;0];
                 subs(x2,s1,0)==0; subs(x2,s1,1)==0;
                 subs(x2,s2,0)==Cm*x1; subs(x2,s2,1)==0];
        S.svars = {{'s1'},{'s1','s2'}};
        S.xs  = @(c,W) {[dd(s1,c(1:3));dd(s1,c(8:10))], (1-s2)*(dd(s1,c(1:3))+2*dd(s1,c(8:10)))+s1*(1-s1)*s2*(1-s2)*(c(4)+c(5)*s1+c(6)*s2+c(7)*s1*s2^2)};
        S.rhs = @(Xp,W) {D(Xp{1},s1,2)+Am*Xp{1}+[1;0]*subs(D(Xp{2},s2,1),s2,0), D(Xp{2},s1,2)+D(Xp{2},s2,2)};
        S.bc  = @(Xp,W) {subs(Xp{1},s1,0), subs(Xp{1},s1,1), subs(Xp{2},s1,0), subs(Xp{2},s1,1), subs(Xp{2},s2,0)-Cm*Xp{1}, subs(Xp{2},s2,1)};
    case 'S23_hinf_indomain'
        % S6 with a disturbance w(s1,s2) in the 2-D equation, z = int x1
        g = getp(p,'g',1);
        x1 = pde_var(1,s1,[0,1]);  x2 = pde_var(1,[s1;s2],[0,1;0,1]);
        w = pde_var('in',1,[s1;s2],[0,1;0,1]);  z = pde_var('out',1);
        S.PDE = [diff(x1,t)==diff(x1,s1,2)+g*int(x2,s2,[0,1]);
                 diff(x2,t)==diff(x2,s1,2)+diff(x2,s2,2)+g*x1+w;
                 z==int(x1,s1,[0,1]);
                 subs(x1,s1,0)==0; subs(x1,s1,1)==0;
                 subs(x2,s1,0)==0; subs(x2,s1,1)==0;
                 subs(diff(x2,s2),s2,0)==0; subs(diff(x2,s2),s2,1)==0];
        S.svars = {{'s1'},{'s1','s2'}};  S.wvars = {{'s1','s2'}};
        S.xs  = @(c,W) {dd(s1,c(1:3)), dd(s1,c(4:6))*nd(s2,c(7:10))+dd(s1,c(11:13))*nd(s2,c(14:17))};
        S.rhs = @(Xp,W) {D(Xp{1},s1,2)+g*int(Xp{2},s2,0,1), D(Xp{2},s1,2)+D(Xp{2},s2,2)+g*Xp{1}+W{1}};
        S.bc  = @(Xp,W) {subs(Xp{1},s1,0), subs(Xp{1},s1,1), subs(Xp{2},s1,0), subs(Xp{2},s1,1), ...
                       subs(D(Xp{2},s2,1),s2,0), subs(D(Xp{2},s2,1),s2,1)};
        S.z   = @(Xp,W) {int(Xp{1},s1,0,1)};
    case 'S24_hinf_ode_2d_edge'
        % S2 with w into the ODE and z = int int x2
        a = getp(p,'a',1);
        X = pde_var(1);  x2 = pde_var(1,[s1;s2],[0,1;0,1]);
        w = pde_var('in',1);  z = pde_var('out',1);
        S.PDE = [diff(X,t)==-a*X+w; diff(x2,t)==diff(x2,s1,2)+diff(x2,s2,2);
                 z==int(x2,[s1;s2],[0,1;0,1]);
                 subs(diff(x2,s1),s1,0)==0; subs(diff(x2,s1),s1,1)==0;
                 subs(x2,s2,0)==X; subs(x2,s2,1)==0];
        S.svars = {{},{'s1','s2'}};  S.wvars = {{}};
        S.xs  = @(c,W) {c(1), c(1)*(1-s2)+s2*(1-s2)*(nd(s1,c(2:5))*(c(6)+c(7)*s2)+nd(s1,c(8:11))*s2^2)};
        S.rhs = @(Xp,W) {-a*Xp{1}+W{1}, D(Xp{2},s1,2)+D(Xp{2},s2,2)};
        S.bc  = @(Xp,W) {subs(D(Xp{2},s1,1),s1,0), subs(D(Xp{2},s1,1),s1,1), subs(Xp{2},s2,0)-Xp{1}, subs(Xp{2},s2,1)};
        S.z   = @(Xp,W) {int(int(Xp{2},s1,0,1),s2,0,1)};
    case 'S25_hinf_wentzell'
        % S4 with a disturbance w(s1) on the edge state, z = int x1
        c0 = getp(p,'c',1);  r = getp(p,'r',0);
        x1 = pde_var(1,s1,[0,1]);  x2 = pde_var(1,[s1;s2],[0,1;0,1]);
        w = pde_var('in',1,s1,[0,1]);  z = pde_var('out',1);
        S.PDE = [diff(x1,t)==diff(x1,s1,2)+c0*subs(diff(x2,s2),s2,0)+r*x1+w;
                 diff(x2,t)==diff(x2,s1,2)+diff(x2,s2,2)+r*x2;
                 z==int(x1,s1,[0,1]);
                 subs(x1,s1,0)==0; subs(x1,s1,1)==0;
                 subs(x2,s1,0)==0; subs(x2,s1,1)==0;
                 subs(x2,s2,0)==x1; subs(x2,s2,1)==0];
        S.svars = {{'s1'},{'s1','s2'}};  S.wvars = {{'s1'}};
        S.xs  = @(c,W) {dd(s1,c(1:3)), (1-s2)*dd(s1,c(1:3))+s1*(1-s1)*s2*(1-s2)*(c(4)+c(5)*s1+c(6)*s2+c(7)*s1*s2^2)};
        S.rhs = @(Xp,W) {D(Xp{1},s1,2)+c0*subs(D(Xp{2},s2,1),s2,0)+r*Xp{1}+W{1}, D(Xp{2},s1,2)+D(Xp{2},s2,2)+r*Xp{2}};
        S.bc  = @(Xp,W) {subs(Xp{1},s1,0), subs(Xp{1},s1,1), subs(Xp{2},s1,0), subs(Xp{2},s1,1), subs(Xp{2},s2,0)-Xp{1}, subs(Xp{2},s2,1)};
        S.z   = @(Xp,W) {int(Xp{1},s1,0,1)};
    case 'S27_derivative_source'
        % the derivative of a 1-D state as Dirichlet data on the edge s2=0 of
        % a 2-D state: x2 in H^{2,2} needs d^3 x1/ds1^3, beyond the order 2
        % of x1, so the converter must refuse it
        S.expect = 'differentiable';
        x1 = pde_var(1,s1,[0,1]);  x2 = pde_var(1,[s1;s2],[0,1;0,1]);
        S.PDE = [diff(x1,t)==diff(x1,s1,2)-x1; diff(x2,t)==diff(x2,s1,2)+diff(x2,s2,2);
                 subs(diff(x1,s1),s1,0)==0; subs(diff(x1,s1),s1,1)==0;
                 subs(x2,s1,0)==0; subs(x2,s1,1)==0;
                 subs(x2,s2,0)==diff(x1,s1); subs(x2,s2,1)==0];
        S.svars = {{'s1'},{'s1','s2'}};
        S.xs  = @(c,W) {nd(s1,c(1:4)), (1-s2)*diff(nd(s1,c(1:4)),s1)+s1*(1-s1)*s2*(1-s2)*(c(5)+c(6)*s1+c(7)*s2)};
        S.rhs = @(Xp,W) {D(Xp{1},s1,2)-Xp{1}, D(Xp{2},s1,2)+D(Xp{2},s2,2)};
        S.bc  = @(Xp,W) {subs(D(Xp{1},s1,1),s1,0), subs(D(Xp{1},s1,1),s1,1), subs(Xp{2},s1,0), subs(Xp{2},s1,1), ...
                         subs(Xp{2},s2,0)-D(Xp{1},s1,1), subs(Xp{2},s2,1)};
    case 'S28_robin_edge_source'
        % Robin data from a 1-D state on the edge s2=0:
        % x2_s2(s1,0) - x2(s1,0) = -x1(s1)
        x1 = pde_var(1,s1,[0,1]);  x2 = pde_var(1,[s1;s2],[0,1;0,1]);
        S.PDE = [diff(x1,t)==diff(x1,s1,2); diff(x2,t)==diff(x2,s1,2)+diff(x2,s2,2);
                 subs(x1,s1,0)==0; subs(x1,s1,1)==0;
                 subs(x2,s1,0)==0; subs(x2,s1,1)==0;
                 subs(diff(x2,s2),s2,0)-subs(x2,s2,0)==-x1; subs(x2,s2,1)==0];
        S.svars = {{'s1'},{'s1','s2'}};
        S.xs  = @(c,W) {dd(s1,c(1:3)), dd(s1,c(1:3))*(1-s2)/2+s1*(1-s1)*s2^2*(1-s2)*(c(4)+c(5)*s1+c(6)*s2)};
        S.rhs = @(Xp,W) {D(Xp{1},s1,2), D(Xp{2},s1,2)+D(Xp{2},s2,2)};
        S.bc  = @(Xp,W) {subs(Xp{1},s1,0), subs(Xp{1},s1,1), subs(Xp{2},s1,0), subs(Xp{2},s1,1), ...
                         subs(D(Xp{2},s2,1),s2,0)-subs(Xp{2},s2,0)+Xp{1}, subs(Xp{2},s2,1)};
    case 'S29_piesim2d_ex8'
        % PIESIM 2-D library example 8: ODE + L2[x] + L2[y] + L2[x,y] coupled
        % through integral BCs, on [0,1] x [0,3], visc 0.5
        a = 0; b = 1; c1 = 0; d1 = 3; visc = 0.5;
        x1 = pde_var();  x2 = pde_var(s1,[a,b]);  x3 = pde_var(s2,[c1,d1]);  x4 = pde_var([s1;s2],[a,b;c1,d1]);
        S.PDE = [diff(x1,t)==0;
                 diff(x2,t)==visc*diff(x2,s1,2);
                 diff(x3,t)==visc*diff(x3,s2,2);
                 diff(x4,t)==visc*(diff(x4,s1,2)+diff(x4,s2,2));
                 int(x2,s1,a,b)==x1;  subs(diff(x2,s1),s1,b)==0;
                 int(x3,s2,c1,d1)==x1;  subs(diff(x3,s2),s2,d1)==0;
                 int(x4,s1,a,b)==x3;  int(x4,s2,c1,d1)==x2;
                 subs(diff(x4,s1),s1,b)==0;   subs(diff(x4,s2),s2,d1)==0];
        S.svars = {{},{'s1'},{'s2'},{'s1','s2'}};
        S.xs  = @(c,W) states_S29(c,s1,s2);
        S.rhs = @(Xp,W) {0*Xp{1}, visc*D(Xp{2},s1,2), visc*D(Xp{3},s2,2), visc*(D(Xp{4},s1,2)+D(Xp{4},s2,2))};
        S.bc  = @(Xp,W) {int(Xp{2},s1,a,b)-Xp{1}, subs(D(Xp{2},s1,1),s1,b), int(Xp{3},s2,c1,d1)-Xp{1}, ...
                         subs(D(Xp{3},s2,1),s2,d1), int(Xp{4},s1,a,b)-Xp{3}, int(Xp{4},s2,c1,d1)-Xp{2}, ...
                         subs(D(Xp{4},s1,1),s1,b), subs(D(Xp{4},s2,1),s2,d1)};
    case 'S30_piesim2d_ex9'
        % PIESIM 2-D library example 9: L2[y] + L2[x,y], an integral and a
        % periodic condition in x, on [-1,1] x [0,1], visc 0.1
        a = -1; b = 1; c1 = 0; d1 = 1; visc = 0.1;
        x1 = pde_var(s2,[c1,d1]);  x2 = pde_var([s1;s2],[a,b;c1,d1]);
        S.PDE = [diff(x1,t)==visc*diff(x1,s2,2);
                 diff(x2,t)==visc*(diff(x2,s1,2)+diff(x2,s2,2));
                 subs(x1,s2,c1)==0;    subs(diff(x1,s2),s2,d1)==0;
                 int(x2,s1,[a,b])==x1;    subs(x2,s1,b)==subs(x2,s1,a);
                 subs(x2,s2,c1)==0;        subs(diff(x2,s2),s2,d1)==0];
        S.svars = {{'s2'},{'s1','s2'}};
        S.xs  = @(c,W) states_S30(c,s1,s2);
        S.rhs = @(Xp,W) {visc*D(Xp{1},s2,2), visc*(D(Xp{2},s1,2)+D(Xp{2},s2,2))};
        S.bc  = @(Xp,W) {subs(Xp{1},s2,c1), subs(D(Xp{1},s2,1),s2,d1), int(Xp{2},s1,a,b)-Xp{1}, ...
                         subs(Xp{2},s1,b)-subs(Xp{2},s1,a), subs(Xp{2},s2,c1), subs(D(Xp{2},s2,1),s2,d1)};
    case 'S31_integral_bc_2d'
        % one 2-D state with an integral condition in s1: int x ds1 = 0,
        % x_s1(1,s2) = 0, Dirichlet in s2 (reproduces PIESIM library ex 8, 9)
        x = pde_var(1,[s1;s2],[0,1;0,1]);
        S.PDE = [diff(x,t)==diff(x,s1,2)+diff(x,s2,2);
                 int(x,s1,[0,1])==0; subs(diff(x,s1),s1,1)==0;
                 subs(x,s2,0)==0; subs(x,s2,1)==0];
        S.svars = {{'s1','s2'}};
        % p'(1) = 0 and int p = 0: p = (s1^2 (3 - 2 s1) - 1/2) c1 + ...
        S.xs  = @(c,W) {s2*(1-s2)*((s1^2*(3-2*s1)-1/2)*(c(1)+c(2)*s2)+(s1*(2-s1)-2/3)*(c(3)+c(4)*s2^2))};
        S.rhs = @(Xp,W) {D(Xp{1},s1,2)+D(Xp{1},s2,2)};
        S.bc  = @(Xp,W) {int(Xp{1},s1,0,1), subs(D(Xp{1},s1,1),s1,1), subs(Xp{1},s2,0), subs(Xp{1},s2,1)};
    case 'C1_heat2d_scalar_w'
        % control for the H-infinity and H2 executives: one 2-D heat state,
        % Dirichlet on every edge, a scalar w as a uniform source, z = int int
        % x; no lower-dimensional and no auxiliary state
        x = pde_var(1,[s1;s2],[0,1;0,1]);  w = pde_var('in',1);  z = pde_var('out',1);
        S.PDE = [diff(x,t)==diff(x,s1,2)+diff(x,s2,2)+w;
                 z==int(x,[s1;s2],[0,1;0,1]);
                 subs(x,s1,0)==0; subs(x,s1,1)==0; subs(x,s2,0)==0; subs(x,s2,1)==0];
        S.svars = {{'s1','s2'}};  S.wvars = {{}};
        S.xs  = @(c,W) {dd(s1,c(1:3))*dd(s2,c(4:6))+s2*dd(s1,c(7:9))*dd(s2,c(10:12))};
        S.rhs = @(Xp,W) {D(Xp{1},s1,2)+D(Xp{1},s2,2)+W{1}};
        S.bc  = @(Xp,W) {subs(Xp{1},s1,0), subs(Xp{1},s1,1), subs(Xp{1},s2,0), subs(Xp{1},s2,1)};
        S.z   = @(Xp,W) {int(int(Xp{1},s1,0,1),s2,0,1)};
    case 'S26_lib42_ode_heat_delay'
        % library PDE 42 (ODE + heat, both with a delay, as transport states
        % x3(s2), x4(s1,s2) on s2 in [1,1+tau]); reaction a in x2
        tau = 1;  A0 = -1;  A1 = 0.5;  B0 = -1;  a = getp(p,'a',1);  a2 = 0.5;
        x1 = pde_var(1);  x2 = pde_var(1,s1,[0,1]);  x3 = pde_var(1,s2,[1,1+tau]);
        x4 = pde_var(1,[s1;s2],[0,1;1,1+tau]);
        S.PDE = [diff(x1,t)==A0*x1+B0*subs(x2,s1,0)+A1*subs(x3,s2,1);
                 diff(x2,t)==diff(x2,s1,2)+a*x2+a2*subs(x4,s2,1);
                 diff(x3,t)==diff(x3,s2);
                 diff(x4,t)==diff(x4,s2);
                 subs(diff(x2,s1),s1,0)==0;     subs(x2,s1,1)==0;
                 subs(x3,s2,1+tau)==x1;
                 subs(diff(x4,s1),s1,0)==0;     subs(x4,s1,1)==0;
                 subs(x4,s2,1+tau)==x2];
        S.svars = {{},{'s1'},{'s2'},{'s1','s2'}};
        S.xs  = @(c,W) states_S26(c,s1,s2,nD,tau);
        S.rhs = @(Xp,W) {A0*Xp{1}+B0*subs(Xp{2},s1,0)+A1*subs(Xp{3},s2,1), D(Xp{2},s1,2)+a*Xp{2}+a2*subs(Xp{4},s2,1), ...
                         D(Xp{3},s2,1), D(Xp{4},s2,1)};
        S.bc  = @(Xp,W) {subs(D(Xp{2},s1,1),s1,0), subs(Xp{2},s1,1), subs(Xp{3},s2,1+tau)-Xp{1}, ...
                         subs(D(Xp{4},s1,1),s1,0), subs(Xp{4},s1,1), subs(Xp{4},s2,1+tau)-Xp{2}};
    otherwise
        error('mixed_domain_systems:name','Unknown system %s.',name);
end
end



function C = catalogue()
% C = CATALOGUE() returns the table of systems.
%
% OUTPUT
% - C:  struct array with fields name, spaces (the PIE state spaces), what
%       (the coupling exercised) and baseline (the reference used).
R = 'R';  L1 = 'L2[s1]';  L2 = 'L2[s2]';  L12 = 'L2[s1,s2]';
row = @(n,sp,wh,bl) struct('name',n,'spaces',sp,'what',wh,'baseline',bl);
C = [row('S1_ode_heat1d_bc',[R,' x ',L1],'ODE in a Dirichlet BC of a 1-D heat','spectrum {-a} u {r - n^2 pi^2}');
     row('S2_ode_heat2d_edge',[R,' x ',L2,' x ',L12],'ODE as data on an edge of a 2-D heat; Neumann direction (auxiliary state, K)','spectrum {-a} u {-(m^2+n^2) pi^2, n>=1}');
     row('S3_edge1d_heat2d',[L1,' x ',L12],'1-D state as Dirichlet data on an edge (cascade)','margin r* = pi^2');
     row('S4_wentzell',[L1,' x ',L12],'dynamic (Wentzell) edge, two-way 1-D <-> 2-D','margin r* = pi^2 + k0^2, k0 tan k0 = c');
     row('S5_chain_ode_1d_2d',[R,' x ',L1,' x ',L12],'ODE -> 1-D (Neumann data) -> 2-D edge, ODE in a 2-D BC with an s2 coefficient','spectrum {-a} u ... (cascade)');
     row('S6_indomain_1d_2d',[L1,' x ',L1,' x ',L12],'in-domain 1-D <-> 2-D (integral, multiplier); Neumann in s2 (aux, K)','margin g* = pi^2');
     row('S7_ode_2d_indomain',[R,' x ',R,' x ',L1,' x ',L2,' x ',L12],'ODE <-> 2-D in-domain; Neumann on every edge (3 aux, K)','margin g* = sqrt(a b)');
     row('S8_two_edges',[L1,' x ',L2,' x ',L12],'1-D states on two edges (lower s2=0, upper s1=1)','margin r3* = pi^2');
     row('S9_ode_2d_flux',[R,' x ',L2,' x ',L12],'ODE fixes an edge value, the edge flux drives the ODE (two-way, K)','margin a* = c; H-inf from 1/(s - a + c sqrt(s) coth sqrt(s))');
     row('S10_incompatible_corner',[R,' x ',L12],'corner data clash','must be rejected');
     row('S11_incompatible_neumann_corner',[L1,' x ',L12],'derivative corner clash','must be rejected');
     row('S12_domain_mismatch',[L1,' x ',L12],'one variable on two intervals','must be rejected');
     row('S13_chain_1d_2d_3d_incompatible',[L1,' x ',L12,' x L2[s1,s2,s3]'],'edge clash in 3-D','must be rejected');
     row('S13_chain_1d_2d_3d',[L1,' x ',L12,' x L2[s1,s2,s3]'],'1-D edge -> 2-D face -> 3-D','semantics only');
     row('S14_input_in_source_group_2d',[L1,' x ',L12],'input in a 2-D BC group with a lower-dimensional source','must be rejected (documented)');
     row('S15_ode_input_same_bc',[R,' x ',L1],'input and ODE in one 1-D BC (input lifted into the source)','semantics only');
     row('S16_io_indomain',[L1,' x ',L1,' x ',L12],'L2[s1,s2] input, scalar and L2[s1] (trace) outputs','semantics only');
     row('S17_transport_edge_diffusion',[L1,' x ',L12],'transport edge state, target of order 0 in s1','margin r2* = pi^2');
     row('S18_lib43_interior_delay',[L1,' x ',L12],'library PDE 43 (delay as 2-D transport), both signs','tau* = 1.034726 (sgn -1; library header: 1.0347, citing Caliskan 2009); unstable for every tau (sgn +1)');
%    row('S19_2d_integral_in_1d_bc',[L1,' x ',L12],'higher-dimensional state in a 1-D BC','PIE exists; converter rejects'); % MMP, 11/14/2026 (was)
     row('S19_2d_integral_in_1d_bc',[L1,' x ',L12],'higher-dimensional state in a 1-D BC','PIE exists; converted since dbf458a2'); % MMP, 11/14/2026
     row('S20_reversed_var_order',[L1,' x ',L12],'2-D state declared on [s2;s1]','semantics only');
     row('S21_nonunit_domains',[L1,' x ',L12],'domains [-1,2] x [0.5,3], s1-dependent edge coefficient','semantics only');
     row('S22_vector_states',['L2^2[s1] x ',L12],'vector-valued edge state, matrix coupling','semantics only');
     row('S23_hinf_indomain',[L1,' x ',L1,' x ',L12],'H-inf, L2[s1,s2] disturbance','gain by series over s1 modes');
     row('S24_hinf_ode_2d_edge',[R,' x ',L2,' x ',L12],'H-inf, ODE -> 2-D edge (K)','gain 1/(2a)');
     row('S25_hinf_wentzell',[L1,' x ',L12],'H-inf, L2[s1] disturbance on the Wentzell edge','gain by series over s1 modes');
     row('S26_lib42_ode_heat_delay',[R,' x ',L1,' x ',L2,' x ',L12],'library PDE 42 (Kang & Fridman 2017): ODE + heat with delays as transport states','margin a* = pi^2/4 - 0.5 (mode-0 delay equation)');
     row('S27_derivative_source',[R,' x ',L1,' x ',L12],'derivative of a 1-D state as edge data (needs d^3 x1/ds1^3)','must be rejected');
     row('S28_robin_edge_source',[L1,' x ',L12],'Robin edge condition with 1-D data','semantics only');
     row('S29_piesim2d_ex8',[R,' x ',L1,' x ',L2,' x ',L12],'PIESIM 2-D library ex 8: integral BCs coupling all four spaces','semantics; PIESIM vs exact solution');
     row('S30_piesim2d_ex9',[L2,' x ',L12],'PIESIM 2-D library ex 9: integral and periodic BCs in x','semantics; PIESIM vs exact solution');
     row('S31_integral_bc_2d',L12,'one 2-D state with an integral BC in s1 (PIESIM reproduction)','semantics; PIESIM vs exact solution');
     row('C1_heat2d_scalar_w',L12,'control: one 2-D heat state, Dirichlet on every edge, scalar w as a uniform source, z = int int x (no lower-dimensional or auxiliary state)','gain 0.0351442537, H2 norm 0.1150629240 (series over the modes)')];
end



function v = getp(p,f,d)
% V = GETP(P,F,D) returns the field F of the 'struct' P, or D if absent.
if isfield(p,f), v = p.(f); else, v = d; end
end

function q = dk(q,v,k)
% Q = DK(Q,V,K) differentiates the 'polynomial' Q K times in the variable V.
for i=1:k
    q = diff(q,v);
end
end

function X = states_S5(c,s1,s2,nD)
% X = STATES_S5(C,S1,S2,ND) returns states of S5 satisfying its BCs:
% x1'(0) = X, x1(1) = 0; x2(s1,0) = x1, x2_s1(0,s2) = (1-s2) X, x2 = 0 on
% s1=1 and s2=1.
Xo = c(1);
x1 = Xo*(s1-1)+nD(s1,c(2:4));
g  = nD(s1,c(5:7))*(c(8)+c(9)*s2)+nD(s1,c(10:12))*s2^2;
x2 = (1-s2)*x1+s2*(1-s2)*g;
X = {Xo,x1,x2};
end

function X = states_S19(c,s1,s2,dd)
% X = STATES_S19(C,S1,S2,DD) returns states of S19 satisfying its BCs:
% x2 Dirichlet on the square, x1(1) = 0, x1(0) = int int x2.
x2 = s1*(1-s1)*s2*(1-s2)*(c(4)+c(5)*s1+c(6)*s2^2);
m  = int(int(x2,s1,0,1),s2,0,1);
x1 = m*(1-s1)+dd(s1,c(1:3));
X = {x1,x2};
end

function X = states_S26(c,s1,s2,nD,tau)
% X = STATES_S26(C,S1,S2,ND,TAU) returns states of S26 satisfying its BCs:
% x2'(0) = 0, x2(1) = 0, x3(1+tau) = x1, x4 = x2 on s2 = 1+tau and of the
% form of x2 in s1.
x1 = c(1);
x2 = nD(s1,c(2:4));
x3 = x1+(s2-1-tau)*(c(5)+c(6)*s2);
x4 = x2+(s2-1-tau)*(nD(s1,c(7:9))+s2*nD(s1,c(10:12)));
X = {x1,x2,x3,x4};
end

function X = states_S29(c,s1,s2)
% X = STATES_S29(C,S1,S2) returns states of S29 satisfying its BCs: u with
% u_x(1,y) = 0, u_y(x,3) = 0, and its moments U1 = int u dy,
% U2 = int u dx, U0 = int int u.
p1 = c(1)+c(2)*s1*(2-s1);  p2 = c(3)+c(4)*s1^2*(3-2*s1);
q1 = c(5)+c(6)*s2*(6-s2);  q2 = c(7)+c(8)*s2^2*(9-2*s2);
u  = p1*q1+p2*q2;
U1 = int(u,s2,0,3);  U2 = int(u,s1,0,1);  U0 = int(U1,s1,0,1);
X = {U0,U1,U2,u};
end

function X = states_S30(c,s1,s2)
% X = STATES_S30(C,S1,S2) returns states of S30 satisfying its BCs: u with
% u(1,y) = u(-1,y), u(x,0) = 0, u_y(x,1) = 0, and U1 = int u dx.
p1 = c(1)+c(2)*s1^2+c(3)*(s1^3-s1);  p2 = c(4)+c(5)*s1*(s1^2-1);
q1 = s2*(2-s2);  q2 = s2^2*(3-2*s2);
u  = p1*q1+p2*q2;
U1 = int(u,s1,-1,1);
X = {U1,u};
end

function X = states_S13(c,s1,s2,s3,dd)
% X = STATES_S13(C,S1,S2,S3,DD) returns states of S13 satisfying its BCs.
x1 = dd(s1,c(1:3));
x2 = (1-s2)*x1+s1*(1-s1)*s2*(1-s2)*(c(4)+c(5)*s1*s2);
x3 = (1-s3)*x2+s1*(1-s1)*s2*(1-s2)*s3*(1-s3)*(c(6)+c(7)*s1*s2*s3);
X = {x1,x2,x3};
end
