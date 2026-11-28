function S = high_dim_systems(name)
% S = HIGH_DIM_SYSTEMS(NAME) returns one test system on 4 or 5 spatial
% variables: its PDE, and its semantics written by hand with polynomial
% operations, independently of the parser. S = HIGH_DIM_SYSTEMS('list')
% returns the names, as a cellstr, and their dimensions.
%
% INPUT
% - name:   char, the name of a system, or 'list';
%
% OUTPUT
% - S:      for a system, 'struct' with fields
%             name     char, NAME;
%             dim      4 or 5, the largest number of variables of a state;
%             PDE      'pde_struct', the PDE (not initialized);
%             xid, wid, uid, zid, yid
%                      arrays, the IDs of the declared states, inputs and
%                      outputs, in the order of the cells below;
%             svars    1 x nx cell, svars{k} the cellstr of the variables of
%                      state k;
%             rhs      @(X,W,U): 1 x nx cell, the right-hand side of each
%                      state equation at the 'polynomial' states X and inputs
%                      W (disturbances), U (controls); [] for a PDE without
%                      dynamics;
%             bc       @(X,W,U): cell of 'polynomial' residuals that vanish
%                      iff X satisfies the BCs;
%             z, y     @(X,W,U): cell, the outputs ({} if none);
%             wgen, ugen
%                      @(c): cell of inputs compatible with the BCs, from a
%                      vector c of 20 coefficients;
%             xs       @(c): 1 x nx cell of states satisfying the BCs, for a
%                      PIE with constraints K ([] otherwise);
%             expect   char, '' for a system the converter must accept, else
%                      a fragment of the error it must raise;
%             check    false to check only that the PDE converts;
%           for 'list', a 'struct' array with fields name and dim.
%
% NOTES
% Domains are [0,1] unless stated. Inputs in the BCs vanish on the edges
% shared with Dirichlet faces, as the BCs require.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PIETOOLS - high_dim_systems
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
% For support, contact M. Peet, Arizona State University at mpeet@asu.edu
% Initial coding DJ - 11/25/2026

names = {'H4_dir',4; 'H4_mixed',4; 'T4_transport',4; 'V4_vector',4; 'V4_noncommute',4; 'O4_orders',4;
         'VO4_varorder',4; 'M4_face',4; 'M4_product',4; 'W4_inputs',4; 'W4_input_coef',4;
         'P4_periodic',4; 'P4_periodic2',4; 'E4_corner',4; 'BVP4',4; 'D3_delay',4;
         'E4_tangential',4; 'E4_fewBCs',4;
         'H5_dir',5; 'H5_coupled',5; 'M5_face',5; 'BVP5',5; 'E5_trivial',5};
if strcmp(name,'list')
    S = struct('name',names(:,1),'dim',names(:,2));
    return
end
k = find(strcmp(names(:,1),name));
if isempty(k), error('Unknown system %s.',name); end
pvar t s1 s2 s3 s4 s5 th
V = [s1;s2;s3;s4;s5];   nm = {'s1','s2','s3','s4','s5'};
D = @(q,v,k) dk(q,v,k);
sb = @(q,v,a) subs(polynomial(q),v,a);
I = @(q,v,a,b) int(polynomial(q),v,a,b);
dom01 = @(n) repmat([0,1],n,1);
S = struct('name',name,'dim',names{k,2},'PDE',[],'xid',[],'wid',[],'uid',[],'zid',[],'yid',[], ...
    'svars',{{}},'rhs',[],'bc',[],'z',@(X,W,U) {},'y',@(X,W,U) {}, ...
    'wgen',@(c) {},'ugen',@(c) {},'xs',[],'expect','','check',true);
switch name
    case {'H4_dir','H5_dir'}
        % Heat equation, Dirichlet on every face
        N = S.dim;  vv = V(1:N);
        x = pde_var(vv,dom01(N),2*ones(N,1));
        S.PDE = [diff(x,t)==lap(x,vv); dirichlet(x,vv)];
        S.xid = id(x);  S.svars = {nm(1:N)};
        S.rhs = @(X,W,U) {plap(X{1},vv)};
        S.bc  = @(X,W,U) pdir(X{1},vv);

    case 'H4_mixed'
        % Domains other than [0,1]; Dirichlet, Neumann and Robin BCs;
        % variable advection, reaction, Volterra integrals (one with a kernel
        % in the variable of integration) and a full integral; a scalar
        % output and an output on the face s1=1
        dom = [0,1; -1,1; 0,2; 1,3];
        x = pde_var(V(1:4),dom,[2;2;2;2]);
        z1 = pde_var('out',1);  z2 = pde_var('out',1,[s2;s3;s4],dom(2:4,:));
        P = diff(x,s1,2)+2*diff(x,s2,2)+0.5*diff(x,s3,2)+diff(x,s4,2)+s1*diff(x,s2)+(1+s4)*x;
        P = P+int(x,s3,[0,s3])+int(s3*x,s3,[s3,2])+int(x,s4,[1,3]);
        S.PDE = [diff(x,t)==P;  z1==int(x,V(1:4),dom);  z2==subs(x,s1,1);
                 subs(x,s1,0)==0; subs(diff(x,s1),s1,1)==0;
                 subs(diff(x,s2),s2,-1)-subs(x,s2,-1)==0; subs(x,s2,1)==0;
                 subs(x,s3,0)==0; subs(diff(x,s3),s3,2)+2*subs(x,s3,2)==0;
                 subs(diff(x,s4),s4,1)==0; subs(x,s4,3)==0];
        S.xid = id(x);  S.zid = [id(z1),id(z2)];  S.svars = {nm(1:4)};
        rhs = @(X) D(X,s1,2)+2*D(X,s2,2)+0.5*D(X,s3,2)+D(X,s4,2)+s1*D(X,s2,1)+(1+s4)*X + ...
                   I(sb(X,s3,th),th,0,s3)+I(th*sb(X,s3,th),th,s3,2)+I(X,s4,1,3);
        S.rhs = @(X,W,U) {rhs(X{1})};
        S.z   = @(X,W,U) {I(I(I(I(X{1},s1,0,1),s2,-1,1),s3,0,2),s4,1,3), sb(X{1},s1,1)};
        S.bc  = @(X,W,U) {sb(X{1},s1,0), sb(D(X{1},s1,1),s1,1), sb(D(X{1},s2,1),s2,-1)-sb(X{1},s2,-1), ...
                    sb(X{1},s2,1), sb(X{1},s3,0), sb(D(X{1},s3,1),s3,2)+2*sb(X{1},s3,2), ...
                    sb(D(X{1},s4,1),s4,1), sb(X{1},s4,3)};

    case 'T4_transport'
        % First order in every direction, inflow BCs, variable speed
        x = pde_var(V(1:4),dom01(4),ones(4,1));
        S.PDE = [diff(x,t)==-(1+s2)*diff(x,s1)-2*diff(x,s2)-0.5*diff(x,s3)-diff(x,s4)+x;
                 subs(x,s1,0)==0; subs(x,s2,0)==0; subs(x,s3,0)==0; subs(x,s4,0)==0];
        S.xid = id(x);  S.svars = {nm(1:4)};
        S.rhs = @(X,W,U) {-(1+s2)*D(X{1},s1,1)-2*D(X{1},s2,1)-0.5*D(X{1},s3,1)-D(X{1},s4,1)+X{1}};
        S.bc  = @(X,W,U) {sb(X{1},s1,0), sb(X{1},s2,0), sb(X{1},s3,0), sb(X{1},s4,0)};

    case {'V4_vector','V4_noncommute'}
        % 2-vector state of orders (2,2,1,0), matrix coefficients and BCs.
        % The BC matrices along s1 and s3 commute in V4_vector; in
        % V4_noncommute they do not, so that the BCs imply
        % [R1,P3]*x(1,s2,1,s4) = 0, and the converter refuses
        x = pde_var(2,V(1:4),dom01(4),[2;2;1;0]);
        A1 = [1 0;0 2];  B3 = [1 1;0 1];  C0 = [0 1;-1 0];  R1 = [2 1;1 2];  P3 = 0.5*[0 1;1 0];
        if strcmp(name,'V4_noncommute'), R1 = [1 2;0 1];  S.expect = 'do not commute'; end
        S.PDE = [diff(x,t)==A1*diff(x,s1,2)+diff(x,s2,2)+B3*diff(x,s3)+C0*x+int(x,s4,[0,1]);
                 subs(x,s1,0)==0; R1*subs(x,s1,1)+subs(diff(x,s1),s1,1)==0;
                 subs(x,s2,0)==0; subs(x,s2,1)==0;
                 subs(x,s3,0)-P3*subs(x,s3,1)==0];
        S.xid = id(x);  S.svars = {nm(1:4)};
        S.rhs = @(X,W,U) {A1*D(X{1},s1,2)+D(X{1},s2,2)+B3*D(X{1},s3,1)+C0*X{1}+I(X{1},s4,0,1)};
        S.bc  = @(X,W,U) {sb(X{1},s1,0), R1*sb(X{1},s1,1)+sb(D(X{1},s1,1),s1,1), sb(X{1},s2,0), sb(X{1},s2,1), ...
                    sb(X{1},s3,0)-P3*sb(X{1},s3,1)};

    case 'O4_orders'
        % Orders (2,1,0,3): third order in s4, no derivative in s3
        x = pde_var(V(1:4),dom01(4),[2;1;0;3]);
        S.PDE = [diff(x,t)==diff(x,s1,2)-(1+s3)*diff(x,s2)+diff(x,s4,3)+s3*x;
                 subs(x,s1,0)==0; subs(x,s1,1)==0; subs(x,s2,0)==0;
                 subs(x,s4,0)==0; subs(x,s4,1)==0; subs(diff(x,s4),s4,0)==0];
        S.xid = id(x);  S.svars = {nm(1:4)};
        S.rhs = @(X,W,U) {D(X{1},s1,2)-(1+s3)*D(X{1},s2,1)+D(X{1},s4,3)+s3*X{1}};
        S.bc  = @(X,W,U) {sb(X{1},s1,0), sb(X{1},s1,1), sb(X{1},s2,0), sb(X{1},s4,0), sb(X{1},s4,1), sb(D(X{1},s4,1),s4,0)};

    case 'VO4_varorder'
        % Variables declared out of order, with an order per declared
        % variable and different domains; a second state y(s4,s3)
        x = pde_var([s3;s1;s4;s2],[0,2; 0,1; -1,1; 0,1],[2;1;2;1]);
        y = pde_var([s4;s3],[-1,1; 0,2],[2;2]);
        S.PDE = [diff(x,t)==diff(x,s3,2)-diff(x,s1)+diff(x,s4,2)-2*diff(x,s2)+y;
                 diff(y,t)==diff(y,s4,2)+diff(y,s3,2)+int(x,[s1;s2],dom01(2));
                 subs(x,s3,0)==0; subs(x,s3,2)==0; subs(x,s1,0)==0;
                 subs(x,s4,-1)==0; subs(diff(x,s4),s4,1)==0; subs(x,s2,0)==0;
                 subs(y,s3,0)==0; subs(y,s3,2)==0; subs(y,s4,-1)==0; subs(y,s4,1)==0];
        S.xid = [id(y),id(x)];  S.svars = {{'s3','s4'},nm(1:4)};
        S.rhs = @(X,W,U) {D(X{1},s4,2)+D(X{1},s3,2)+I(I(X{2},s1,0,1),s2,0,1), ...
                          D(X{2},s3,2)-D(X{2},s1,1)+D(X{2},s4,2)-2*D(X{2},s2,1)+X{1}};
        S.bc  = @(X,W,U) {sb(X{2},s3,0), sb(X{2},s3,2), sb(X{2},s1,0), sb(X{2},s4,-1), sb(D(X{2},s4,1),s4,1), ...
                    sb(X{2},s2,0), sb(X{1},s3,0), sb(X{1},s3,2), sb(X{1},s4,-1), sb(X{1},s4,1)};

    case 'M4_face'
        % An ODE, 2D, 3D and 4D state; the 3D state is the value of the 4D
        % state on the face s4=0; distributed couplings both ways
        X0 = pde_var(1);  y = pde_var([s1;s2],dom01(2),[2;2]);
        z3 = pde_var(V(1:3),dom01(3),[2;2;2]);  x = pde_var(V(1:4),dom01(4),2*ones(4,1));
        S.PDE = [diff(X0,t)==-X0+int(x,V(1:4),dom01(4));
                 diff(y,t)==lap(y,V(1:2))+int(x,[s3;s4],dom01(2))+X0;
                 diff(z3,t)==lap(z3,V(1:3))+int(x,s4,[0,1])+y;
                 diff(x,t)==lap(x,V(1:4))+z3+y+s3*X0;
                 dirichlet(y,V(1:2)); dirichlet(z3,V(1:3)); dirichlet(x,V(1:3));
                 subs(x,s4,0)==z3; subs(x,s4,1)==0];
        S.xid = [id(X0),id(y),id(z3),id(x)];  S.svars = {{},nm(1:2),nm(1:3),nm(1:4)};
        S.rhs = @(X,W,U) {-X{1}+I(I(I(I(X{4},s1,0,1),s2,0,1),s3,0,1),s4,0,1), ...
                          plap(X{2},V(1:2))+I(I(X{4},s3,0,1),s4,0,1)+X{1}, ...
                          plap(X{3},V(1:3))+I(X{4},s4,0,1)+X{2}, ...
                          plap(X{4},V(1:4))+X{3}+X{2}+s3*X{1}};
        S.bc  = @(X,W,U) [pdir(X{2},V(1:2)), pdir(X{3},V(1:3)), pdir(X{4},V(1:3)), {sb(X{4},s4,0)-X{3}, sb(X{4},s4,1)}];

    case 'M5_face'
        % An ODE, 1D, 4D and 5D state; the 4D state is the value of the 5D
        % state on the face s5=0
        X0 = pde_var(1);  x1 = pde_var(s5,[0,1],2);
        z4 = pde_var(V(1:4),dom01(4),2*ones(4,1));  x = pde_var(V,dom01(5),2*ones(5,1));
        S.PDE = [diff(X0,t)==-X0+int(x,V,dom01(5));
                 diff(x1,t)==diff(x1,s5,2)+int(x,V(1:4),dom01(4));
                 diff(z4,t)==lap(z4,V(1:4))+int(x,s5,[0,1])+X0;
                 diff(x,t)==lap(x,V)+z4+x1;
                 dirichlet(x1,s5); dirichlet(z4,V(1:4)); dirichlet(x,V(1:4));
                 subs(x,s5,0)==z4; subs(x,s5,1)==0];
        S.xid = [id(X0),id(x1),id(z4),id(x)];  S.svars = {{},{'s5'},nm(1:4),nm};
        int4 = @(X) I(I(I(I(X,s1,0,1),s2,0,1),s3,0,1),s4,0,1);
        S.rhs = @(X,W,U) {-X{1}+I(int4(X{4}),s5,0,1), D(X{2},s5,2)+int4(X{4}), ...
                          plap(X{3},V(1:4))+I(X{4},s5,0,1)+X{1}, plap(X{4},V)+X{3}+X{2}};
        S.bc  = @(X,W,U) [pdir(X{2},s5), pdir(X{3},V(1:4)), pdir(X{4},V(1:4)), {sb(X{4},s5,0)-X{3}, sb(X{4},s5,1)}];

    case 'M4_product'
        % States on (s1,s2), (s3,s4), (s2,s4) and (s1,...,s4), coupled
        xa = pde_var([s1;s2],dom01(2),[2;2]);  xb = pde_var([s3;s4],dom01(2),[2;2]);
        xd = pde_var([s2;s4],dom01(2),[2;2]);  xc = pde_var(V(1:4),dom01(4),2*ones(4,1));
        S.PDE = [diff(xa,t)==lap(xa,[s1;s2])+int(xc,[s3;s4],dom01(2));
                 diff(xb,t)==lap(xb,[s3;s4])+int(xc,[s1;s2],dom01(2));
                 diff(xd,t)==lap(xd,[s2;s4])+int(xc,[s1;s3],dom01(2));
                 diff(xc,t)==lap(xc,V(1:4))+xa+xb+s1*s3*xd;
                 dirichlet(xa,[s1;s2]); dirichlet(xb,[s3;s4]); dirichlet(xd,[s2;s4]); dirichlet(xc,V(1:4))];
        S.xid = [id(xa),id(xb),id(xd),id(xc)];  S.svars = {{'s1','s2'},{'s3','s4'},{'s2','s4'},nm(1:4)};
        S.rhs = @(X,W,U) {plap(X{1},[s1;s2])+I(I(X{4},s3,0,1),s4,0,1), ...
                          plap(X{2},[s3;s4])+I(I(X{4},s1,0,1),s2,0,1), ...
                          plap(X{3},[s2;s4])+I(I(X{4},s1,0,1),s3,0,1), ...
                          plap(X{4},V(1:4))+X{1}+X{2}+s1*s3*X{3}};
        S.bc  = @(X,W,U) [pdir(X{1},[s1;s2]), pdir(X{2},[s3;s4]), pdir(X{3},[s2;s4]), pdir(X{4},V(1:4))];

    case {'W4_inputs','W4_input_coef'}
        % Finite, distributed and face disturbances; a finite control in the
        % PDE and a face control in the Neumann BC; outputs with
        % feedthrough, and a face sensor. W4_input_coef: the control in the
        % BC multiplied by a function of s1,s2,s3, which is refused
        g = s1*(1-s1)*s2*(1-s2)*s3*(1-s3);
        x = pde_var(V(1:4),dom01(4),2*ones(4,1));
        w1 = pde_var('in',1);  w2 = pde_var('in',1,[s1;s2],dom01(2));  w3 = pde_var('in',1,V(1:3),dom01(3));
        u = pde_var('control',1);  u2 = pde_var('control',1,V(1:3),dom01(3));
        z = pde_var('out',1);  y = pde_var('sense',1,V(1:3),dom01(3));
        if strcmp(name,'W4_input_coef')
            ubc = g*u;  S.expect = 'single element of a single input';
        else
            ubc = u2;
        end
        S.PDE = [diff(x,t)==lap(x,V(1:4))+s1*w1+w2+s2*u;
                 z==int(x,V(1:4),dom01(4))+w1;  y==subs(x,s4,1);
                 dirichlet(x,V(1:3)); subs(x,s4,0)==w3; subs(diff(x,s4),s4,1)==ubc];
        S.xid = id(x);  S.wid = [id(w1),id(w2),id(w3)];  S.uid = [id(u),id(u2)];
        S.zid = id(z);  S.yid = id(y);  S.svars = {nm(1:4)};
        S.wgen = @(c) {c(1), c(2)+c(3)*s1+c(4)*s2^2+c(5)*s1*s2, g*(c(6)+c(7)*s1+c(8)*s3^2)};
        S.ugen = @(c) {c(9), g*(c(10)+c(11)*s2+c(12)*s1*s3)};
        S.rhs = @(X,W,U) {plap(X{1},V(1:4))+s1*W{1}+W{2}+s2*U{1}};
        S.z   = @(X,W,U) {I(I(I(I(X{1},s1,0,1),s2,0,1),s3,0,1),s4,0,1)+W{1}};
        S.y   = @(X,W,U) {sb(X{1},s4,1)};
        S.bc  = @(X,W,U) [pdir(X{1},V(1:3)), {sb(X{1},s4,0)-W{3}, sb(D(X{1},s4,1),s4,1)-U{2}}];

    case {'P4_periodic','P4_periodic2'}
        % Periodic in s1 (and in s3 for P4_periodic2), Dirichlet otherwise:
        % rank-deficient directions, so auxiliary states and K. The states
        % are built from q(s), periodic in value and slope
        x = pde_var(V(1:4),dom01(4),2*ones(4,1));
        per = @(v) [subs(x,v,0)-subs(x,v,1)==0; subs(diff(x,v),v,0)-subs(diff(x,v),v,1)==0];
        pres = @(X,v) {sb(X,v,0)-sb(X,v,1), sb(D(X,v,1),v,0)-sb(D(X,v,1),v,1)};
        q = @(s,c) c(1)+c(2)*s^2*(1-s)^2+c(3)*s^3*(1-s)^3;
        if strcmp(name,'P4_periodic')
            S.PDE = [diff(x,t)==lap(x,V(1:4))+x; per(s1); dirichlet(x,V(2:4))];
            S.bc  = @(X,W,U) [pres(X{1},s1), pdir(X{1},V(2:4))];
            S.xs  = @(c) {q(s1,c(1:3))*s2*(1-s2)*s3*(1-s3)*s4*(1-s4)*(c(4)+c(5)*s2+c(6)*s3*s4) + ...
                          q(s1,c(7:9))*s2*(1-s2)*s3*(1-s3)*s4*(1-s4)*s2^2};
        else
            S.PDE = [diff(x,t)==lap(x,V(1:4))+x; per(s1); per(s3); dirichlet(x,[s2;s4])];
            S.bc  = @(X,W,U) [pres(X{1},s1), pres(X{1},s3), pdir(X{1},[s2;s4])];
            S.xs  = @(c) {q(s1,c(1:3))*q(s3,c(4:6))*s2*(1-s2)*s4*(1-s4)*(c(7)+c(8)*s2+c(9)*s4^2)};
        end
        S.xid = id(x);  S.svars = {nm(1:4)};
        S.rhs = @(X,W,U) {plap(X{1},V(1:4))+X{1}};

    case 'E4_corner'
        % A redundant condition at the corner x(0,0,s3,s4)=0, implied by the
        % Dirichlet faces: accepted
        x = pde_var(V(1:4),dom01(4),2*ones(4,1));
        S.PDE = [diff(x,t)==lap(x,V(1:4)); dirichlet(x,V(1:4)); subs(subs(x,s1,0),s2,0)==0];
        S.xid = id(x);  S.svars = {nm(1:4)};
        S.rhs = @(X,W,U) {plap(X{1},V(1:4))};
        S.bc  = @(X,W,U) [pdir(X{1},V(1:4)), {sb(sb(X{1},s1,0),s2,0)}];

    case {'BVP4','BVP5'}
        % No dynamics: Dirichlet, Neumann and Robin BCs and a face input
        % (BVP4), Dirichlet (BVP5)
        N = S.dim;  vv = V(1:N);
        x = pde_var(vv,dom01(N),2*ones(N,1));
        if N==4
            w = pde_var('in',1,V(1:3),dom01(3));
            S.PDE = [dirichlet(x,[s1;s2]); subs(diff(x,s3),s3,0)==0; subs(x,s3,1)==0;
                     subs(x,s4,0)==w; subs(diff(x,s4),s4,1)+3*subs(x,s4,1)==0];
            S.wid = id(w);
            % w vanishes where x is fixed by Dirichlet BCs, w_s3(s3=0)=0
            S.wgen = @(c) {s1*(1-s1)*s2*(1-s2)*(1-s3)^2*(1+2*s3)*(c(1)+c(2)*s1+c(3)*s2)};
            S.bc = @(X,W,U) [pdir(X{1},[s1;s2]), {sb(D(X{1},s3,1),s3,0), sb(X{1},s3,1), ...
                             sb(X{1},s4,0)-W{1}, sb(D(X{1},s4,1),s4,1)+3*sb(X{1},s4,1)}];
        else
            S.PDE = dirichlet(x,vv);
            S.bc = @(X,W,U) pdir(X{1},vv);
        end
        S.xid = id(x);  S.svars = {nm(1:N)};

    case 'D3_delay'
        % 3D heat with a delayed state: a 4D history state is added (only
        % the conversion and the dimensions of the blocks are checked)
        x = pde_var(V(1:3),dom01(3),2*ones(3,1));
        S.PDE = [diff(x,t)==lap(x,V(1:3))+0.5*subs(x,t,t-0.5); dirichlet(x,V(1:3))];
        S.check = false;

    case 'H5_coupled'
        % 5D heat coupled to a 1D state in s5 and an ODE
        X0 = pde_var(1);  x1 = pde_var(s5,[0,1],2);  x = pde_var(V,dom01(5),2*ones(5,1));
        S.PDE = [diff(X0,t)==-X0+int(x,V,dom01(5));
                 diff(x1,t)==diff(x1,s5,2)+int(x,V(1:4),dom01(4));
                 diff(x,t)==lap(x,V)+x1+X0;
                 dirichlet(x1,s5); dirichlet(x,V)];
        S.xid = [id(X0),id(x1),id(x)];  S.svars = {{},{'s5'},nm};
        int4 = @(X) I(I(I(I(X,s1,0,1),s2,0,1),s3,0,1),s4,0,1);
        S.rhs = @(X,W,U) {-X{1}+I(int4(X{3}),s5,0,1), D(X{2},s5,2)+int4(X{3}), plap(X{3},V)+X{2}+X{1}};
        S.bc  = @(X,W,U) [pdir(X{2},s5), pdir(X{3},V)];

    case 'E4_tangential'
        x = pde_var(V(1:4),dom01(4),2*ones(4,1));
        S.PDE = [diff(x,t)==lap(x,V(1:4)); subs(diff(x,s2),s1,0)==0; subs(x,s1,1)==0; dirichlet(x,V(2:4))];
        S.expect = 'derivative of a state along a direction other than';
    case 'E4_fewBCs'
        x = pde_var(V(1:4),dom01(4),2*ones(4,1));
        S.PDE = [diff(x,t)==lap(x,V(1:4)); dirichlet(x,V([1 2 4])); subs(x,s3,0)==0];
        S.expect = 'should be at least the total order';
    case 'E5_trivial'
        x = pde_var(V,dom01(5),2*ones(5,1));
        S.PDE = [diff(x,t)==lap(x,V); dirichlet(x,V(1:4)); subs(x,s5,0)==0; subs(x,s5,1)-subs(x,s5,1)==0];
        S.expect = 'reduces to 0=0';
end
end



%% ------------------------------------------------------------------------
function E = dirichlet(x,vv)
% E = DIRICHLET(X,VV) stacks the conditions x=0 at s=0 and s=1 for each s
% in VV.
%
% INPUT
% - x:   'pde_struct', a state declared with pde_var;
% - vv:  n x 1 'polynomial', the variables;
%
% OUTPUT
% - E:   'pde_struct', the 2n conditions.
E = [subs(x,vv(1),0)==0; subs(x,vv(1),1)==0];
for i=2:numel(vv)
    E = [E; subs(x,vv(i),0)==0; subs(x,vv(i),1)==0];                       %#ok<AGROW>
end
end

function L = lap(x,vv)
% L = LAP(X,VV) returns the sum of the second derivatives of the state X
% (a 'pde_struct') along the variables VV (n x 1 'polynomial').
L = diff(x,vv(1),2);
for i=2:numel(vv), L = L+diff(x,vv(i),2); end
end

function L = plap(q,vv)
% L = PLAP(Q,VV) returns the same for the 'polynomial' Q.
L = 0*polynomial(q);
for i=1:numel(vv), L = L+dk(q,vv(i),2); end
end

function R = pdir(q,vv)
% R = PDIR(Q,VV) returns the 1 x 2n cell of the values of the 'polynomial'
% Q at s=0 and s=1 for each s in VV (n x 1 'polynomial').
R = {};
for i=1:numel(vv)
    R = [R, {subs(polynomial(q),vv(i),0), subs(polynomial(q),vv(i),1)}];   %#ok<AGROW>
end
end

function q = dk(q,v,k)
% Q = DK(Q,V,K) differentiates the 'polynomial' Q K times in the variable V.
q = polynomial(q);
for i=1:k, q = diff(q,v); end
end

function i = id(p)
% I = ID(P) returns the ID of the single state, input or output declared by
% pde_var as P (a 'pde_struct').
f = {'x','w','u','z','y'};
for k=1:numel(f)
    if ~isempty(p.(f{k})), i = p.(f{k}){1}.ID; return, end
end
end
