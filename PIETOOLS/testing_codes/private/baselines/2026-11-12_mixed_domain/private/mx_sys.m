function S = mx_sys(name,p)
% S = MX_SYS(NAME,P) returns one mixed-domain test system: its PDE, and its
% semantics written by hand, independently of the parser. The scratch
% version of MIXED_DOMAIN_SYSTEMS that N3-N10 ran, with S32 (library PDE 43
% by the delay syntax) in addition.
% INPUT
% - name:   char, the system name;
% - p:      (optional) 'struct' of parameters overriding the defaults.
% OUTPUT
% - S:      'struct': PDE ('pde_struct', not initialized); svars{k},
%           wvars{j} (variables of state k, input j); xs @(c,W) (states
%           satisfying the BCs); rhs @(X,W) (the right-hand side of each
%           state equation); bc @(X,W) (residuals that vanish iff the BCs
%           hold); z, y @(X,W) (outputs); expect ('' for a supported
%           system, else the expected error).

% Initial coding MMP, 11/12/2026.
if nargin<2, p = struct(); end
pvar s1 s2 s3 t
D = @(q,v,k) dk(q,v,k);
nd  = @(s,c) c(1)+c(2)*(2*s^3-3*s^2)+s^2*(1-s)^2*(c(3)+c(4)*s);        % q'(0)=q'(1)=0
nD  = @(s,c) c(1)*(1-s)^2*(1+2*s)+s^2*(1-s)*(c(2)+c(3)*s);            % q'(0)=0, q(1)=0
dd  = @(s,c) s*(1-s)*(c(1)+c(2)*s+c(3)*s^2);                          % q(0)=q(1)=0
S = struct('name',name,'PDE',[],'svars',{{}},'wvars',{{}},'xs',[],'rhs',[],'bc',[], ...
           'z',[],'y',[],'expect','','par',p);
none = @(X,W) {};
S.z = none;  S.y = none;
clear stateNameGenerator
switch name
    case 'S1_ode_heat1d_bc'
        % ODE -> 1-D heat through a Dirichlet BC (cascade)
        a = getp(p,'a',1);  r = getp(p,'r',0);
        X = pde_var(1);  x1 = pde_var(1,s1,[0,1]);
        S.PDE = [diff(X,t)==-a*X; diff(x1,t)==diff(x1,s1,2)+r*x1;
                 subs(x1,s1,0)==X; subs(x1,s1,1)==0];
        S.svars = {{},{'s1'}};
        S.xs  = @(c,W) {c(1), c(1)*(1-s1)+dd(s1,c(2:4))};
        S.rhs = @(Xp,W) {-a*Xp{1}, D(Xp{2},s1,2)+r*Xp{2}};
        S.bc  = @(Xp,W) {subs(Xp{2},s1,0)-Xp{1}, subs(Xp{2},s1,1)};
    case 'S2_ode_heat2d_edge'
        % ODE -> 2-D heat through the edge s2=0; Neumann in s1 (aux state, K)
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
        % 1-D heat on the edge s2=0 of a 2-D heat (cascade); reaction r in x1
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
        % 2-D heat with a dynamic (Wentzell) BC on s2=0: x1 = x2(s1,0),
        % x1_t = x1_s1s1 + c*x2_s2(s1,0) + r*x1, x2_t = lap x2 + r*x2
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
        % ODE -> 1-D (Neumann BC) -> 2-D (edge s2=0, and the ODE in the s1=0 BC)
        a = getp(p,'a',1);
        X = pde_var(1);  x1 = pde_var(1,s1,[0,1]);  x2 = pde_var(1,[s1;s2],[0,1;0,1]);
        S.PDE = [diff(X,t)==-a*X; diff(x1,t)==diff(x1,s1,2);
                 diff(x2,t)==diff(x2,s1,2)+diff(x2,s2,2);
                 subs(diff(x1,s1),s1,0)==X; subs(x1,s1,1)==0;
                 subs(diff(x2,s1),s1,0)==(1-s2)*X; subs(x2,s1,1)==0;
                 subs(x2,s2,0)==x1; subs(x2,s2,1)==0];
        S.svars = {{},{'s1'},{'s1','s2'}};
        S.xs  = @(c,W) mx_s5(c,s1,s2,nD);
        S.rhs = @(Xp,W) {-a*Xp{1}, D(Xp{2},s1,2), D(Xp{3},s1,2)+D(Xp{3},s2,2)};
        S.bc  = @(Xp,W) {subs(D(Xp{2},s1,1),s1,0)-Xp{1}, subs(Xp{2},s1,1), ...
                       subs(D(Xp{3},s1,1),s1,0)-(1-s2)*Xp{1}, subs(Xp{3},s1,1), ...
                       subs(Xp{3},s2,0)-Xp{2}, subs(Xp{3},s2,1)};
    case 'S6_indomain_1d_2d'
        % 1-D <-> 2-D in-domain coupling (integral over s2, multiplier),
        % Neumann in s2: exact spectrum -m^2 pi^2 +- g (and -(m^2+n^2)pi^2)
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
        % ODE <-> 2-D in-domain coupling, Neumann on all edges (3 aux blocks)
        % mean mode: [-a g; g -b]
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
        % 1-D states on two edges (s2=0 lower, s1=1 upper) of a 2-D state
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
        % ODE <-> 2-D: the ODE fixes the edge value, the edge flux drives the
        % ODE; w into the ODE, z = X. Mean mode: lam + c*mu*coth(mu) - a = 0
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
        % ODE in the Dirichlet edge of a Dirichlet square: corner data clash
        X = pde_var(1);  x2 = pde_var(1,[s1;s2],[0,1;0,1]);
        S.PDE = [diff(X,t)==-X; diff(x2,t)==diff(x2,s1,2)+diff(x2,s2,2);
                 subs(x2,s1,0)==0; subs(x2,s1,1)==0;
                 subs(x2,s2,0)==X; subs(x2,s2,1)==0];
        S.svars = {{},{'s1','s2'}};
        S.expect = 'compatibility';
    case 'S11_incompatible_neumann_corner'
        % x2_s1(0,s2)=0 and x2(s1,0)=x1(s1) need x1'(0)=0; x1 is Dirichlet,
        % so x1'(0) is free: no H^{2,2} state unless x1'(0)=0
        x1 = pde_var(1,s1,[0,1]);  x2 = pde_var(1,[s1;s2],[0,1;0,1]);
        S.PDE = [diff(x1,t)==diff(x1,s1,2); diff(x2,t)==diff(x2,s1,2)+diff(x2,s2,2);
                 subs(x1,s1,0)==0; subs(x1,s1,1)==0;
                 subs(diff(x2,s1),s1,0)==0; subs(x2,s1,1)==0;
                 subs(x2,s2,0)==x1; subs(x2,s2,1)==0];
        S.svars = {{'s1'},{'s1','s2'}};
        S.rhs = @(Xp,W) {D(Xp{1},s1,2), D(Xp{2},s1,2)+D(Xp{2},s2,2)};
        S.bc  = @(Xp,W) {subs(Xp{1},s1,0), subs(Xp{1},s1,1), subs(D(Xp{2},s1,1),s1,0), subs(Xp{2},s1,1), subs(Xp{2},s2,0)-Xp{1}, subs(Xp{2},s2,1)};
        S.expect = 'compatibility';
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
        S.expect = 'compatibility';
    case 'S13_chain_1d_2d_3d'
        % 1-D edge -> 2-D face -> 3-D heat; the 1-D state also enters the
        % 3-D face s2=0 as (1-s3)*x1(s1), so the edge s2=s3=0 is compatible
        x1 = pde_var(1,s1,[0,1]);  x2 = pde_var(1,[s1;s2],[0,1;0,1]);  x3 = pde_var(1,[s1;s2;s3],[0,1;0,1;0,1]);
        S.PDE = [diff(x1,t)==diff(x1,s1,2); diff(x2,t)==diff(x2,s1,2)+diff(x2,s2,2);
                 diff(x3,t)==diff(x3,s1,2)+diff(x3,s2,2)+diff(x3,s3,2);
                 subs(x1,s1,0)==0; subs(x1,s1,1)==0;
                 subs(x2,s1,0)==0; subs(x2,s1,1)==0; subs(x2,s2,0)==x1; subs(x2,s2,1)==0;
                 subs(x3,s1,0)==0; subs(x3,s1,1)==0; subs(x3,s2,0)==(1-s3)*x1; subs(x3,s2,1)==0;
                 subs(x3,s3,0)==x2; subs(x3,s3,1)==0];
        S.svars = {{'s1'},{'s1','s2'},{'s1','s2','s3'}};
        S.xs  = @(c,W) mx_s13(c,s1,s2,s3,dd);
        S.rhs = @(Xp,W) {D(Xp{1},s1,2), D(Xp{2},s1,2)+D(Xp{2},s2,2), D(Xp{3},s1,2)+D(Xp{3},s2,2)+D(Xp{3},s3,2)};
        S.bc  = @(Xp,W) {subs(Xp{1},s1,0), subs(Xp{1},s1,1), subs(Xp{2},s1,0), subs(Xp{2},s1,1), ...
                       subs(Xp{2},s2,0)-Xp{1}, subs(Xp{2},s2,1), subs(Xp{3},s1,0), subs(Xp{3},s1,1), ...
                       subs(Xp{3},s2,0)-(1-s3)*Xp{1}, subs(Xp{3},s2,1), subs(Xp{3},s3,0)-Xp{2}, subs(Xp{3},s3,1)};
    case 'S19_2d_integral_in_1d_bc'
        % a HIGHER-dimensional state in the BC of a 1-D state:
        % x1(0) = int int x2, and x2 driven in-domain by x1
        x1 = pde_var(1,s1,[0,1]);  x2 = pde_var(1,[s1;s2],[0,1;0,1]);
        S.PDE = [diff(x1,t)==diff(x1,s1,2); diff(x2,t)==diff(x2,s1,2)+diff(x2,s2,2)+x1;
                 subs(x1,s1,0)==int(x2,[s1;s2],[0,1;0,1]); subs(x1,s1,1)==0;
                 subs(x2,s1,0)==0; subs(x2,s1,1)==0; subs(x2,s2,0)==0; subs(x2,s2,1)==0];
        S.svars = {{'s1'},{'s1','s2'}};
        S.xs  = @(c,W) mx_s19(c,s1,s2,dd);
        S.rhs = @(Xp,W) {D(Xp{1},s1,2), D(Xp{2},s1,2)+D(Xp{2},s2,2)+Xp{1}};
        S.bc  = @(Xp,W) {subs(Xp{1},s1,0)-int(int(Xp{2},s1,0,1),s2,0,1), subs(Xp{1},s1,1), ...
                       subs(Xp{2},s1,0), subs(Xp{2},s1,1), subs(Xp{2},s2,0), subs(Xp{2},s2,1)};
    case 'S20_reversed_var_order'
        % S4 with the 2-D state declared on [s2;s1] (reversed order)
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
        % S4 on [a,b] x [c,d] = [-1,2] x [0.5,3], with an s1-dependent
        % coefficient on the edge source: x2(s1,c) = (2+s1)*x1(s1)
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
        S.dom = [a b; c1 d1];
    case 'S22_vector_states'
        % a 2-vector 1-D state on the edge of a scalar 2-D state, coupled by
        % matrices: x1_t = x1_ss + [0 1;-1 0] x1, x2(s1,0) = [1 2] x1
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
    case 'S14_input_in_source_group_2d'
        % boundary input w in the s1=1 BC of a 2-D state whose s2=0 edge is
        % a 1-D state: documented as unsupported (solve_PIETOOLS_BVP NOTES)
        x1 = pde_var(1,s1,[0,1]);  x2 = pde_var(1,[s1;s2],[0,1;0,1]);  w = pde_var('in',1);
        S.PDE = [diff(x1,t)==diff(x1,s1,2); diff(x2,t)==diff(x2,s1,2)+diff(x2,s2,2);
                 subs(x1,s1,0)==0; subs(diff(x1,s1),s1,1)==w;
                 subs(x2,s1,0)==0; subs(diff(x2,s1),s1,1)==(1-s2)*w;
                 subs(x2,s2,0)==x1; subs(x2,s2,1)==0];
        S.svars = {{'s1'},{'s1','s2'}};
        S.expect = 'input';
    case 'S15_ode_input_same_bc'
        % 1-D: x1(0) = X + w (input lifted into the ODE source), z = int x1
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
        % S6 with a distributed disturbance in the 2-D equation, a scalar
        % output and an L2[s1] output (trace of x2 on s2=1)
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
        % (order 0 in s1, no BCs in s1): orders of source and target differ
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
        % library PDE 43 (Caliskan 2009): heat with interior delay, the delay
        % as a 2-D transport state. sgn=+1 is the library file as coded
        % (+a1*x(t-tau): unstable for every tau); sgn=-1 is the system of its
        % header and of the reference, stable for tau < 1.0347 (c=a1=1, a0=1.9)
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
    case 'S23_hinf_indomain'
        % S6 (g) with a distributed disturbance w(s1,s2) in the 2-D
        % equation and z = int x1: exact gain by a series over s1 modes
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
        % S2 (ODE a -> 2-D edge, Neumann in s1) with w into the ODE and
        % z = int int x2: only the s1-mean mode is excited
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
        % S4 (c, r) with a disturbance w(s1) on the edge state, z = int x1
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
    case 'S32_delay_syntax'
        % library PDE 43 written with the delay syntax: the converter's
        % expand_delays builds the transport state itself
        tau = getp(p,'tau',0.5);  sg = getp(p,'sgn',-1);
        x1 = pde_var(1,s1,[0,pi]);
        S.PDE = [diff(x1,t)==diff(x1,s1,2)+1.9*x1+sg*subs(x1,t,t-tau);
                 subs(x1,s1,0)==0; subs(x1,s1,pi)==0];
        S.svars = {{'s1'}};
    otherwise
        error('mx_sys:name','Unknown system %s',name);
end
end

function v = getp(p,f,d)
if isfield(p,f), v = p.(f); else, v = d; end
end

function q = dk(q,v,k)
for i=1:k
    q = diff(q,v);
end
end

function X = mx_s5(c,s1,s2,nD)
Xo = c(1);
x1 = Xo*(s1-1)+nD(s1,c(2:4));
g  = nD(s1,c(5:7))*(c(8)+c(9)*s2)+nD(s1,c(10:12))*s2^2;
x2 = (1-s2)*x1+s2*(1-s2)*g;
X = {Xo,x1,x2};
end

function X = mx_s19(c,s1,s2,dd)
% x2 Dirichlet on the square; x1(1)=0 and x1(0) = int int x2
x2 = s1*(1-s1)*s2*(1-s2)*(c(4)+c(5)*s1+c(6)*s2^2);
m  = int(int(x2,s1,0,1),s2,0,1);
x1 = m*(1-s1)+dd(s1,c(1:3));
X = {x1,x2};
end

function X = mx_s13(c,s1,s2,s3,dd)
x1 = dd(s1,c(1:3));
x2 = (1-s2)*x1+s1*(1-s1)*s2*(1-s2)*(c(4)+c(5)*s1*s2);
x3 = (1-s3)*x2+s1*(1-s1)*s2*(1-s2)*s3*(1-s3)*(c(6)+c(7)*s1*s2*s3);
X = {x1,x2,x3};
end
