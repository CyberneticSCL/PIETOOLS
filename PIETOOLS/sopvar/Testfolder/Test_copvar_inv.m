function Test_copvar_inv(quick)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% TEST_COPVAR_INV([QUICK]) asserts on the container inverse (@sopvar/inv,
% @copvar/inv) and the gain reconstruction (getController_sop,
% getObserver_sop), 1-D:
%   (a) hand-built 3-PI blocks (the cases of the 10/09/2026 probe) and a
%       4-PI operator with an R^2 space: the operator residuals
%       max |P*Pinv - I| and max |Pinv*P - I| (parameters on a grid, kernels
%       on their triangles) are below 1e-6, computed with the CONTAINER
%       composition and, independently, with the opvar composition of the
%       converted inverse; the two agree;
%   (b) the container inverse agrees with the stock inv_opvar_2 on the same
%       operators to 1e-5 on the grid (two implementations of one formula);
%   (c) the identity and a pure matrix invert exactly;
%   (d) QUICK false (default true): the synthesis executives on the test
%       plant of test_executives_sop, Hinf_control_sop / Hinf_estimator_sop
%       (MOSEK, light): the container gains agree with the stock gains
%       (getController / getObserver on the same solved P, Z) to 1e-4
%       relative, and the closed-loop numerical gain of the container K is
%       below the certified gamma.
% Each check prints its line; the function errors at the first failure.
%
% Initial coding MMP, 10/09/2026
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if nargin<1,    quick = true;   end
warning('off','sopvar:noncanonicalMultiplier');
pvar s s_dum
n = 0;
cases = {};
opvar Pc; Pc.I = [0,1]; Pc.var1 = s; Pc.var2 = s_dum;
Pc.R.R0 = 1+0.5*s;  Pc.R.R1 = 0.3*s*s_dum;  Pc.R.R2 = 0.3*(s+s_dum);
cases{end+1} = {'c1 R0=1+0.5s, R1=0.3 s t, R2=0.3(s+t)',Pc};
opvar Pc; Pc.I = [0,1]; Pc.var1 = s; Pc.var2 = s_dum;
Pc.R.R0 = polynomial(1);  Pc.R.R1 = 0.3*s;  Pc.R.R2 = 0.3*s_dum;
cases{end+1} = {'c2 R0=1, R1=0.3 s, R2=0.3 t',Pc};
opvar Pc; Pc.I = [0,1]; Pc.var1 = s; Pc.var2 = s_dum;
Pc.R.R0 = 1+0.5*s;  Pc.R.R1 = 0.3*s;  Pc.R.R2 = 0.3*s_dum;
cases{end+1} = {'c3 R0=1+0.5s, R1=0.3 s, R2=0.3 t',Pc};
opvar Pc; Pc.I = [0,1]; Pc.var1 = s; Pc.var2 = s_dum;
Pc.R.R0 = [2+s, 0.2*s; 0.2*s, 1.5];  Pc.R.R1 = [0.2*s*s_dum, 0.1; 0.1*s_dum, 0.3*s];  Pc.R.R2 = [0.1, 0.2*s; 0.3*s_dum, 0.1*s*s_dum];
cases{end+1} = {'c4 2x2 non-separable',Pc};
opvar Pc; Pc.I = [0,1]; Pc.var1 = s; Pc.var2 = s_dum;
Pc.R.R0 = 1+0.5*s;  Pc.R.R1 = 0.3*s*s_dum;  Pc.R.R2 = 0.3*s*s_dum;
cases{end+1} = {'c5 separable R1=R2',Pc};
opvar Pc; Pc.I = [0,1]; Pc.var1 = s; Pc.var2 = s_dum;
Pc.R.R0 = 1+0.5*s;  Pc.R.R1 = 0.3*s*s_dum;
cases{end+1} = {'c6 lower kernel only',Pc};
opvar Pc; Pc.I = [0,1]; Pc.var1 = s; Pc.var2 = s_dum;
Pc.P = [2, 0.1; 0, 1.5];   Pc.Q1 = [0.3*s, 0.1; 0.2, 0.1*s^2];   Pc.Q2 = [0.2*s^2, 0.1*s; 0.1, 0.3];
Pc.R.R0 = [2+s, 0.2*s; 0.2*s, 1.5];  Pc.R.R1 = [0.2*s*s_dum, 0.1; 0.1*s_dum, 0.3*s];  Pc.R.R2 = [0.1, 0.2*s; 0.3*s_dum, 0.1*s*s_dum];
cases{end+1} = {'c7 4-PI, R^2 x L2^2',Pc};
sg = linspace(0,1,21);
for c = 1:numel(cases)
    nm = cases{c}{1};   Pc = cases{c}{2};
    Iop = mat2opvar(eye(sum(Pc.dim(:,1))),Pc.dim(:,1),[s,s_dum],[0,1]);
    Pm = opvar2copvar(Pc);
    t0 = tic;   [Pinv,info] = inv(Pm);     tk = toc(t0);
    Pinv_op = copvar2opvar(Pinv);
    % (a) residuals: container composition, and opvar composition
    E1 = copvar2opvar(Pm*Pinv) - Iop;       E2 = copvar2opvar(Pinv*Pm) - Iop;
    E3 = Pc*Pinv_op - Iop;
    r1 = opres(E1,sg);  r2 = opres(E2,sg);  r3 = opres(E3,sg);
    n = check(n,r1<1e-6 && r2<1e-6 && abs(r1-r3)<1e-8, ...
        sprintf('%s: %.2f s, d %d, ranks [%d %d], fit relrms %.1e; |P*Pinv-I| %.2e (opvar %.2e), |Pinv*P-I| %.2e', ...
        nm,tk,info.d,info.r1,info.r2,max(info.relrms),r1,r3,r2));
    % (b) against inv_opvar_2 on the same operator
    if all(Pc.dim(1,:)==0)
        P2 = inv_opvar_2(Pc);
        dR = max([pdiff(Pinv_op.R.R0,P2.R.R0,s,s_dum,sg), pdiff(Pinv_op.R.R1,P2.R.R1,s,s_dum,sg), pdiff(Pinv_op.R.R2,P2.R.R2,s,s_dum,sg)]);
        n = check(n,dR<1e-5,sprintf('   vs inv_opvar_2: max parameter difference on the grid %.2e',dR));
    end
end
% (e) the direct construction (no inverse object) against the inverse route,
%     on the 4-PI operator c7 with a free operator Z: R^2 x L2^2 -> R^2
opvar Zc; Zc.I = [0,1]; Zc.var1 = s; Zc.var2 = s_dum;
Zc.P = [1, 0.5; 0, 2];   Zc.Q1 = [s, 1; 0.5*s^2, 2*s];
Pc = cases{7}{2};
t0 = tic;   [Kd,infod] = getController_direct_sop(Pc,Zc);    tk = toc(t0);
Ki = copvar2opvar(getController_sop(Pc,Zc));
dP = max(abs(Kd.P - double(Ki.P)),[],'all');
q = zeros(size(Kd.Q1));
for i = 1:numel(Kd.s),  q(:,:,i) = double(subs(Ki.Q1,s,Kd.s(i)));    end
dQ = max(abs(q(:) - Kd.Q1(:)))/max(abs(Kd.Q1(:)));
zmax = opres(Zc,sg);
rKd = opres(Kd.op*Pc - Zc,sg)/zmax;
x1 = [0.3; -0.2];   x2 = @(t) [sin(pi*t); t.^2];
u1 = Kd.apply(x1,x2);
u2 = double(Ki.P)*x1;
for i = 1:numel(Kd.s),  u2 = u2 + Kd.w(i)*(double(subs(Ki.Q1,s,Kd.s(i)))*x2(Kd.s(i)));   end
n = check(n,dP<1e-8 && dQ<1e-7 && rKd<1e-6 && max(abs(u1-u2))<1e-8 && isa(Kd.c,'copvar'), ...
    sprintf('direct gains on c7: %.2f s; K1 vs Z*inv(P) %.1e, K2 on the grid %.1e (rel), |K P - Z|/|Z| of the fitted opvar %.1e (fit d %d relrms %.1e), apply vs quadrature of the inverse route %.1e, condT %.2g', ...
    tk,dP,dQ,rKd,Kd.fit.d,Kd.fit.relrms,max(abs(u1-u2)),infod.condT));
% (c) identity and a pure matrix
Pm = opvar2copvar(mat2opvar(eye(3),[1;2],[s,s_dum],[0,1]));
Pinv = inv(Pm);     E = copvar2opvar(Pm*Pinv) - mat2opvar(eye(3),[1;2],[s,s_dum],[0,1]);
n = check(n,opres(E,sg)<1e-12,sprintf('identity on R x L2^2: residual %.2e',opres(E,sg)));
opvar Pc; Pc.I = [0,1]; Pc.var1 = s; Pc.var2 = s_dum;  Pc.P = [2 1; 0 3];
Pm = opvar2copvar(Pc);   Pinv = inv(Pm);
n = check(n,norm(full(Pinv.C{1,1}.params{1})-inv([2 1;0 3]))<1e-14,'pure matrix on R^2');
% (d) the executives
if ~quick
    st = lpisettings('light');
    if ~isempty(which('mosekopt')),     st.sos_opts.solver = 'mosek';   end
    st.sop.verbose = false;
    pvar t
    x = pde_var(1,s,[0,1]);     w = pde_var('in',1);     u = pde_var('control',1);
    y = pde_var('sense',1);     z = pde_var('out',2);
    PDE = [diff(x,t,1)==diff(x,s,2)+0.5*pi^2*x+(s*(1-s))*w+u;
           z==[int(x,s,[0,1]); u];     y==int(x,s,[0,1]);
           subs(x,s,0)==0;  subs(x,s,1)==0];
    evalc('PIE = initialize(convert(PDE));');
    % the gains are judged by their defining equations, K P = Z and P L = Z
    % (residual on the grid relative to max |Z|), the stock gains alongside
    [~,K,gam,P,Z,info] = PIETOOLS_Hinf_control_sop(PIE,st);
    Kst = getController(P,Z);
    zmax = opres(Z,sg);
    rK = opres(K*P - Z,sg)/zmax;     rKst = opres(Kst*P - Z,sg)/zmax;
    dK = pdiff(K.Q1,Kst.Q1,s,s_dum,sg)/max(1e-300,pmax(Kst.Q1,s,s_dum,sg));
    W = pie_witness_sop(closedLoopPIE(PIE,K),'gain',struct('N_cheb',24));
    n = check(n,isa(info.G,'copvar') && isa(K,'opvar') && rK<1e-6 && W.gain<=gam*(1+1e-6), ...
        sprintf('Hinf_control_sop: gamma %.7g, closed-loop numerical gain %.7g; |K P - Z|/|Z| %.1e (stock getController %.1e), K vs stock rel %.1e, inverse d %d relrms %.1e', ...
        gam,W.gain,rK,rKst,dK,info.inv.d,max(info.inv.relrms)));
    % the direct construction on the solved P, Z against the inverse route
    Kd = getController_direct_sop(P,Z);
    q = zeros(size(Kd.Q1));
    for i = 1:numel(Kd.s),  q(:,:,i) = double(subs(K.Q1,s,Kd.s(i)));    end
    dQ = max(abs(q(:) - Kd.Q1(:)))/max(abs(Kd.Q1(:)));
    rKd = opres(Kd.op*P - Z,sg)/zmax;
    Wd = pie_witness_sop(closedLoopPIE(PIE,Kd.op),'gain',struct('N_cheb',24));
    n = check(n,dQ<1e-7 && rKd<1e-6 && Wd.gain<=gam*(1+1e-6), ...
        sprintf('   direct gains: K2 vs inverse route %.1e (rel), |K P - Z|/|Z| %.1e (fit d %d relrms %.1e), closed-loop numerical gain %.7g', ...
        dQ,rKd,Kd.fit.d,Kd.fit.relrms,Wd.gain));
    [~,L,gam,P,Z,info] = PIETOOLS_Hinf_estimator_sop(PIE,st);
    Lst = getObserver(P,Z);
    zmax = opres(Z,sg);
    rL = opres(P*L - Z,sg)/zmax;     rLst = opres(P*Lst - Z,sg)/zmax;
    dL = max([pdiff(L.Q2,Lst.Q2,s,s_dum,sg), max(abs(double(L.P)-double(Lst.P)))])/max(1e-300,max([pmax(Lst.Q2,s,s_dum,sg), max(abs(double(Lst.P)))]));
    n = check(n,isa(info.G,'copvar') && isa(L,'opvar') && rL<1e-6, ...
        sprintf('Hinf_estimator_sop: gamma %.7g; |P L - Z|/|Z| %.1e (stock getObserver %.1e), L vs stock rel %.1e, inverse d %d relrms %.1e, cond R0 %.1e', ...
        gam,rL,rLst,dL,info.inv.d,max(info.inv.relrms),info.inv.condR0));
end
fprintf('Test_copvar_inv: %d checks passed\n',n);
end

function n = check(n,ok,msg)
if ok,  fprintf('  ok   %s\n',msg);   n = n+1;
else,   error('Test_copvar_inv:fail','FAILED %s',msg);
end
end
function r = opres(E,sg)
r = max([maxabs(double(E.P)), pmax(E.Q1,E.var1,E.var2,sg), pmax(E.Q2,E.var1,E.var2,sg), pmax(E.R.R0,E.var1,E.var2,sg), ...
         kmax(E.R.R1,E.var1,E.var2,sg,true), kmax(E.R.R2,E.var1,E.var2,sg,false)]);
end
function m = maxabs(v)
if isempty(v),  m = 0;  else,   m = max(abs(v(:)));   end
end
function m = pmax(p,v1,v2,sg)
m = kmax(p,v1,v2,sg,[]);
end
function d = pdiff(p,q,v1,v2,sg)
d = kmax(p-q,v1,v2,sg,[]);
end
function m = kmax(p,v1,v2,sg,lower)
% sup over the grid of a parameter: p(s) on sg, or p(s,t) on the lower / upper triangle
if isempty(p),  m = 0;  return,     end
if isa(p,'double'),     m = maxabs(p);  return,     end
m = 0;
if ~any(strcmp(p.varname,v2.varname{1}))
    for i = 1:numel(sg),    m = max(m,maxabs(double(subs(p,v1,sg(i)))));   end
else
    for i = 1:numel(sg)
        pi_ = subs(p,v1,sg(i));
        for j = 1:numel(sg)
            if isempty(lower) || (lower && sg(j)<=sg(i)) || (~lower && sg(j)>=sg(i))
                m = max(m,maxabs(double(subs(pi_,v2,sg(j)))));
            end
        end
    end
end
end
