function pf_check_fam(storage,solver,which,cfg,dw)
% PF_CHECK_FAM(STORAGE,SOLVER,WHICH,CFG,DW) solves, at light and dw DW
% (default 1) with the storage STORAGE and the solver SOLVER (default
% mosek), the term set CFG (a PF_CFG element; [] = the tier's own 'auto'),
% the LPIs of the systems WHICH (default all three)
% ('stock' or 'def'), the stability LPI of D4 (S17 with a Dirichlet 1-D heat
% state as the edge data; the target x2 of order 0 in s1; r2 = 0.9 pi^2,
% the record's certified case) and of S3 (control: x2 diffuses in s1 too;
% r = 0.9 pi^2), and evaluates each certified derivative operator
% D = T'*P*A + A'*P*T on the two-direction family of the obstruction
% argument: e1 = (L_n(s1), 0), e2 = (0, L_n(s1) psi''(s2)), L_n the
% shifted Legendre polynomial of unit norm, psi = s2 (1-s2)^k. With
% g(c) = <e1 + c e2, D (e1 + c e2)> = a + 2 b c + d c^2, the certificate
% -D >= 0 requires a <= 0, d <= 0 and b^2 <= a d; b^2 > a d refutes it. The
% obstruction argument predicts a -> 0 like n^-4 and b bounded away from 0
% when x2 has no s1-diffusion (D4). V = <T e, P T e> is printed for scale.
% SCRATCH, MMP, 11/16/2026.
if nargin<1 || isempty(storage),    storage = 'stock';  end
if nargin<2 || isempty(solver),     solver = 'mosek';   end
if nargin<3 || isempty(which),      which = 1:3;        end
if nargin<4,                        cfg = [];           end
if nargin<5 || isempty(dw),         dw = 1;             end
cname = 'auto';     if ~isempty(cfg),   cname = cfg.name;   end
if strcmp(storage,'def'),   pf_setup('route');   else,   pf_setup(false);  end
pvar s1 s2 t
r2 = 0.9*pi^2;
names = {'D4','S3','S17'};    pdes = cell(1,3);
clear stateNameGenerator
x1 = pde_var(1,s1,[0,1]);   x2 = pde_var(1,[s1;s2],[0,1;0,1],[0;2]);
pdes{1} = [diff(x1,t)==diff(x1,s1,2); diff(x2,t)==diff(x2,s2,2)+r2*x2; ...
           subs(x1,s1,0)==0; subs(x1,s1,1)==0; subs(x2,s2,0)==x1; subs(x2,s2,1)==0];
Sys = mixed_domain_systems('S3_edge1d_heat2d',struct('r',0.9*pi^2));    pdes{2} = Sys.PDE;
Sys = mixed_domain_systems('S17_transport_edge_diffusion',struct('r2',0.5*pi^2));    pdes{3} = Sys.PDE;
for i = which
    PIE = initialize(convert(pdes{i},{'pie','silent'}));
    st = lpisettings_sop('light');  st.sos_opts.solver = solver;    st.sop.verbose = false;
    st.sop.dw_feas = dw*[1 1 1];    st.sop.max_raise_feas = [0 0 0];
    if ~isempty(cfg)
        st.lpi.ineq.psatz = 0;  st.lpi.ineq.psatz_offset = 0;
        st.lpi.ineq.psatz_rows = cfg.rows;  st.lpi.ineq.psatz_rows_offset = cfg.offs;
        st.lpi.ineq.psatz_rows_perdir = cfg.pd;
    end
    t0 = tic;   [~,P,info] = PIETOOLS_stability_sop(PIE,st,'P');    ts = toc(t0);
    fprintf('@@ CHK %s storage=%s solver=%s cfg=%s dw=%d %s cert=%d ndv=%d m=%d t=%.0f\n',names{i},storage,solver,cname,dw,info.status,info.certified,info.shape.ndv,info.shape.m,ts);
    if isempty(P),  continue;   end
    PIEc = exec_pie_sop(PIE);   T = PIEc.T;     A = PIEc.A;
    D = T'*P*A + A'*P*T;
    [sp,~] = copvar_space_list(T,'in');     % the fundamental-state spaces, in the container's order
    k1 = find(cellfun(@(c) isequal(sort(c),{'s1'}),sp));
    k2 = find(cellfun(@(c) isequal(sort(c),{'s1','s2'}),sp));
    assert(numel(sp)==2 && ~isempty(k1) && ~isempty(k2),'pf_check_d4: unexpected spaces');
    for n = [4 6 8]
        Ln = shifted_legendre(n,s1);
        for k = [1 3]
            psi = s2*(1-s2)^k;    psi2 = diff(diff(psi,s2),s2);
            e1 = cell(1,2);     e2 = cell(1,2);
            e1{k1} = Ln;        e1{k2} = polynomial(0);
            e2{k1} = polynomial(0);     e2{k2} = Ln*psi2;
            y1 = apply_cop(D,e1);   y2 = apply_cop(D,e2);
            a = ip(e1,y1,k1,k2);    b = ip(e1,y2,k1,k2);    b2 = ip(e2,y1,k1,k2);   d = ip(e2,y2,k1,k2);
            w1 = apply_cop(T,e1);   w2 = apply_cop(T,e2);
            va = ip(w1,apply_cop(P,w1),k1,k2);  vd = ip(w2,apply_cop(P,w2),k1,k2);
            cstar = NaN;    gmax = NaN;
            if d<0,     cstar = -b/d;   gmax = a + 2*b*cstar + d*cstar^2;    end
            fprintf('@@ FAM %s storage=%s solver=%s n=%d k=%d a=%.3e b=%.3e (b''=%.3e) d=%.3e b^2-ad=%.3e Va=%.3e Vd=%.3e gmax=%.3e at c=%.3g %s\n', ...
                    names{i},storage,solver,n,k,a,b,b2,d,b^2-a*d,va,vd,gmax,cstar,ternary(b^2>a*d || a>0 || d>0,'VIOLATED','ok'));
        end
    end
    % Family B (Step 1 of the argument): a boundary layer in s2 fed by the
    % edge datum, x1 = g, x2 = g(s1) (1-s2)^K; the fundamental state by the
    % orders of each system (D4: x1 order 2, x2 order (0,2); S3: 2, (2,2);
    % S17: 1, (0,2)); g = s1 L_m (S17) or s1 (1-s1) L_m (D4, S3). The
    % certified form f = <xh, -D xh> must be >= 0; Terr checks T xh = (g, x2).
    o1 = [2 2 1];   o2 = [0 2; 2 2; 0 2];
    for K = [4 8 16]
        for m = [0 2 4 8]
            Lm = shifted_legendre(m,s1);
            if i==3,    g = s1*Lm;  else,   g = s1*(1-s1)*Lm;  end
            X2 = g*(1-s2)^K;
            xh1 = g;    for j = 1:o1(i),    xh1 = diff(xh1,s1);  end
            xh2 = X2;   for j = 1:o2(i,1),  xh2 = diff(xh2,s1);  end
            for j = 1:o2(i,2),  xh2 = diff(xh2,s2);  end
            e = cell(1,2);  e{k1} = xh1;    e{k2} = xh2;
            f = -ip(e,apply_cop(D,e),k1,k2);
            w = apply_cop(T,e);     v = ip(w,apply_cop(P,w),k1,k2);
            dc = cell(1,2);     dc{k1} = polynomial(w{k1})-g;     dc{k2} = polynomial(w{k2})-X2;
            terr = sqrt(abs(ip(dc,dc,k1,k2)));
            fprintf('@@ FAMB %s storage=%s solver=%s K=%d m=%d f=%.3e V=%.3e f/V=%.3e Terr=%.1e %s\n', ...
                    names{i},storage,solver,K,m,f,v,f/v,terr,ternary(f<-1e-9,'VIOLATED','ok'));
        end
    end
end
disp('## CHK DONE');
end

function L = shifted_legendre(n,s)
% L = SHIFTED_LEGENDRE(N,S) the shifted Legendre polynomial of degree N on
% [0,1] in the variable S, scaled to unit L2 norm.
u = 2*s-1;  L0 = polynomial(1);  L1 = u;
if n==0,    L = L0*sqrt(1);     return;     end
for m = 1:n-1
    L2 = ((2*m+1)*u*L1 - m*L0)/(m+1);   L0 = L1;   L1 = L2;
end
L = L1*sqrt(2*n+1);
end

function y = apply_cop(P,v)
% Y = APPLY_COP(P,V) applies the 'copvar' P to the cell V of 'polynomial'
% blocks: y{i} = sum_j P.C{i,j} v{j} (apply_sopvar), as test_mixed_domain_pde.
y = cell(size(P.C,1),1);
for i = 1:size(P.C,1)
    acc = polynomial(zeros(P.dim_out(i),1));
    for j = 1:size(P.C,2)
        if isempty(P.C{i,j}),   continue,   end
        acc = acc + apply_sopvar(P.C{i,j},polynomial(v{j}));
    end
    y{i} = acc;
end
end

function v = ip(x,y,k1,k2)
% V = IP(X,Y,K1,K2) the L2 inner product of the state cells X and Y: block
% K1 over s1, block K2 over (s1,s2), each on [0,1].
p1 = polynomial(x{k1})*polynomial(y{k1});     p2 = polynomial(x{k2})*polynomial(y{k2});
v = double(int(p1,polynomial({'s1'}),0,1)) + double(int(int(p2,polynomial({'s2'}),0,1),polynomial({'s1'}),0,1));
end

function s = ternary(c,a,b)
if c,   s = a;  else,   s = b;  end
end
