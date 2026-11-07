% gaterank2.m -- the within-program test with a NON-DEGENERATE candidate set.
%
% gaterank returned op/row EXACTLY constant (6.204, 6.459) and 0/10 rank
% disagreements.  That is an artefact: solving at five different f_src with
% P=I produces five residuals that are PARALLEL, so every degree-1 homogeneous
% measure agrees by construction and the test has no power.  Here the solved
% points are augmented with random perturbations, which break parallelism, and
% with points from a second solver.  Only a candidate set spanning more than
% one direction can rank anything.
cuadmm_path;
rng(7);                                   % fixed seed: these scripts are otherwise seed-flaky
LAMSTAR = pi^2;
SRC = [0.30 0.60 0.90];
DST = 1.20;
st = lpisettings('heavy'); st.eppos=1e-2; st.eppos2=1e-2;
CAND = {};  LBL = {};
for slv = {'mosek','sedumi'}
    for f = SRC
        s2 = st; s2.sos_opts.solver = slv{1};
        p = build_fixedP1(f*LAMSTAR,s2);
        evalc('sk = lpisolve(p,s2.sos_opts);');
        CAND{end+1} = sk.solinfo.RRx(:);  LBL{end+1} = sprintf('%s_f%.2f',slv{1},f);
    end
end
base = CAND{3};
for j = 1:6                                % perturbations at several magnitudes
    r = randn(numel(base),1);
    CAND{end+1} = base + (10^(-j/2))*norm(base)*r/norm(r);
    LBL{end+1} = sprintf('pert1e-%.1f',j/2);
end
[pd,Dd,Ed] = build_fixedP1(DST*LAMSTAR,st);
Atf=[];bf=[];
for i=1:pd.expr.num, Atf=[Atf,pd.expr.At{i}]; bf=[bf;pd.expr.b{i}]; end
nb = norm(full(bf));
n = numel(CAND);  rb=zeros(1,n); oi=rb; oa=rb;
fprintf('G2 f_dst=%.2f cand|normx|rel_b|op_ind|op_abs|op_over_row\n',DST);
for k = 1:n
    pd.solinfo.RRx = CAND{k};
    rb(k) = norm(full(Atf'*CAND{k}-bf))/nb;
    Ds = getsol_lpivar(pd,Dd);  Es = getsol_lpivar(pd,Ed);
    GR = cuadmm_private('opnorm_pi',Ds+Es,16);  GD = cuadmm_private('opnorm_pi',Ds,16);
    oa(k)=norm(GR);  oi(k)=oa(k)/max(norm(GD),eps);
    fprintf('G2 %s|%.4e|%.4e|%.4e|%.4e|%.4f\n',LBL{k},norm(CAND{k}),rb(k),oi(k),oa(k),oi(k)/max(rb(k),eps));
end
d_i = nnz(triu(sign(rb-rb.')~=sign(oi-oi.'),1));
d_a = nnz(triu(sign(rb-rb.')~=sign(oa-oa.'),1));
fprintf('G2SUM pairs=%d|disagree_op_ind=%d|disagree_op_abs=%d|ratio_min=%.4f|ratio_max=%.4f|spread=%.1fx\n', ...
    nchoosek(n,2),d_i,d_a,min(oi./max(rb,eps)),max(oi./max(rb,eps)),max(oi./max(rb,eps))/min(oi./max(rb,eps)));
fprintf('G2DONE\n');

function [prog,Dop,Deop,Pop] = build_fixedP1(lam,st)
pvar s t
x   = pde_var('state',1,s,[0,1]);
PIE = convert([diff(x,t,1)==diff(x,s,2)+lam*x; subs(x,s,0)==0; subs(x,s,1)==0]);
PIE = initialize(PIE);
Top = PIE.T;  Aop = PIE.A;
prog = lpiprogram(PIE.vars(:,1),PIE.vars(:,2),PIE.dom);
Iop = opvar();  Iop.I = PIE.dom;  Iop.var1 = PIE.vars(1,1); Iop.var2 = PIE.vars(1,2);
Iop.dim = Top.dim;
Iop.P = eye(Top.dim(1,1));  Iop.R.R0 = eye(Top.dim(2,1));
Pop = Iop;
Dop = Aop'*(Iop*Top) + (Iop*Top)'*Aop;
[prog,De1op] = poslpivar(prog,Dop.dim,st.dd2,st.options2);
Deop = De1op;
if st.override2 ~= 1
    [prog,De2op] = poslpivar(prog,Dop.dim,st.dd3,st.options3);
    Deop = Deop + De2op;
end
prog = lpi_eq(prog,Dop+Deop,'symmetric');
end
