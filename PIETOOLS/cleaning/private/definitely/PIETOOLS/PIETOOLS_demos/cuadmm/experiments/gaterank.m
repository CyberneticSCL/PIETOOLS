% gaterank.m -- the WITHIN-PROGRAM ranking test.
%
% gatexfer varied the program (f) at a fixed point, so its ordering claim is
% across programs and weaker than the peer's.  Here ONE program (f_dst, known
% infeasible by Poincare) is scored at SEVERAL distinct non-trivial candidate
% points, obtained by solving at different feasible f_src.  That is the exact
% configuration in which a gate either ranks correctly or inverts: one linear
% system, several x, two measures.
cuadmm_path;
LAMSTAR = pi^2;
SRC = [0.30 0.50 0.70 0.80 0.90];
DST = [1.05 1.20];
st = lpisettings('heavy'); st.sos_opts.solver = 'mosek'; st.eppos=1e-2; st.eppos2=1e-2;

X = cell(size(SRC));
for k = 1:numel(SRC)
    p = build_fixedP1(SRC(k)*LAMSTAR,st);
    evalc('sk = lpisolve(p,st.sos_opts);');
    X{k} = sk.solinfo.RRx(:);
end
for f = DST
    [pd,Dd,Ed,Pd] = build_fixedP1(f*LAMSTAR,st);
    Atf=[];bf=[];
    for i=1:pd.expr.num, Atf=[Atf,pd.expr.At{i}]; bf=[bf;pd.expr.b{i}]; end
    nb = norm(full(bf));
    fprintf('GR f_dst=%.2f src|normx|rel_b|op_ind|op_abs|op_over_row|rank_row|rank_op\n',f);
    rb=zeros(1,numel(SRC)); oi=rb; oa=rb;
    for k = 1:numel(SRC)
        pd.solinfo.RRx = X{k};
        rb(k) = norm(full(Atf'*X{k}-bf))/nb;
        Ds = getsol_lpivar(pd,Dd);  Es = getsol_lpivar(pd,Ed);
        GR = cuadmm_private('opnorm_pi',Ds+Es,16);  GD = cuadmm_private('opnorm_pi',Ds,16);
        oa(k) = norm(GR);  oi(k) = oa(k)/max(norm(GD),eps);
    end
    [~,ir] = sort(rb);  [~,io] = sort(oi);
    rr = zeros(1,numel(SRC)); ro = rr;  rr(ir)=1:numel(SRC); ro(io)=1:numel(SRC);
    for k = 1:numel(SRC)
        fprintf('GR %.2f|%.2f|%.4e|%.4e|%.4e|%.4e|%.3f|%d|%d\n', ...
            f,SRC(k),norm(X{k}),rb(k),oi(k),oa(k),oi(k)/max(rb(k),eps),rr(k),ro(k));
    end
    fprintf('GRINV f=%.2f|kendall_disagree=%d of %d pairs\n', f, ...
        nnz(triu(sign(rb-rb.')~=sign(oi-oi.'),1)), nchoosek(numel(SRC),2));
end
fprintf('GRDONE\n');

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
