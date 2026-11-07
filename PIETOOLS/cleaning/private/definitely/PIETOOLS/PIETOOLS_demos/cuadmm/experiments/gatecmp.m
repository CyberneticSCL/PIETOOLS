% gatecmp.m -- both residual quantities at points where INFEASIBILITY IS ANALYTIC.
%
% The 1-D Dirichlet heat PIE with P fixed to the identity: d/dt<x,x> <= 0 holds
% iff lam <= lam* = pi^2 (Poincare).  So f = lam/lam* > 1 is infeasible by
% theorem, not by a solver verdict -- which is what makes this a valid test of
% a GATE.  Reported side by side, on the same solved point:
%    rel_b   = ||At'x-b||/||b||             the SDP row residual
%    op_ind  = ||Dop+Deop|| / ||Dop||       operator, solution-dependent denom
%    op_indP = ||Dop+Deop|| / ||Pop||       operator, non-collapsing denom
% If a gate on the operator quantity ACCEPTS an f>1 row, the gate is wrong.
cuadmm_path;
LAMSTAR = pi^2;
FR = [0.50 0.90 0.99 1.01 1.05 1.20 2.00];
fprintf('GC f|lam|m|rel_b|op_ind|op_indP|psd_min|normDop|normx|feasratio|numerr|t\n');
for f = FR
    lam = f*LAMSTAR;
    st = lpisettings('heavy');  st.sos_opts.solver = 'mosek';
    st.eppos = 1e-2; st.eppos2 = 1e-2;
    t0 = tic;
    [prog,Dop,Deop,Pop] = build_fixedP1(lam,st);
    evalc('sol = lpisolve(prog,st.sos_opts);');
    tw = toc(t0);
    Atf=[];bf=[];
    for i=1:sol.expr.num, Atf=[Atf,sol.expr.At{i}]; bf=[bf;sol.expr.b{i}]; end
    xv = sol.solinfo.RRx(:);
    rb = norm(full(Atf'*xv-bf))/norm(full(bf));
    Ds = getsol_lpivar(sol,Dop);  Es = getsol_lpivar(sol,Deop);
    GR = cuadmm_private('opnorm_pi',Ds+Es,16);  GD = cuadmm_private('opnorm_pi',Ds,16);  GP = cuadmm_private('opnorm_pi',Pop,16);
    S = cuadmm_private('sdpshape',sol);  off = S.Kf;  pmin = inf;
    for k = 1:numel(S.Ks)
        N = S.Ks(k);  Xk = reshape(xv(off+(1:N^2)),N,N);
        pmin = min(pmin, min(eig((Xk+Xk')/2)));  off = off + N^2;
    end
    I = sol.solinfo.info;
    fprintf('GC %.2f|%.4f|%d|%.4e|%.4e|%.4e|%+.3e|%.4e|%.4e|%+.4f|%d|%.1f\n', ...
        f,lam,S.m,rb,norm(GR)/max(norm(GD),eps),norm(GR)/max(norm(GP),eps), ...
        pmin,norm(GD),norm(xv),I.feasratio,I.numerr,tw);
end
fprintf('GCDONE\n');

function [prog,Dop,Deop,Pop] = build_fixedP1(lam,st)
pvar s t
x   = pde_var('state',1,s,[0,1]);
PIE = convert([diff(x,t,1)==diff(x,s,2)+lam*x; subs(x,s,0)==0; subs(x,s,1)==0]);
PIE = initialize(PIE);
Top = PIE.T;  Aop = PIE.A;
prog = lpiprogram(PIE.vars(:,1),PIE.vars(:,2),PIE.dom);
Iop = opvar();  Iop.I = PIE.dom;  Iop.var1 = PIE.vars(1,1); Iop.var2 = PIE.vars(1,2);
Iop.dim = Top.dim;
n0 = Top.dim(1,1);  n1 = Top.dim(2,1);
Iop.P = eye(n0);  Iop.R.R0 = eye(n1);
Pop = Iop;                                   % P = I, the operator the gate scales by
Dop = Aop'*(Iop*Top) + (Iop*Top)'*Aop;
[prog,De1op] = poslpivar(prog,Dop.dim,st.dd2,st.options2);
Deop = De1op;
if st.override2 ~= 1
    [prog,De2op] = poslpivar(prog,Dop.dim,st.dd3,st.options3);
    Deop = Deop + De2op;
end
prog = lpi_eq(prog,Dop+Deop,'symmetric');
end
