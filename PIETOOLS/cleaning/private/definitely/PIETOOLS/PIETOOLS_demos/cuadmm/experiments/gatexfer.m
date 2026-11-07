% gatexfer.m -- the ranking test needs a NON-TRIVIAL infeasible point.
%
% gatecmp showed that solving directly at lam > lam* returns X=0, where both
% measures agree by construction and neither ranks anything.  To get a point
% that is infeasible but NOT trivial, take the certificate solved at f_src
% (feasible) and evaluate it against the program built at f_dst > 1, where
% Poincare says no certificate exists.  Same dimensions, same structure, only
% lam differs, so the vector transplants exactly.
%
% Both quantities are then computed at ONE point that is known-infeasible by
% theorem and known-nonzero by construction -- the configuration in which a
% gate either ranks correctly or inverts.
cuadmm_path;
LAMSTAR = pi^2;
SRC = 0.90;
DST = [0.90 0.95 1.01 1.05 1.20 2.00];
st = lpisettings('heavy'); st.sos_opts.solver = 'mosek'; st.eppos=1e-2; st.eppos2=1e-2;

[p0,D0,E0,P0] = build_fixedP1(SRC*LAMSTAR,st);
evalc('s0 = lpisolve(p0,st.sos_opts);');
x0 = s0.solinfo.RRx(:);
fprintf('GX src f=%.2f normx=%.4e\n',SRC,norm(x0));

fprintf('GX f_dst|rel_b|op_ind|op_indP|normDop|normRes|ratio_op_over_row\n');
for f = DST
    [pd,Dd,Ed,Pd] = build_fixedP1(f*LAMSTAR,st);
    Atf=[];bf=[];
    for i=1:pd.expr.num, Atf=[Atf,pd.expr.At{i}]; bf=[bf;pd.expr.b{i}]; end
    % transplant: give the destination program the source point
    pd.solinfo.RRx = x0;  pd.solinfo.x = x0;
    rb = norm(full(Atf'*x0-bf))/norm(full(bf));
    Ds = getsol_lpivar(pd,Dd);  Es = getsol_lpivar(pd,Ed);
    GR = cuadmm_private('opnorm_pi',Ds+Es,16); GD = cuadmm_private('opnorm_pi',Ds,16); GP = cuadmm_private('opnorm_pi',Pd,16);
    oi = norm(GR)/max(norm(GD),eps);  oip = norm(GR)/max(norm(GP),eps);
    fprintf('GX %.2f|%.4e|%.4e|%.4e|%.4e|%.4e|%.3f\n', ...
        f,rb,oi,oip,norm(GD),norm(GR),oi/max(rb,eps));
end
fprintf('GXDONE\n');

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
