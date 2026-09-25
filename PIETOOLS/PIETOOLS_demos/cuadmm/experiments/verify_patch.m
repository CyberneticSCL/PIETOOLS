% verify_patch.m -- the repo patch must reproduce the shadowed function-handle
% result exactly.  NO shadow on the path: this exercises the committed code.
cuadmm_path;
fprintf('VP poslpivar_2d resolves to %s\n', which('poslpivar_2d'));
LAMSTAR = 2*pi^2;
for f = [0.99 0.999 0.9999]
    lam = f*LAMSTAR;
    prog = build3(lam,[3 4 5 6]);
    so.solver='mosek'; evalc('sol = lpisolve(prog,so);');
    S = cuadmm_private('sdpshape',sol); I = sol.solinfo.info; xv = sol.solinfo.RRx(:);
    Atf=[];bf=[];
    for i=1:sol.expr.num, Atf=[Atf,sol.expr.At{i}]; bf=[bf;sol.expr.b{i}]; end
    rb = norm(full(Atf'*xv-bf))/norm(full(bf));
    cert = (I.numerr==0)&&(I.feasratio>0.5)&&(rb<1e-6)&&abs(rb-1)>1e-6;
    fprintf('VP frac=%.4f m=%d nblk=%d cert=%d rel_b=%.3e feasratio=%+.4f\n', ...
            f,S.m,S.nblk,cert,rb,I.feasratio);
end
fprintf('VPDONE\n');

function prog = build3(lam,codes)
pvar s1 s2 t
x = pde_var('state',1,[s1;s2],[0,1;0,1]);
PIE = convert([diff(x,t,1)==diff(x,s1,2)+diff(x,s2,2)+lam*x;
               subs(x,s1,0)==0; subs(x,s1,1)==0;
               subs(x,s2,0)==0; subs(x,s2,1)==0]);
PIE = initialize(PIE);
Top = PIE.T; Aop = PIE.A;
st = lpisettings('heavy'); st.settings_2d.eppos = 1e-2*[1;1;1;1];
s2d = st.settings_2d;
dx=s2d.LF_deg.dx; dy=s2d.LF_deg.dy; d2=s2d.LF_deg.d2;
s2d.eq_deg.dx = {1+dx{1}; 1+dx{2}; 1+dx{3}};
s2d.eq_deg.dy = {1+dy{1}, 1+dy{2}, 1+dy{3}};
s2d.eq_deg.d2 = cellfun(@(c) 1+c, d2, 'UniformOutput', false);
prog = lpiprogram(PIE.vars(:,1),PIE.vars(:,2),PIE.dom);
np = Top.dim(:,1);
Iop = opvar2d(eye(sum(np)),[np np],PIE.dom,PIE.vars);
Qop = clean_opvar((Aop'*(Iop*Top))' + Aop'*(Iop*Top),1e-12);
eq_opts = get_eq_opts_2D(Qop,s2d.eq_opts,1e-12);
[prog,Qeop] = poslpivar_2d(prog,Qop.dim,s2d.eq_deg,eq_opts);
for p = codes
    o = eq_opts; o.psatz = p;
    [prog,Qe] = poslpivar_2d(prog,Qop.dim,s2d.eq_deg,o);
    Qeop = Qeop + Qe;
end
prog = lpi_eq_2d(prog,Qop+Qeop,'symmetric');
end
