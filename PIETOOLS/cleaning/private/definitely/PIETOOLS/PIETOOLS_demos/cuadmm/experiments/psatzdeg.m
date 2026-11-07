% psatzdeg.m -- CAN the psatz monomial degree be lowered?
%
% PIETOOLS' own convention gives psatz terms a degree one BELOW the base
% (eq_deg_psatz = eq_deg - 1), because the multiplier itself contributes
% degree.  My earlier 2-D grid gave every psatz term the FULL eq_deg, which
% inflates m.  Sweep the psatz degree offset and measure both m and reach, so
% the cost of the linear-generator fix is not overstated.
%
% checkdeg_lpi_eq_2d is also run as a diagnostic: it is the only routine that
% tests monomial sufficiency, and PIETOOLS_stability_2D ships with it DISABLED
% (toggle = 0).
cuadmm_path;
LAMSTAR = 2*pi^2;
FR = [0.90 0.99 0.999 0.9999];
OFFS = [0 -1 -2];

fprintf('PD offset|m|nblk|Ns|frac|cert|rel_b|feasratio|numerr|checkdeg|t\n');
for off = OFFS
    best = 0; Sm = [];
    for f = FR
        lam = f*LAMSTAR;
        try
            t0=tic; [prog,okdeg] = build4(lam,[3 4 5 6],off);
            so.solver='mosek'; evalc('sol = lpisolve(prog,so);'); tw=toc(t0);
            S = cuadmm_private('sdpshape',sol); if isempty(Sm), Sm=S; end
            I = sol.solinfo.info; xv = sol.solinfo.RRx(:);
            Atf=[];bf=[];
            for i=1:sol.expr.num, Atf=[Atf,sol.expr.At{i}]; bf=[bf;sol.expr.b{i}]; end
            rb = norm(full(Atf'*xv-bf))/norm(full(bf));
            triv = (abs(rb-1)<=1e-6)||norm(xv)<=1e-12;
            cert = (I.numerr==0)&&(I.feasratio>0.5)&&(rb<1e-6)&&~triv;
            fprintf('PD %+d|%d|%d|[%s]|%.4f|%d|%.3e|%+.4f|%d|%d|%.1f\n', ...
                off,S.m,S.nblk,strtrim(num2str(S.Ks)),f,cert,rb,I.feasratio, ...
                I.numerr,okdeg,tw);
            if cert, best=f; else, break; end
        catch ME
            fprintf('PD %+d|ERR|%.4f|%s\n',off,f,strrep(ME.message,newline,' ')); break
        end
    end
    if isempty(Sm), Sm.m=-1; Sm.nblk=-1; end
    fprintf('PDREACH offset=%+d|reach=%.4f|m=%d|nblk=%d\n',off,best,Sm.m,Sm.nblk);
end
fprintf('PDDONE\n');

function [prog,okdeg] = build4(lam,codes,psatz_off)
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
eqd.dx = {1+dx{1}; 1+dx{2}; 1+dx{3}};
eqd.dy = {1+dy{1}, 1+dy{2}, 1+dy{3}};
eqd.d2 = cellfun(@(c) 1+c, d2, 'UniformOutput', false);
pd = bump(eqd,psatz_off);

prog = lpiprogram(PIE.vars(:,1),PIE.vars(:,2),PIE.dom);
np = Top.dim(:,1);
Iop = opvar2d(eye(sum(np)),[np np],PIE.dom,PIE.vars);
Qop = clean_opvar((Aop'*(Iop*Top))' + Aop'*(Iop*Top),1e-12);
eq_opts = get_eq_opts_2D(Qop,s2d.eq_opts,1e-12);
[prog,Qeop] = poslpivar_2d(prog,Qop.dim,eqd,eq_opts);
okdeg = checkdeg_lpi_eq_2d(Qop,Qeop,eqd);      % the DISABLED-by-default check
for p = codes
    o = eq_opts; o.psatz = p;
    [prog,Qe] = poslpivar_2d(prog,Qop.dim,pd,o);
    Qeop = Qeop + Qe;
end
prog = lpi_eq_2d(prog,Qop+Qeop,'symmetric');
end

function d = bump(d,k)
% Shift every degree by k, floored at 0.
d.dx = cellfun(@(c) max(c+k,0), d.dx, 'UniformOutput', false);
d.dy = cellfun(@(c) max(c+k,0), d.dy, 'UniformOutput', false);
d.d2 = cellfun(@(c) max(c+k,0), d.d2, 'UniformOutput', false);
end
