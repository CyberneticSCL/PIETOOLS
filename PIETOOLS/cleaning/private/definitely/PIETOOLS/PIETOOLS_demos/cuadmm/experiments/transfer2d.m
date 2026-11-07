% transfer2d.m -- does the tight-Poincare generator set survive when the
% Lyapunov operator is SEARCHED rather than pinned to the identity?
%
% Maintainer's plan: get the Poincare inequality tight first (done -- the four
% LINEAR box generators {s1,1-s1,s2,1-s2} reach 0.9999*lam* at m=3728 with
% P=I), then check the same certificate works for the heat equation.
%
% Here P is a genuine poslpivar_2d decision variable, exactly as
% PIETOOLS_stability_2D builds it, and only the NEGATIVITY generator set
% changes.  If reach stays ~0.9999 the fix transfers; if it drops, the
% Lyapunov search interacts with the generators and that is worth knowing.
cuadmm_path;
HERE = fileparts(mfilename('fullpath'));
cuadmm_shadow('poslpivar');

LAMSTAR = 2*pi^2;
FR = [0.25 0.50 0.90 0.99 0.999 0.9999];
A = @(x) x;  B = @(x) 1-x;

SETS = { 'P searched, stock psatz OFF' , {}
         'P searched, linear4'         , { @(x,y) A(x), @(x,y) B(x), ...
                                           @(x,y) A(y), @(x,y) B(y) }
         'P=I,       linear4 (ref)'    , { @(x,y) A(x), @(x,y) B(x), ...
                                           @(x,y) A(y), @(x,y) B(y) } };
PFIX = [false false true];

fprintf('TR set|m|nblk|Ns|frac|lam|cert|rel_b|feasratio|numerr|t\n');
for k = 1:size(SETS,1)
    lab = SETS{k,1}; gl = SETS{k,2}; best = 0; Sm = [];
    for f = FR
        lam = f*LAMSTAR;
        try
            t0=tic; prog = build3(lam,gl,PFIX(k));
            so.solver='mosek'; evalc('sol = lpisolve(prog,so);'); tw=toc(t0);
            S = cuadmm_private('sdpshape',sol); if isempty(Sm), Sm=S; end
            I = sol.solinfo.info; xv = sol.solinfo.RRx(:);
            Atf=[];bf=[];
            for i=1:sol.expr.num, Atf=[Atf,sol.expr.At{i}]; bf=[bf;sol.expr.b{i}]; end
            rb = norm(full(Atf'*xv-bf))/norm(full(bf));
            triv = (abs(rb-1)<=1e-6) || norm(xv)<=1e-12;
            cert = (I.numerr==0) && (I.feasratio>0.5) && (rb<1e-6) && ~triv;
            fprintf('TR %s|%d|%d|[%s]|%.4f|%.4f|%d|%.3e|%+.4f|%d|%.1f\n', ...
               lab,S.m,S.nblk,strtrim(num2str(S.Ks)),f,lam,cert,rb,I.feasratio,I.numerr,tw);
            if cert, best=f; else, break; end
        catch ME
            fprintf('TR %s|ERR|%.4f|%s\n',lab,f,strrep(ME.message,newline,' ')); break
        end
    end
    if isempty(Sm), Sm.m=-1; Sm.nblk=-1; end
    fprintf('TRREACH %s|reach=%.4f|m=%d|nblk=%d\n',lab,best,Sm.m,Sm.nblk);
end
fprintf('TRDONE\n');

cuadmm_shadow('off');   % remove the shadow: it stays active for the whole session otherwise

function prog = build3(lam,glist,fixP)
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
if fixP
    Pop = opvar2d(eye(sum(np)),[np np],PIE.dom,PIE.vars);
else
    [prog,Pop] = poslpivar_2d(prog,Top.dim,s2d.LF_deg,s2d.LF_opts);
    Ip  = blkdiag(s2d.eppos(1)*eye(np(1)),s2d.eppos(2)*eye(np(2)), ...
                  s2d.eppos(3)*eye(np(3)),s2d.eppos(4)*eye(np(4)));
    Pop = Pop + opvar2d(Ip,Pop.dim,PIE.dom,PIE.vars);
end
PTop = Pop*Top;  APTop = Aop'*PTop;
Qop  = clean_opvar(APTop' + APTop, 1e-12);

eq_opts = get_eq_opts_2D(Qop,s2d.eq_opts,1e-12);
[prog,Qeop] = poslpivar_2d(prog,Qop.dim,s2d.eq_deg,eq_opts);
for i = 1:numel(glist)
    o = eq_opts; o.psatz = glist{i};
    [prog,Qe] = poslpivar_2d(prog,Qop.dim,s2d.eq_deg,o);
    Qeop = Qeop + Qe;
end
prog = lpi_eq_2d(prog,Qop+Qeop,'symmetric');
end
