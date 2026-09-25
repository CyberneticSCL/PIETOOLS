% psatz2d.m -- WHICH Positivstellensatz generators does the 2-D negativity
% test actually need?  P is pinned to the identity, so the negativity
% certificate is the only unknown.
%
% [0,1]^2 is cut out by TWO polynomials, g1 = s1(1-s1) and g2 = s2(1-s2).
% Stock PIETOOLS offers only psatz=1 (the PRODUCT g1*g2) and psatz=2 (a ball).
% A shadowed poslpivar_2d adds psatz=3 (g1 alone) and psatz=4 (g2 alone), so
% the Putinar certificate  S0 + g1*S1 + g2*S2  and the Schmudgen extension
% (+ g1*g2*S3) can be tested for the first time.
%
% Run at Dup=1 so each solve is ~20 s: the 1-D result says the GENERATOR SET
% matters more than the degree, and that is what this isolates.
cuadmm_path;
HERE = fileparts(mfilename('fullpath'));
cuadmm_shadow('poslpivar');
fprintf('PZ poslpivar_2d resolves to %s\n', which('poslpivar_2d'));

LAMSTAR = 2*pi^2;
FR = [0.10 0.25 0.50 0.75 0.90 0.99];

V = { 'none (stock default)'   , []
      'g1g2 (stock psatz1)'    , 1
      'ball (stock psatz2)'    , 2
      'g1 alone (NEW)'         , 3
      'g2 alone (NEW)'         , 4
      'g1,g2 (Putinar, NEW)'   , [3 4]
      'g1,g2,g1g2 (Schmudgen)' , [3 4 1] };

fprintf('PZ variant|m|nblk|Ns|frac|lam|cert|rel_b|feasratio|numerr|t\n');
for v = 1:size(V,1)
    lab = V{v,1};  pl = V{v,2};  best = 0;  Sm = [];
    for f = FR
        lam = f*LAMSTAR;
        try
            t0 = tic;
            prog = build2(lam,pl);
            so.solver='mosek';
            evalc('sol = lpisolve(prog,so);');
            tw = toc(t0);
            S = cuadmm_private('sdpshape',sol); if isempty(Sm), Sm=S; end
            I = sol.solinfo.info; xv = sol.solinfo.RRx(:);
            Atf=[];bf=[];
            for i=1:sol.expr.num, Atf=[Atf,sol.expr.At{i}]; bf=[bf;sol.expr.b{i}]; end
            rb = norm(full(Atf'*xv-bf))/norm(full(bf));
            triv = (abs(rb-1)<=1e-6) || norm(xv)<=1e-12;
            cert = (I.numerr==0) && (I.feasratio>0.5) && (rb<1e-6) && ~triv;
            fprintf('PZ %s|%d|%d|[%s]|%.2f|%.4f|%d|%.3e|%+.4f|%d|%.1f\n', ...
                lab,S.m,S.nblk,strtrim(num2str(S.Ks)),f,lam,cert,rb,I.feasratio,I.numerr,tw);
            if cert, best=f; else, break; end
        catch ME
            fprintf('PZ %s|ERR|%.2f|%s\n',lab,f,strrep(ME.message,newline,' '));
            break
        end
    end
    if isempty(Sm), Sm.m=-1; Sm.Ks=[]; end
    fprintf('PZREACH %s|reach=%.2f|m=%d|Ns=[%s]\n',lab,best,Sm.m,strtrim(num2str(Sm.Ks)));
end
fprintf('PZDONE\n');

cuadmm_shadow('off');   % remove the shadow: it stays active for the whole session otherwise

function prog = build2(lam,psatz_list)
pvar s1 s2 t
x = pde_var('state',1,[s1;s2],[0,1;0,1]);
PIE = convert([diff(x,t,1)==diff(x,s1,2)+diff(x,s2,2)+lam*x;
               subs(x,s1,0)==0; subs(x,s1,1)==0;
               subs(x,s2,0)==0; subs(x,s2,1)==0]);
PIE = initialize(PIE);
Top = PIE.T; Aop = PIE.A;

st = lpisettings('heavy');
st.settings_2d.eppos = 1e-2*[1;1;1;1];
s2d = st.settings_2d;
% Dup = 1 on the negativity degrees, to keep each solve ~20 s.
dx=s2d.LF_deg.dx; dy=s2d.LF_deg.dy; d2=s2d.LF_deg.d2;
s2d.eq_deg.dx = {1+dx{1}; 1+dx{2}; 1+dx{3}};
s2d.eq_deg.dy = {1+dy{1}, 1+dy{2}, 1+dy{3}};
s2d.eq_deg.d2 = cellfun(@(c) 1+c, d2, 'UniformOutput', false);

prog = lpiprogram(PIE.vars(:,1),PIE.vars(:,2),PIE.dom);
np  = Top.dim(:,1);
Iop = opvar2d(eye(sum(np)),[np np],PIE.dom,PIE.vars);
Qop = (Aop'*(Iop*Top))' + Aop'*(Iop*Top);
Qop = clean_opvar(Qop,1e-12);

eq_opts = get_eq_opts_2D(Qop,s2d.eq_opts,1e-12);
[prog,Qeop] = poslpivar_2d(prog,Qop.dim,s2d.eq_deg,eq_opts);
for p = psatz_list(:)'
    o = eq_opts;  o.psatz = p;
    [prog,Qe] = poslpivar_2d(prog,Qop.dim,s2d.eq_deg,o);
    Qeop = Qeop + Qe;
end
prog = lpi_eq_2d(prog,Qop+Qeop,'symmetric');
end
