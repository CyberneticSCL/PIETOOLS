% poinc2d.m -- transfer the 1-D finding to 2-D, with P FIXED TO THE IDENTITY.
%
% 1-D result (measured, P=I): the Positivstellensatz term on the NEGATIVITY
% operator is REQUIRED -- without it nothing certifies even at 0.5*lam*, at any
% degree up to Dup3 -- and with it, Dup=2 certifies at 0.9999*lam*.
%
% The prediction this makes for 2-D: in 1-D psatz on the negativity is ON by
% default (override2 = 0), but in 2-D `eq_use_psatz = [0;0]`, i.e. OFF.  That
% would explain why 2-D stalls at 0.25*lam* while 1-D reaches 0.9999.
%
% lam* = 2*pi^2 = 19.7392 for the unit square with Dirichlet on all four edges.

cuadmm_path;
LAMSTAR = 2*pi^2;
FR = [0.10 0.25 0.50 0.75 0.90 0.99];

V = { 'Dup3 psatzOFF(default)' , @(st) st
      'Dup3 psatz1'            , @(st) eqp(st,1)
      'Dup3 psatz12'           , @(st) eqp(st,[1;2])
      'Dup2 psatz12'           , @(st) eqp(dup(st,2),[1;2])
      'Dup1 psatz12'           , @(st) eqp(dup(st,1),[1;2]) };

fprintf('P2 variant|m|nblk|Ns|eigcost|frac|lam|cert|rel_b|feasratio|numerr|t\n');
for v = 1:size(V,1)
    lab = V{v,1};  best = 0;  Sm = [];
    for f = FR
        lam = f*LAMSTAR;
        try
            st = lpisettings('heavy');
            st.settings_2d.eppos = 1e-2*[1;1;1;1]; st.eppos=1e-2; st.eppos2=1e-2;
            st.sos_opts.solver = 'mosek';        % 8.3x faster and far more
                                                 % accurate than SeDuMi here
            st = V{v,2}(st);
            t0 = tic;
            [prog,~] = build_fixedP2(lam,st);
            evalc('sol = lpisolve(prog,st.sos_opts);');
            tw = toc(t0);
            S = cuadmm_private('sdpshape',sol);  if isempty(Sm), Sm = S; end
            I = sol.solinfo.info;  xv = sol.solinfo.RRx(:);
            Atf=[];bf=[];
            for i=1:sol.expr.num, Atf=[Atf,sol.expr.At{i}]; bf=[bf;sol.expr.b{i}]; end
            rb = norm(full(Atf'*xv-bf))/norm(full(bf));
            triv = (abs(rb-1)<=1e-6) || norm(xv)<=1e-12;
            cert = (I.numerr==0) && (I.feasratio>0.5) && (rb<1e-6) && ~triv;
            fprintf('P2 %s|%d|%d|[%s]|%g|%.2f|%.4f|%d|%.3e|%+.4f|%d|%.1f\n', ...
                lab,S.m,S.nblk,strtrim(num2str(S.Ks)),S.eigcost,f,lam,cert,rb, ...
                I.feasratio,I.numerr,tw);
            if cert, best = f; else, break; end
        catch ME
            fprintf('P2 %s|ERR|%.2f|%s\n',lab,f,strrep(ME.message,newline,' '));
            break
        end
    end
    if isempty(Sm), Sm.m = -1; end
    fprintf('P2REACH %s|reach=%.2f|m=%d\n',lab,best,Sm.m);
end
fprintf('P2DONE\n');


function [prog,Qop] = build_fixedP2(lam,st)
% The 2-D stability LPI with the Lyapunov operator PINNED to the identity, so
% the negativity certificate is the only unknown.  Mirrors
% PIETOOLS_stability_2D's STEP 3/4 with Pop replaced by I.
pvar s1 s2 t
x = pde_var('state',1,[s1;s2],[0,1;0,1]);
PIE = convert([diff(x,t,1)==diff(x,s1,2)+diff(x,s2,2)+lam*x;
               subs(x,s1,0)==0; subs(x,s1,1)==0;
               subs(x,s2,0)==0; subs(x,s2,1)==0]);
PIE = initialize(PIE);
Top = PIE.T;  Aop = PIE.A;
s2d = st.settings_2d;

prog = lpiprogram(PIE.vars(:,1),PIE.vars(:,2),PIE.dom);

% P = I
np  = Top.dim(:,1);
Iop = opvar2d(eye(sum(np)),[np np],PIE.dom,PIE.vars);

PTop  = Iop*Top;
APTop = Aop'*PTop;
Qop   = APTop' + APTop;
Qop   = clean_opvar(Qop,1e-12);

eq_opts = get_eq_opts_2D(Qop,s2d.eq_opts,1e-12);
[prog,Qeop] = poslpivar_2d(prog,Qop.dim,s2d.eq_deg,eq_opts);
for j = 1:numel(s2d.eq_use_psatz)
    if s2d.eq_use_psatz(j) ~= 0
        o = s2d.eq_opts_psatz{j};
        o.exclude = o.exclude | eq_opts.exclude;
        o.sep     = o.sep     | eq_opts.sep;
        [prog,Qe2] = poslpivar_2d(prog,Qop.dim,s2d.eq_deg_psatz{j},o);
        Qeop = Qeop + Qe2;
    end
end
prog = lpi_eq_2d(prog,Qop+Qeop,'symmetric');
end

function st = eqp(st,v)
v = v(:);
st.settings_2d.eq_use_psatz = v;
for j = 1:numel(v), st.settings_2d.eq_opts_psatz{j}.psatz = v(j); end
end

function st = dup(st,Dup)
s2 = st.settings_2d;
dx=s2.LF_deg.dx; dy=s2.LF_deg.dy; d2=s2.LF_deg.d2;
s2.eq_deg.dx = {Dup+dx{1};  Dup+dx{2};  Dup+dx{3}};
s2.eq_deg.dy = {Dup+dy{1},  Dup+dy{2},  Dup+dy{3}};
s2.eq_deg.d2 = cellfun(@(c) Dup+c, d2, 'UniformOutput', false);
for j = 1:numel(s2.eq_deg_psatz)
    s2.eq_deg_psatz{j}.dx = {Dup-1+dx{1}; Dup-1+dx{2}; Dup-1+dx{3}};
    s2.eq_deg_psatz{j}.dy = {Dup-1+dy{1}, Dup-1+dy{2}, Dup-1+dy{3}};
    s2.eq_deg_psatz{j}.d2 = cellfun(@(c) Dup-1+c, d2, 'UniformOutput', false);
end
st.settings_2d = s2;
end
