% t11_dup2d.m -- does the 98x-cheaper Dup=1 still CERTIFY in 2-D?
% The sep=1 result is the cautionary precedent: a cheaper shape is worthless if
% the restricted class cannot prove anything.  Algebraic + cone measures only --
% the operator-level instrument is 1-D (opvar); an opvar2d version is separate
% work and is NOT claimed here.
cuadmm_path;
HERE = fileparts(mfilename('fullpath'));
cuadmm_shadow('lpisolve');
global CENSUS_PROG
cuadmm_shadow('off');   % real solves from here on

pvar s1 s2 t
x = pde_var('state',1,[s1;s2],[0,1;0,1]);
PIE = convert([diff(x,t,1)==diff(x,s1,2)+diff(x,s2,2)+2*x;
               subs(x,s1,0)==0; subs(x,s1,1)==0;
               subs(x,s2,0)==0; subs(x,s2,1)==0]);

fprintf('D2 eppos|Dup|m|Nmax|eigcost|pinf|dinf|numerr|feasratio|rel_b|rel_cu|psd_min|psd_relmin|trivial|t_solve\n');
for epv = [1e-6 1e-2 1]
for Dup = [1 2 3]
    st = lpisettings('light');
    if Dup ~= 3, st = setDup2(st,Dup); end
    st.settings_2d.eppos = epv*[1;1;1;1];  st.eppos = epv;  st.eppos2 = epv;
    try
        t0 = tic;
        evalc('prog = PIETOOLS_stability_2D(PIE,st);');
        tw = toc(t0);
        S = cuadmm_private('sdpshape',prog);
        I = prog.solinfo.info;
        xv = prog.solinfo.RRx(:);
        Atf=[];bf=[];
        for i=1:prog.expr.num, Atf=[Atf,prog.expr.At{i}]; bf=[bf;prog.expr.b{i}]; end
        res = norm(full(Atf'*xv - bf));  nb = norm(full(bf));
        pmin=inf; pmax=-inf; off=S.Kf;
        for k=1:numel(S.Ks)
            N=S.Ks(k); Xk=reshape(xv(off+(1:N^2)),N,N); Xk=(Xk+Xk')/2;
            ev=eig(Xk); pmin=min(pmin,min(ev)); pmax=max(pmax,max(ev)); off=off+N^2;
        end
        triv = (abs(res/nb - 1) <= 1e-6) || norm(xv)<=1e-12;
        fprintf('D2 %.0e|%d|%d|%d|%g|%d|%d|%d|%+.4f|%.3e|%.3e|%+.3e|%+.3e|%d|%.1f\n', ...
            epv,Dup,S.m,S.Nmax,S.eigcost,I.pinf,I.dinf,I.numerr,I.feasratio, ...
            res/nb, res/(1+nb), pmin, pmin/max(pmax,eps), triv, tw);
    catch ME
        fprintf('D2 %d|ERROR|%s\n',Dup,strrep(ME.message,newline,' '));
    end
end
end
fprintf('D2DONE\n');

function st = setDup2(st,Dup)
s2 = st.settings_2d;
dx = s2.LF_deg.dx;  dy = s2.LF_deg.dy;  d2 = s2.LF_deg.d2;
s2.eq_deg.dx = {Dup+dx{1};  Dup+dx{2};  Dup+dx{3}};
s2.eq_deg.dy = {Dup+dy{1},  Dup+dy{2},  Dup+dy{3}};
s2.eq_deg.d2 = cellfun(@(c) Dup+c, d2, 'UniformOutput', false);
if isfield(s2,'eq_deg_psatz')
    for j = 1:numel(s2.eq_deg_psatz)
        s2.eq_deg_psatz{j}.dx = {Dup-1+dx{1};  Dup-1+dx{2};  Dup-1+dx{3}};
        s2.eq_deg_psatz{j}.dy = {Dup-1+dy{1},  Dup-1+dy{2},  Dup-1+dy{3}};
        s2.eq_deg_psatz{j}.d2 = cellfun(@(c) Dup-1+c, d2, 'UniformOutput', false);
    end
end
st.settings_2d = s2;
end
