function R = test_heatNd_lpi()
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% R = TEST_HEATND_LPI() asserting tests of the benchmark LPI pipeline
% (HEATND_LPI -> HEATND_SDP -> HEATND_SOLVE / HEATND_BISECT), N = 1, 2:
%
%  (1) basis counts from the definitions: P 'paper' has mu(d)^N decision
%      variables (mu = d^2+4d+3, paper eq. 13-14), 'tensor' ((d+1)+2(d+1)^2)^N,
%      'paper' without multiplier ((d+1)(d+2))^N; R, Q 'paper' unpruned Gram
%      mu(d')^N; 'Tspan' Gram = product of per-direction counts (DD: 2 mult
%      + 2x4 integral {1,s,th,s th}; DN: 2 + 2x3 {1,s,th}), and pruned R =
%      prod(2 n_int), Q = R + one-multiplier blocks;
%  (2) 'split' (15a) has the same row space as 'full' (numerical rank of
%      the stacked system, and of each) and the same verdicts;
%  (3) the SDP is affine in k: At(k3) from two builds equals a direct build;
%  (4) the dump route (Sedumi2Mosek + mosekopt) and lpisolve agree on the
%      verdict;
%  (5) a certified point is a certificate SEMANTICALLY: P, R, Q evaluated at
%      the returned decision values satisfy (15a), (15b) to the solver's
%      accuracy (coefficient residual), R >= 0 and -(P*A+A*P+2kP*T) >= 0 on
%      random test functions and on the lowest mode, by quadrature
%      (HEATND_APPLY), and P*T >= ep^2 T*T there;
%  (6) TEST_COPQUADVAR_FACES (the face-weighted positive variables; was
%      HEATND_TEST_POSW, retired with heatNd_posw);
%  (7) kappa = r + k: the LITERAL Cor. 35 (15b), P*A + A*P + 2kP*T with A =
%      rT + A0 composed by the classes, gives the same SDP as HEATND_LPI's
%      kappa build from A0, and as the kappa build of heatNd_pie(N,0) from
%      its OWN base (N = 1, 2);
%  (8) HEATND_BISECT's bracket logic on synthetic verdicts (an oracle, no
%      solver): an uncertain hole below the threshold does not cap the
%      search; an uncertain band never gives 'resolved'; the retry is
%      taken; the lambda1 cap keeps trials at or below lambda1 + tol and
%      stops with its exact label; a certified +1 above lambda1 is rejected
%      from lo and stops the run (with the cap), or is flagged by
%      gap_clipped (without); the upward doubling stops at kmax, with the
%      uncertain points above lo in B.unc; and the input guards (k* <= 0,
%      ep <= 0, k < 0) error;
%  (9) the Farkas gate: a certified -1 has absolute cert_viol <= 1e-8 and a
%      finite certified radius, and the same ray is 'not verified' (0) when
%      the gate is set below it;
%  (10) a dump file stores its row set and reference route, and
%      HEATND_SOLVE(FILE) repeats that route (rows 'keep') by default; a
%      given OPTS.keep overrides a reference route on all rows.
%
% Initial coding MMP, 09/27/2026
% MMP, 09/27/2026 (review fixes): checks (7)-(10).
% MMP, 09/27/2026 (final reviews): (7) compares with an independently built
%   r = 0 base (the shared base made it 0 by construction); (8e) asserts the
%   exact stop label; (8f), (8g) new; (9) on the absolute violation; (10)
%   OPTS.keep against a file's route.
% MMP, 09/27/2026 (library face codes): (6) runs test_copquadvar_faces, the
%   test of copquadvar's new face codes that the benchmark now uses, in
%   place of heatNd_test_posw (retired with the copy heatNd_posw). The RNG
%   state is restored after it, so checks after (6) draw from the rng(1)
%   stream directly (before, from wherever heatNd_test_posw left it).
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

info = heatNd_path();
fprintf('test_heatNd_lpi: maxNumCompThreads = %d\n',info.threads);
rng(1);
R = struct('check',{},'val',{},'ok',{});
% heatNd_test_posw();                                                       % MMP, 09/27/2026 (was)
% R(end+1) = struct('check','(6) heatNd_test_posw','val',0,'ok',true);      % MMP, 09/27/2026 (was)
% RNG saved/restored: later draws no longer depend on (6)'s draw count.     % MMP, 09/27/2026
rs = rng;   test_copquadvar_faces();    rng(rs);                            % MMP, 09/27/2026
R(end+1) = struct('check','(6) test_copquadvar_faces','val',0,'ok',true);   % MMP, 09/27/2026
ep = 0.1;
% % % (1) counts.
for N = 1:2
    pie = heatNd_pie(N,0);
    for d = 0:2
        mu = d^2+4*d+3;
        [~,m1] = heatNd_lpi(pie,d,[],ep,struct('Pbasis','paper','RQ','Tspan','dp',0));
        [~,m2] = heatNd_lpi(pie,d,[],ep,struct('Pbasis','tensor','RQ','Tspan','dp',0));
        [~,m3] = heatNd_lpi(pie,d,[],ep,struct('Pbasis','paper','Pmult',false,'RQ','Tspan','dp',0));
        R(end+1) = chk(sprintf('(1) N=%d d=%d nP paper = mu(d)^N = %d',N,d,mu^N),m1.nP,m1.nP==mu^N); %#ok<AGROW>
        R(end+1) = chk(sprintf('(1) N=%d d=%d nP tensor',N,d),m2.nP,m2.nP==((d+1)+2*(d+1)^2)^N);   %#ok<AGROW>
        R(end+1) = chk(sprintf('(1) N=%d d=%d nP paper, no multiplier',N,d),m3.nP,m3.nP==((d+1)*(d+2))^N); %#ok<AGROW>
    end
    for dp = 0:2
        [~,m] = heatNd_lpi(pie,0,[],ep,struct('RQ','paper','dp',dp,'prune',false,'psatz','none'));
        mu = dp^2+4*dp+3;
        R(end+1) = chk(sprintf('(1) N=%d dp=%d R,Q paper Gram = mu^N = %d',N,dp,mu^N),m.sdp.Ks, ...
                       isequal(m.sdp.Ks,[mu mu].^N));                           %#ok<AGROW>
    end
    [~,m] = heatNd_lpi(pie,0,[],ep,struct('RQ','Tspan','dp',0,'prune',false,'psatz','none'));
    nint = 4*strcmp(pie.bc,'DD') + 3*~strcmp(pie.bc,'DD');
    g = prod(2+2*nint);
    R(end+1) = chk(sprintf('(1) N=%d Tspan unpruned Gram %d',N,g),m.sdp.Ks,isequal(m.sdp.Ks,[g g])); %#ok<AGROW>
    [~,m] = heatNd_lpi(pie,0,[],ep,struct('RQ','Tspan','dp',0,'prune',true,'psatz','none'));
    gR = prod(2*nint);  gQ = gR;
    for i = 1:N,    gQ = gQ + 2*prod(2*nint([1:i-1,i+1:N]));   end             % 2 mult monomials
    R(end+1) = chk(sprintf('(1) N=%d Tspan pruned Gram R %d Q %d',N,gR,gQ),m.sdp.Ks,isequal(m.sdp.Ks,[gR gQ])); %#ok<AGROW>
end
% % % (2)-(5) per N at the benchmark settings.
ob = struct('RQ','Tspan','dp',0,'psatz','linear');
for N = 1:2
    pie = heatNd_pie(N,0);  ks = pie.exact.kstar;
    of = ob;    of.eq15a = 'full';
    [~,mb] = heatNd_lpi(pie,0,[],ep,ob);    [~,mf] = heatNd_lpi(pie,0,[],ep,of);
    ob2 = ob;   ob2.base = mb.base;     of2 = of;   of2.base = mf.base;
    kk = [0.5 0.95 1.05]*ks;
    for k = kk
        ps = heatNd_lpi(pie,0,k,ep,ob2);    pf = heatNd_lpi(pie,0,k,ep,of2);
        Ds = heatNd_sdp(ps);    Df = heatNd_sdp(pf);
        rs = rk(Ds.At);     rf = rk(Df.At);     rj = rk([Ds.At, Df.At]);
        R(end+1) = chk(sprintf('(2) N=%d k=%.3fk* rank split %d = full %d = joint %d (m %d vs %d)', ...
                       N,k/ks,rs,rf,rj,Ds.m,Df.m),[rs rf rj],rs==rf && rf==rj);   %#ok<AGROW>
        vs = heatNd_solve(Ds);  vf = heatNd_solve(Df);
        R(end+1) = chk(sprintf('(2) N=%d k=%.3fk* verdict split %+d = full %+d',N,k/ks,vs.st,vf.st), ...
                       [vs.st vf.st],vs.st==vf.st);                             %#ok<AGROW>
        vl = heatNd_solve(ps,struct('route','lpisolve'));
        R(end+1) = chk(sprintf('(4) N=%d k=%.3fk* dump %+d = lpisolve %+d (rel_b %.1e / %.1e)', ...
                       N,k/ks,vs.st,vl.st,vs.rel_b,vl.rel_b),[vs.st vl.st],vs.st==vl.st); %#ok<AGROW>
    end
    % (3) affine in k.
    k1 = 0.5*ks;    k2 = ks;    k3 = 0.8*ks;
    D1 = heatNd_sdp(heatNd_lpi(pie,0,k1,ep,ob2));   D2 = heatNd_sdp(heatNd_lpi(pie,0,k2,ep,ob2));
    D3 = heatNd_sdp(heatNd_lpi(pie,0,k3,ep,ob2));
    Ai = D1.At + (k3-k1)*(D2.At-D1.At)/(k2-k1);
    e = full(max(abs(Ai(:)-D3.At(:))))/full(max(abs(D3.At(:))));
    R(end+1) = chk(sprintf('(3) N=%d At(k) affine: rel err %.1e, b equal',N,e),e, ...
                   e<1e-13 && isequal(D1.b,D3.b) && isequal(D1.K,D3.K));        %#ok<AGROW>
    % (5) semantic certificate at k = 0.5 k*.
    k = 0.5*ks;     p = heatNd_lpi(pie,0,k,ep,ob2);     D = heatNd_sdp(p);
    v = heatNd_solve(D,struct('keepx',true));
    R(end+1) = chk(sprintf('(5) N=%d k=0.5k* certified',N),v.st,v.st==1);      %#ok<AGROW>
    dv = full(D.RR*v.x)*D.bscl;                         % decvartable order
    nm = cellstr(string(p.decvartable));    nm = nm(1:numel(dv));
    B = mb.base;
    E15a = B.E1 - B.R;      E15b = B.X1 + k*B.X2 + B.Q;
    [a1,s1] = coefres(E15a,nm,dv);  [a2,s2] = coefres(E15b,nm,dv);
    R(end+1) = chk(sprintf('(5) N=%d coef residual (15a) %.1e (15b) %.1e',N,a1/s1,a2/s2), ...
                   [a1/s1 a2/s2],a1/s1<1e-5 && a2/s2<1e-5);                     %#ok<AGROW>
    Lop = subsd(B.X1 + k*B.X2,nm,dv);   Rop = subsd(B.R,nm,dv);
    Pop = subsd(B.E1,nm,dv);            % P*T - ep^2 T*T
    [Xq,wq] = qgrid(pie.dom,10);
    worst = -inf;   worstR = inf;   worstP = inf;   sc = 0;
    for t = 1:6
        if t==1,    f = pie.exact.eigfun;
        else,       Pw = randi([0 3],4,N);  cf = randn(4,1);
                    f = @(X) sum(cf'.*prod(reshape(X,[],1,N).^reshape(Pw,1,[],N),3),2);
        end
        fx = f(Xq);
        qL = sum(wq.*fx.*heatNd_apply(Lop,f,Xq,10));   % <v,(P*A+A*P+2kP*T)v>
        qR = sum(wq.*fx.*heatNd_apply(Rop,f,Xq,10));
        qP = sum(wq.*fx.*heatNd_apply(Pop,f,Xq,10));
        nT = sum(wq.*heatNd_apply(pie.T,f,Xq,10).^2);  % ||T v||^2 as scale
        worst = max(worst,qL/nT);   worstR = min(worstR,qR/nT);  worstP = min(worstP,qP/nT);
        sc = max(sc,abs(qL/nT));
    end
    R(end+1) = chk(sprintf('(5) N=%d <v,(P*A+A*P+2kP*T)v>/||Tv||^2 max %.2e <= 0',N,worst), ...
                   worst,worst<=1e-6*max(sc,1));                                 %#ok<AGROW>
    R(end+1) = chk(sprintf('(5) N=%d <v,Rv>/||Tv||^2 min %.2e >= 0, <v,(P*T-ep^2T*T)v> min %.2e',N,worstR,worstP), ...
                   [worstR worstP],worstR>=-1e-8 && worstP>=-1e-8);              %#ok<AGROW>
end
% % % (7) kappa: the literal Cor. 35 (15b) with A = r T + A0 at rate k is
% the SDP heatNd_lpi builds at kappa = r + k from A0, and the SDP of the
% r = 0 instance at kappa built from ITS OWN base (a shared base would give
% 0 by construction: both calls then evaluate base.X1 + kappa X2 + Q).
for N = 1:2
    p0 = heatNd_pie(N,0);   [~,m0] = heatNd_lpi(p0,0,[],ep,ob);     o0 = ob;    o0.base = m0.base;
    for r = [3.7 12]
        pr = heatNd_pie(N,r);   k = 0.3;
        [~,mb] = heatNd_lpi(pr,0,[],ep,ob);     B = mb.base;
        PA = B.P'*pr.A;                                 % Cor. 35 as printed
        Dl = heatNd_sdp(lpi_eq_cdopvar(B.prog,(PA+PA') + k*B.X2 + B.Q,'symmetric'));
        o2 = ob;    o2.base = B;
        Dk = heatNd_sdp(heatNd_lpi(pr,0,k,ep,o2));      % kappa build, r = r
        Dz = heatNd_sdp(heatNd_lpi(p0,0,r+k,ep,o0));    % kappa build, r = 0, own base
        same = Dl.m==Dk.m && isequal(Dl.K,Dk.K) && Dl.m==Dz.m;
        e1 = NaN;   e2 = NaN;
        if same
            e1 = full(max(abs(Dl.At(:)-Dk.At(:))))/full(max(abs(Dk.At(:))));
            e2 = full(max(abs(Dz.At(:)-Dk.At(:))))/full(max(abs(Dk.At(:))));
            same = same && isequal(Dl.b,Dk.b) && isequal(Dz.b,Dk.b);
        end
        R(end+1) = chk(sprintf('(7) N=%d r=%g: literal (15b) vs kappa build %.1e, vs r=0 own base %.1e',N,r,e1,e2), ...
                       [e1 e2],same && e1<=1e-14 && e2<=1e-14);            %#ok<AGROW>
    end
end
% % % (8) bracket logic on synthetic verdicts (oracle), and input guards.
p1 = heatNd_pie(1,0);
ora = @(f) @(kap,so) struct('st',f(kap,so),'why','oracle','rows',getf(so,'rows','all'),'tight',getf(so,'tight',false));
% (a) an uncertain hole below the threshold must not cap the search
% (kappa0 = [1 7] puts the third trial, 5.5, in the hole).
fa = @(kap,so) (kap<=6.2)*(1 - (kap>=5.2 && kap<5.8)) - (kap>6.2);
B = heatNd_bisect(p1,0,ep,struct(),struct('oracle',ora(fa),'verbose',false,'kappa0',[1 7]));
ia = find(B.trace(:,2)==0,1);
R(end+1) = chk(sprintf('(8a) hole [5.2,5.8): lo %.4f hi %.4f resolved %d, unc %s',B.lo,B.hi,B.resolved,mat2str(B.unc_all,4)), ...
               [B.lo B.hi],B.lo<=6.2 && B.hi>6.2 && B.lo>5.8 && B.resolved && B.hi-B.lo<=B.tol && ...
               ~isempty(ia) && any(B.trace(:,12)==3));                   %#ok<AGROW>
% (b) an uncertain band closes nothing: never 'resolved'.
fb = @(kap,so) (kap<=5) - (kap>=6);
B = heatNd_bisect(p1,0,ep,struct(),struct('oracle',ora(fb),'verbose',false,'retry',{{}},'maxsolve',14));
R(end+1) = chk(sprintf('(8b) band (5,6): lo %.4f hi %.4f resolved %d (%s)',B.lo,B.hi,B.resolved,B.stop), ...
               [B.lo B.hi],B.lo<=5 && B.hi>=6 && ~B.resolved);           %#ok<AGROW>
% (c) the retry: uncertain unless tight -> every final verdict is attempt 2.
fc = @(kap,so) getf(so,'tight',false)*((kap<=6.2) - (kap>6.2));
B = heatNd_bisect(p1,0,ep,struct(),struct('oracle',ora(fc),'verbose',false));
fin = B.trace(B.trace(:,2)~=0,:);
R(end+1) = chk(sprintf('(8c) retry tight: lo %.4f hi %.4f resolved %d',B.lo,B.hi,B.resolved), ...
               [B.lo B.hi],B.resolved && all(fin(:,12)==2) && all(fin(:,11)==1)); %#ok<AGROW>
% (e) the lambda1 cap: uncertain on (lambda1 - 0.01, 12], -1 above 12; no
% trial point above lambda1 + tol except the two upward-expansion points.
lam = p1.exact.lambda1;
fe = @(kap,so) (kap<=lam-0.01) - (kap>12);
B = heatNd_bisect(p1,0,ep,struct(),struct('oracle',ora(fe),'verbose',false,'retry',{{}}));
tried = B.trace(B.trace(:,12)>0,1);     over = tried(tried>lam+B.tol);
R(end+1) = chk(sprintf('(8e) cap: lo %.4f hi %.4f, %d of %d trials above lambda1 + tol (%s)',B.lo,B.hi, ...
               numel(over),numel(tried),B.stop),numel(over), ...
               numel(over)<=2 && ~B.resolved && B.lo<=lam-0.01 && ...
               strcmp(B.stop,'uncertain points fill [lo, lambda1 + tol] (Thm. 34 cap)')); %#ok<AGROW>
% (f) a certified +1 above lambda1 (the gate accepting an infeasible SDP):
% +1 up to 1.05 lambda1. With the cap it never becomes lo and the run stops
% at it; without, lo > lambda1 is flagged and gap(2) clipped to 0.
ff = @(kap,so) (kap<=1.05*lam) - (kap>1.05*lam);
B = heatNd_bisect(p1,0,ep,struct(),struct('oracle',ora(ff),'verbose',false,'retry',{{}}));
R(end+1) = chk(sprintf('(8f) +1 above lambda1, cap: lo %.6f (lambda1 %.6f), above %s, last trial %.6f (%s)', ...
               B.lo,lam,mat2str(B.above_kstar,8),B.trace(end,1),B.stop),B.lo, ...
               strcmp(B.stop,'certified +1 above lambda1: numerical acceptance') && ~B.resolved && ...
               B.lo<=lam && ~isempty(B.above_kstar) && all(B.above_kstar>lam) && ...
               B.trace(end,1)>lam && B.trace(end,2)==1 && B.gap(2)>=0 && ~B.gap_clipped); %#ok<AGROW>
B = heatNd_bisect(p1,0,ep,struct(),struct('oracle',ora(ff),'verbose',false,'retry',{{}},'cap',false));
R(end+1) = chk(sprintf('(8f) +1 above lambda1, no cap: lo %.6f, gap %s, clipped %d (%s)',B.lo, ...
               mat2str(B.gap,4),B.gap_clipped,B.stop),B.lo, ...
               B.lo>lam && B.gap_clipped && B.gap(2)==0 && ~isempty(B.above_kstar)); %#ok<AGROW>
% (g) +1 up to 9, uncertain above: the doubling stops at kmax (4 lambda1),
% and B.unc (hi NaN read as Inf) holds the uncertain points above lo.
fg = @(kap,so) double(kap<=9);
B = heatNd_bisect(p1,0,ep,struct(),struct('oracle',ora(fg),'verbose',false,'retry',{{}}));
R(end+1) = chk(sprintf('(8g) no -1: max trial %.4f (kmax %.4f), unc %s (%s)',max(B.trace(:,1)),4*lam, ...
               mat2str(B.unc,5),B.stop),max(B.trace(:,1)), ...
               strcmp(B.stop,'no certified -1 up to kmax') && max(B.trace(:,1))<=4*lam*(1+eps) && ...
               isnan(B.hi) && numel(B.unc)==3 && isequal(B.unc,B.unc_all) && all(B.unc>B.lo)); %#ok<AGROW>
% (d) guards: k* <= 0, ep <= 0, k < 0.
R(end+1) = chk('(8d) bisect errors for k* <= 0',0,throws(@() heatNd_bisect(heatNd_pie(1,12),0,ep,struct(), ...
               struct('oracle',ora(fa),'verbose',false)),'heatNd_bisect:kstar'));          %#ok<AGROW>
R(end+1) = chk('(8d) heatNd_lpi errors for ep = 0',0,throws(@() heatNd_lpi(p1,0,1,0,ob),'heatNd_lpi:ep')); %#ok<AGROW>
R(end+1) = chk('(8d) heatNd_lpi errors for k < 0',0,throws(@() heatNd_lpi(p1,0,-1,ep,ob),'heatNd_lpi:k'));  %#ok<AGROW>
% % % (9) the Farkas gate: a certified -1, and the same ray with ctol < 0.
[~,mb] = heatNd_lpi(p1,0,[],ep,ob);     o2 = ob;    o2.base = mb.base;
D = heatNd_sdp(heatNd_lpi(p1,0,1.1*p1.exact.lambda1,ep,o2));
v1 = heatNd_solve(D);   v2 = heatNd_solve(D,struct('ctol',-1));
R(end+1) = chk(sprintf('(9) 1.1 kappa*: st %+d cert_viol %.1e radius %.1e; ctol<0: st %+d',v1.st, ...
               v1.cert_viol,v1.cert_radius,v2.st),[v1.st v2.st],v1.st==-1 && v1.cert_viol<=1e-8 && ...
               v1.cert_radius>=1e8 && v2.st==0 && contains(v2.why,'not verified')); %#ok<AGROW>
% % % (10) a dump file repeats its stored route (rows, tolerances).
p2 = heatNd_pie(2,0);   [~,mb] = heatNd_lpi(p2,0,[],ep,ob);     o2 = ob;    o2.base = mb.base;
kap = 0.5*p2.exact.lambda1;     [pk,mk] = heatNd_lpi(p2,0,kap,ep,o2);   D = heatNd_sdp(pk);
keep = heatNd_lindep(D.At);
v1 = heatNd_solve(D,struct('keep',keep));
D.keep = keep;  D.ref = struct('st',v1.st,'rows','keep','tight',false,'rel_b',v1.rel_b);
f = [tempname '.mat'];  heatNd_sdp(D,f,mk);
v2 = heatNd_solve(f);   v3 = heatNd_solve(f,struct('rows','all'));
S = load(f,'info');     delete(f);
R(end+1) = chk(sprintf('(10) dump: ref %+d (%d rows); file %+d rows %s (%d); all rows %+d; kappa %.4f', ...
               v1.st,numel(keep),v2.st,v2.rows,v2.nrows,v3.st,S.info.kappa),[v1.st v2.st v3.st], ...
               v1.st==1 && v2.st==v1.st && strcmp(v2.rows,'keep') && v2.nrows==numel(keep) && ...
               strcmp(v3.rows,'all') && v3.nrows==D.m && S.info.kappa==kap && numel(keep)<D.m); %#ok<AGROW>
% A reference route on all rows is the file default, but a given OPTS.keep
% wins (HEATND_SOLVE header).
D.ref.rows = 'all';     f = [tempname '.mat'];  heatNd_sdp(D,f,mk);
v4 = heatNd_solve(f);   v5 = heatNd_solve(f,struct('keep',keep));   delete(f);
R(end+1) = chk(sprintf('(10) ref route all rows: file default rows %s (%d); OPTS.keep: rows %s (%d)', ...
               v4.rows,v4.nrows,v5.rows,v5.nrows),[v4.nrows v5.nrows], ...
               strcmp(v4.rows,'all') && v4.nrows==D.m && strcmp(v5.rows,'keep') && v5.nrows==numel(keep)); %#ok<AGROW>
ok = true;
fprintf('\n');
for i = 1:numel(R)
    ok = ok && R(i).ok;
    fprintf('%-86s %s\n',R(i).check,char(string(R(i).ok)));
end
assert(ok,'test_heatNd_lpi: a check failed.');
fprintf('test_heatNd_lpi: all %d checks pass.\n',numel(R));
end


function s = chk(c,v,ok),   s = struct('check',c,'val',v,'ok',logical(ok));   end

function x = getf(s,f,d)
% Field F of struct S, or D if absent.
if isfield(s,f),    x = s.(f);  else,   x = d;  end
end

function tf = throws(fh,id)
% True if FH() errors with identifier ID.
tf = false;
try,    fh();   catch ME,   tf = strcmp(ME.identifier,id);  end
end

function r = rk(At)
s = svd(full(At'*At));  r = nnz(s > max(s)*1e-12);
end

function [amax,scale] = coefres(X,names,vals)
% Largest kernel coefficient of the 1x1 'cdopvar' X at the decision values,
% and the no-cancellation scale max(|A| + |B||d|) (lpi_eq_sdopvar NOTES:
% by the canonical multiplier form X == 0 iff every coefficient is 0).
B = X.C{1};     zd = cellstr(string(B.Zd));
[tf,loc] = ismember(zd(:),names(:));    d = zeros(numel(zd),1);    d(tf) = vals(loc(tf));
amax = 0;   scale = 0;
for g = 1:numel(B.params.A)
    a = full(B.params.A{g}(:));     Bg = B.params.B{g};
    if isempty(Bg),     bd = zeros(size(a));    ad = bd;
    else,               bd = full(Bg'*d);       ad = full(abs(Bg)'*abs(d));
    end
    if isempty(a),      a = zeros(size(bd));    end
    if isempty(a),      continue,   end
    amax = max(amax,max(abs(a+bd)));    scale = max([scale;abs(a)+ad]);
end
end

function Y = subsd(X,names,vals)
% Fixed operator X(d) of a 1x1 'cdopvar' (params = unvec(A + B'd), sdopvar.m).
B = X.C{1};     zd = cellstr(string(B.Zd));
[tf,loc] = ismember(zd(:),names(:));    d = zeros(numel(zd),1);    d(tf) = vals(loc(tf));
nL = prod(cellfun(@numel,B.ZL));    nR = prod(cellfun(@numel,B.ZR));
m = B.dims(1)*nL;   n = B.dims(2)*nR;
prm = cell(size(B.params.A));
for g = 1:numel(prm)
    a = B.params.A{g};  Bg = B.params.B{g};
    if isempty(a),  a = sparse(m*n,1);  end
    if ~isempty(Bg),    a = a + Bg'*d;  end
    prm{g} = reshape(full(a),m,n);
end
Y = copvar({sopvar(prm,B.vars,B.ZL,B.ZR,B.dom,B.dims)});
end

function [X,w] = qgrid(dom,n)
b = (1:n-1)./sqrt(4*(1:n-1).^2-1);
[V,L] = eig(diag(b,1)+diag(b,-1));
[x,ix] = sort(diag(L));     wg = 2*V(1,ix)'.^2;     x = (x+1)/2;    wg = wg/2;
N = size(dom,1);    X = zeros(1,0);     w = 1;
for d = 1:N
    L = dom(d,2)-dom(d,1);  p = dom(d,1)+L*x;   q = L*wg;   m = size(X,1);
    X = [repmat(X,n,1), repelem(p,m,1)];    w = repmat(w,n,1).*repelem(q,m,1);
end
end
