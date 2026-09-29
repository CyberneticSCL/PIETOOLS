function B = heatNd_bisect(pie,d,ep,opts,bopts)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% B = HEATND_BISECT(PIE,D,EP,OPTS,BOPTS) largest rate certified by the
% Cor. 35 LPI (HEATND_LPI with PIE, D, EP, OPTS), by bisection over
% CERTIFIED verdicts only (HEATND_SOLVE: +1 feasible, -1 infeasible, 0
% uncertain - reported, never coerced).
%
% KAPPA. The SDP depends on kappa = r + k only (A = A0 + r T, HEATND_PIE),
% so the search is on kappa, over k >= 0 (Thm. 34), i.e. kappa >= r, and B
% reports kappa-hat, k-hat = kappa-hat - r and the r-invariant gap
% k* - k-hat = lambda1 - kappa-hat. k* <= 0 (r >= lambda1) is an ERROR:
% there is no k >= 0 to certify; bisect heatNd_pie(N,0) instead (the same
% SDPs, with k = kappa).
%
% AFFINE. (15a) does not involve kappa, and (15b) is X1 + kappa X2 + Q = 0
% on a fixed set of coefficients with b = 0 there, so
%   At(kappa) = At(kap1) + (kappa - kap1) dA,  dA = (At(kap2) - At(kap1))/(kap2 - kap1),
% with b and K independent of kappa. The program is built ONCE, on the
% r = 0 instance (base plus two finalizations, kap1 = lambda1/2, kap2 =
% lambda1), checked to share m, b and K, and every step is one sparse
% combination and one solve. TEST_HEATND_LPI (3) checks At(kappa) against
% a direct build.
%
% BRACKET. lo = largest kappa with a certified +1, hi = smallest with a
% certified -1. An uncertain point NEVER bounds the search: the next trial
% is the midpoint of the widest gap between consecutive points of {lo,
% uncertain points in (lo,hi), hi}, so the search continues on both sides
% of it. RESOLVED only when a certified -1 closes the bracket: hi - lo <=
% tol = max(rtol*|lo|, atol). Also stops at maxsolve (attempts), at the
% wall budget tmax (priced from the last solve), or when every gap is <=
% tol but hi - lo is not (uncertain points fill the bracket); B.stop says
% which. Start: kappa0 = [klo khi] (default [0.5 1.1] lambda1, klo clamped
% to >= r); an uncertified klo is halved toward kappa = r (k = 0, final);
% without a certified -1, kappa is doubled from khi until one appears, at
% most up to BOPTS.kmax (stop 'no certified -1 up to kmax').
% CAP (BOPTS.cap, default true: k* is exact and Thm. 34 applies). Every
% kappa > lambda1 is infeasible (a certificate would give a rate above the
% exact k*), so trial points are placed at or below min(hi, lambda1 + tol),
% and a certified +1 at kappa > lambda1 is a CONTRADICTION: the gate accepts
% infeasible SDPs within ~1e-8 of the threshold (MEASURED, review
% 09/27/2026: 2-D heavy d = 0, +1 at lambda1 + 1e-8 with rel_b 4.1e-7). It
% is kept out of lo, recorded in B.above_kstar, and every loop (klist
% included) stops: 'certified +1 above lambda1: numerical acceptance',
% resolved false. lo = lambda1 exactly stops with 'lo = lambda1 (Thm. 34
% bound)'. In the recorded runs (2-D listingL d = 1, tight tolerances)
% MOSEK returned UNKNOWN above lambda1 on both row sets, so solves there
% mostly spend the budget. The verdicts stay certified; the theorem only
% decides where to solve and which +1 to reject.
% Verdicts are assumed monotone in kappa (the operator LPI is: (P,R,Q) for
% kappa gives Q + 2(kappa-kappa')P*T for kappa' < kappa; a finite Gram
% basis need not be); B.monotone is false if a certified -1 lies at or
% below lo.
%
% RETRY. A 0 is re-solved per BOPTS.retry, in order, stopping at the first
% certified verdict; the changes accumulate: 'tight' (MOSEK tolerances
% 1e-10, HEATND_SOLVE opts.tight), 'loose' (MOSEK defaults again: tight
% tolerances can end UNKNOWN where the defaults find a Farkas ray), 'rows'
% (the other row set: all rows <-> the HEATND_LINDEP set, computed once if
% needed). Every attempt is a trace row; the bracket uses the last
% attempt's verdict.
%
% BOPTS (defaults): kappa0 [0.5 1.1]*lambda1 (absolute kappa), rtol 1e-3,
%   atol 1e-6, maxsolve 30 (attempts, retries included), tmax Inf (s),
%   trace0 (rows [kappa st] of earlier final verdicts, reused without
%   solving), kaff (a B.kaff from an earlier call or file: skips the build),
%   kaff_file (load if it exists, else save the affine SDP there, -v7.3;
%   files written before the kappa change carry k1 at meta.r: kap1 = k1 +
%   meta.r), solve (opts for HEATND_SOLVE, e.g. tight), retry ({'tight',
%   'rows'}), cap (true), kmax (4 lambda1: bound of the upward doubling),
%   verbose (true), lindep (false: true solves on a maximal
%   independent row set from HEATND_LINDEP at kap1, computed once and stored
%   in kaff; verdicts are still on all rows), klist ([]: absolute kappas to
%   solve in order instead of searching - for costly 3-D points), oracle
%   ([]: v = oracle(kappa, solve_opts) replaces the SDP solve, v with at
%   least st, why; for testing the bracket logic, TEST_HEATND_LPI (8), and
%   as HEATND_POINCARE's solve hook: there kappa is the Poincare lam, PIE
%   is a stub with r = 0 and exact.lambda1 = exact.kstar = the exact
%   constant, and D, EP, OPTS are unused).
%
% OUTPUT struct B: lo, hi (kappa), kappa_hat (= lo), khat (= lo - r), gap
%   ([max(0, lambda1 - hi), max(0, lambda1 - lo)]: k* - k-hat lies in
%   [gap(1), gap(2)]; clipped at 0 since the LPI is sufficient, kappa-hat <=
%   lambda1; NaN if lo is), gap_clipped (true if lambda1 - lo < 0 was
%   clipped: only with cap false), above_kstar (every certified +1 kappa >
%   lambda1), unc (uncertain kappas in (lo,hi), a NaN end read as -Inf /
%   Inf), unc_all (every uncertain kappa),
%   resolved, stop, monotone, r, lambda1, kstar, tol, trace (one row per
%   attempt: [kappa st rel_b psd_relmin t_mosek iter eta mosek_threads
%   cert_viol rows tight attempt]; cert_viol = HEATND_SOLVE's gated
%   absolute Farkas violation; rows 1 = independent set, 0 = all;
%   attempt 0 = reused from trace0), why (cell of verdict texts), nsolve,
%   t_build, t_total, sdp (m, nx, Ks, nnz), meta (HEATND_LPI meta without
%   base), kaff (for reuse).
%
% Initial coding MMP, 09/27/2026
% MMP, 09/27/2026 (review fixes): search on kappa = r + k; an uncertain
%   point no longer caps the search (it capped it before, and 'resolved'
%   could be set by one); retry policy; tolerance relative to |lo|; error
%   for k* <= 0; MOSEK's thread count, row set and tolerances in the trace.
% MMP, 09/27/2026 (poincare): header only - the oracle is also the solve
%   hook of HEATND_POINCARE (b-affine family); no functional change.
% MMP, 09/27/2026 (final reviews): with the cap, a certified +1 above
%   lambda1 is rejected and stops the run (the gate accepts infeasible SDPs
%   just above the threshold; before, it became lo under a benign label and
%   gap(2) went negative); the unreachable 'lo within tol of lambda1' test
%   is replaced by lo >= lambda1; gap(2) clipped and flagged; B.unc with a
%   NaN end; upward doubling bounded by kmax; trace column 9 is cert_viol
%   (HEATND_SOLVE now gates on it; it was cert_rel).
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if nargin<4 || isempty(opts),   opts = struct();    end
if nargin<5 || isempty(bopts),  bopts = struct();   end
if isfield(bopts,'k0')
    error('heatNd_bisect:k0','BOPTS.k0 (rates) was replaced by BOPTS.kappa0 (absolute kappa = r + k).')
end
r = pie.r;  lam1 = pie.exact.lambda1;   ks = pie.exact.kstar;
if ks<=0
    error('heatNd_bisect:kstar',['k* = %g <= 0 (r >= lambda1): no rate k >= 0 to certify (Thm. 34). ' ...
          'The SDP depends on kappa = r + k only: bisect heatNd_pie(N,0) instead.'],ks)
end
df = struct('kappa0',[0.5 1.1]*lam1,'rtol',1e-3,'atol',1e-6,'maxsolve',30,'tmax',Inf, ...
            'trace0',zeros(0,2),'kaff',[],'kaff_file','','solve',struct(),'retry',{{'tight','rows'}}, ...
            'cap',true,'kmax',4*lam1,'verbose',true,'lindep',false,'klist',[],'oracle',[]);
fn = fieldnames(df);
for i = 1:numel(fn),    if ~isfield(bopts,fn{i}),   bopts.(fn{i}) = df.(fn{i});     end,    end
T0 = tic;

% % % The k-affine SDP (not needed with an oracle).
K = bopts.kaff;
if isempty(bopts.oracle)
    if isempty(K) && ~isempty(bopts.kaff_file) && exist(bopts.kaff_file,'file')
        S = load(bopts.kaff_file,'kaff');   K = S.kaff;
    end
    if isempty(K)
        K = build_kaff(pie,d,ep,opts,lam1);
        if ~isempty(bopts.kaff_file),   kaff = K;   save(bopts.kaff_file,'kaff','-v7.3');   end
    end
    if ~isfield(K,'kap1')                       % pre-kappa file: built at meta.r
        K.kap1 = K.k1 + K.meta.r;
    end
    if bopts.lindep && ~isfield(K,'keep'),  K = add_keep(K,bopts);  end
end
tb = toc(T0);

% % % State: final verdicts (kv, sv) and the attempt trace.
tr = zeros(0,12);   why = {};   kv = zeros(0,1);    sv = zeros(0,1);
for i = 1:size(bopts.trace0,1)                      % reuse earlier verdicts
    kv(end+1,1) = bopts.trace0(i,1);    sv(end+1,1) = bopts.trace0(i,2);            %#ok<AGROW>
    tr(end+1,:) = [bopts.trace0(i,1:2) NaN(1,9) 0];  why{end+1} = 'reused';         %#ok<AGROW>
end
% A certified +1 above lambda1 contradicts Thm. 34 (header, CAP): with the
% cap it never enters lo, and ANOM stops every loop through budget().
anom = bopts.cap && any(sv==1 & kv>lam1);
ns = 0;     tlast = 0;
    function [lo,hi,unc] = bracket()
        lo = max([kv(sv==1 & ~(bopts.cap & kv>lam1)); -Inf]);  hi = min([kv(sv==-1); Inf]);
        if isinf(lo),   lo = NaN;   end
        if isinf(hi),   hi = NaN;   end
        unc = kv(sv==0);
    end
    function t = tolw(lo)
        t = max(bopts.rtol*abs(lo),bopts.atol);
    end
    function ok = budget()
        % False also after a contradiction (ANOM): no further solve is trusted.
        ok = ~anom && ns<bopts.maxsolve && toc(T0)+1.1*tlast<bopts.tmax;
    end
    function st = attempt(kap,so,a)
        % One solve at kappa with HEATND_SOLVE options SO; a trace row.
        if isempty(bopts.oracle)
            D = K.D1;   D.At = K.D1.At + (kap-K.kap1)*K.dA;
            v = heatNd_solve(D,so);
        else
            v = vfill(bopts.oracle(kap,so));
        end
        st = v.st;  ns = ns+1;  tlast = v.wall;
        tr(end+1,:) = [kap st v.rel_b v.psd_relmin v.t_mosek v.iter v.eta v.mosek_threads ...
                       v.cert_viol strcmp(v.rows,'keep') v.tight a];
        why{end+1} = v.why;
        if bopts.verbose
            fprintf(['  kappa = %-12.8g (k %-10.6g, %.4f lambda1)  st = %+d  rel_b %.1e  eta %.1e  psd %+.1e  ' ...
                     'mosek %.1f s, %d it, %g thr, rows %s, tight %d, try %d  %s\n'], ...
                    kap,kap-r,kap/lam1,st,v.rel_b,v.eta,v.psd_relmin,v.t_mosek,v.iter,v.mosek_threads, ...
                    v.rows,v.tight,a,v.why);
        end
    end
    function st = solve_at(kap)
        % Certified verdict at kappa, with the retries of BOPTS.retry.
        so = bopts.solve;
        if bopts.lindep,    so.keep = K.keep;   so.rows = 'keep';   end
        a = 1;  st = attempt(kap,so,a);
        for jr = 1:numel(bopts.retry)                % jr: not shared with the parent
            if st~=0 || ~budget(),  break,  end
            switch bopts.retry{jr}
                case 'tight'
                    if isfield(so,'tight') && so.tight,     continue,   end
                    so.tight = true;
                case 'loose'
                    if ~isfield(so,'tight') || ~so.tight,   continue,   end
                    so.tight = false;
                case 'rows'
                    if isfield(so,'rows') && strcmp(so.rows,'keep')
                        so.keep = [];   so.rows = 'all';
                    else
                        if isempty(bopts.oracle)
                            if ~isfield(K,'keep'),  K = add_keep(K,bopts);  end
                            so.keep = K.keep;
                        end
                        so.rows = 'keep';
                    end
                otherwise
                    error('heatNd_bisect:retry','Retry ''%s'' is not ''tight'', ''loose'' or ''rows''.',bopts.retry{jr})
            end
            a = a+1;    st = attempt(kap,so,a);
        end
        kv(end+1,1) = kap;  sv(end+1,1) = st;
        if bopts.cap && st==1 && kap>lam1,  anom = true;    end     % header, CAP
    end
    function st = verdict_at(kap)
        % A known final verdict is reused (search mode), else solved.
        i0 = find(kv==kap,1,'last');
        if ~isempty(i0),    st = sv(i0);    else,   st = solve_at(kap);     end
    end

stop = 'budget';
if ~isempty(bopts.klist)                            % fixed list, no search
    for kk = reshape(bopts.klist,1,[])
        if ~budget(),   break,  end
        solve_at(kk);
    end
    stop = 'klist';
else
    % % % Lower end: a certified +1, halving toward k = 0.
    [lo,hi] = bracket();
    if isnan(lo) && budget()
        klo = max(bopts.kappa0(1),r);
        if ~isnan(hi) && klo>=hi,   klo = r + (hi-r)/2;     end
        st = verdict_at(klo);
        while st~=1 && klo>r && budget()
            klo = r + (klo-r)/2;    if klo-r<1e-3*lam1,     klo = r;    end
            st = verdict_at(klo);
        end
    end
    % % % Upper end: a certified -1 (an uncertain point does not stop this).
    [lo,hi] = bracket();
    % Doubling bounded by kmax: without it, UNKNOWN above lambda1 spends the
    % whole budget at absurd kappa (review: 1.1e4 on an oracle).
    if ~isnan(lo) && isnan(hi) && budget()
        kh = min(max(bopts.kappa0(2),1.05*lo+eps),bopts.kmax);
        st = verdict_at(kh);
        while st~=-1 && budget()
            if kh>=bopts.kmax,  stop = 'no certified -1 up to kmax';   break,  end
            kh = min(2*max(kv),bopts.kmax);     st = verdict_at(kh);
        end
    end
    % % % Bisect the widest gap of {lo, uncertain points in (lo,top), top},
    % top = hi, or min(hi, lambda1 + tol) with the cap.
    while budget()
        [lo,hi,unc] = bracket();
        if isnan(lo) || isnan(hi) || hi<=lo,    break,  end
        if hi-lo<=tolw(lo),     break,  end
        % With the cap lo <= lambda1 (bracket excludes +1 above it), so lo
        % >= lambda1 means lo = lambda1: nothing above it can be feasible.
        if bopts.cap && lo>=lam1,   stop = 'lo = lambda1 (Thm. 34 bound)';  break,  end
        top = hi;   if bopts.cap,   top = min(hi,lam1+tolw(lo));    end
        p = unique([lo; unc(unc>lo & unc<top); top]);   g = diff(p);
        [gm,j] = max(g);
        if gm<=tolw(lo)
            if top<hi,  stop = 'uncertain points fill [lo, lambda1 + tol] (Thm. 34 cap)';
            else,       stop = 'uncertain points fill the bracket';
            end
            break
        end
        solve_at((p(j)+p(j+1))/2);
    end
end
[lo,hi,unc] = bracket();
tl = tolw(lo);
res = ~anom && ~isnan(lo) && ~isnan(hi) && hi>lo && hi-lo<=tl;
if anom,                stop = 'certified +1 above lambda1: numerical acceptance';   % header, CAP
elseif res,             stop = 'resolved';
elseif strcmp(stop,'klist')
elseif isnan(lo),       stop = 'no certified feasible kappa';
elseif isnan(hi) && strcmp(stop,'budget'),  stop = 'no certified infeasible kappa';  % keeps 'kmax'
elseif hi<=lo,          stop = 'non-monotone: certified -1 at or below lo';
end
mono = ~any(kv(sv==-1) <= lo);
if isempty(K) || ~isfield(K,'D1')
    sdp = struct('m',NaN,'nx',NaN,'Ks',[],'Kf',NaN,'nnz',NaN);  mt = [];
else
    sdp = struct('m',K.D1.m,'nx',K.D1.nx,'Ks',K.D1.K.s,'Kf',K.D1.K.f,'nnz',K.D1.nnz);  mt = K.meta;
end
% gap clipped at 0: the LPI is sufficient, so kappa-hat <= lambda1 (a
% certified -1 may lie above lambda1; lo above it only with cap false,
% flagged by gap_clipped). g2 < 0 is false for NaN, so a NaN lo stays NaN.
g2 = lam1 - lo;     g2c = g2;   if g2<0,    g2c = 0;    end
loI = lo;   if isnan(loI),  loI = -Inf;     end     % unc: a NaN end is open
hiI = hi;   if isnan(hiI),  hiI = Inf;      end
B = struct('lo',lo,'hi',hi,'kappa_hat',lo,'khat',lo-r,'gap',[max(0,lam1-hi), g2c], ...
           'gap_clipped',g2<0,'above_kstar',kv(sv==1 & kv>lam1)', ...
           'unc',unc(unc>loI & unc<hiI)','unc_all',unc','resolved',res,'stop',stop, ...
           'monotone',mono,'r',r,'lambda1',lam1,'kstar',ks,'tol',tl, ...
           'trace',tr,'why',{why},'nsolve',ns,'t_build',tb,'t_total',toc(T0), ...
           'sdp',sdp,'meta',mt,'kaff',K);
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function K = build_kaff(pie,d,ep,opts,lam1)
% Base once, two finalizations, difference. Built on the r = 0 instance
% (the SDP depends on kappa only), so k = kappa exactly; kap1, kap2 ~= 0 so
% no coefficient of X2 is structurally absent from either.
pie0 = pie;     pie0.r = 0;     pie0.A = pie.A0;    pie0.exact.kstar = lam1;
kap1 = lam1/2;  kap2 = lam1;
t0 = tic;
[~,meta] = heatNd_lpi(pie0,d,[],ep,opts);
o2 = opts;  o2.base = meta.base;
p1 = heatNd_lpi(pie0,d,kap1,ep,o2);     D1 = heatNd_sdp(p1);    clear p1
p2 = heatNd_lpi(pie0,d,kap2,ep,o2);     D2 = heatNd_sdp(p2);    clear p2
if D1.m~=D2.m || ~isequal(D1.K,D2.K) || D1.nx~=D2.nx
    error('heatNd_bisect:affine','The SDP shape depends on kappa; not affine.')
end
if norm(D1.b*D1.bscl - D2.b*D2.bscl,inf) > 1e-12*max(norm(D1.b*D1.bscl,inf),1)
    error('heatNd_bisect:affine','b depends on kappa; not affine.')
end
% D2's normalization equals D1's (same b), so At is on one scale.
dA = (D2.At - D1.At)/(kap2-kap1);
D1.RR = [];                                         % not needed for solving
mt = meta;  mt = rmfield(mt,'base');    mt.t_kaff = toc(t0);
K = struct('D1',D1,'dA',dA,'kap1',kap1,'kap2',kap2,'meta',mt);
end


function K = add_keep(K,bopts)
% Independent rows at kap1 (HEATND_LINDEP), once; saved with the kaff file.
[K.keep,K.lindep] = heatNd_lindep(K.D1.At);
if bopts.verbose,   fprintf('  lindep: rank %d of %d rows (%.1f s)\n',K.lindep.rank,K.lindep.m,K.lindep.t);  end
if ~isempty(bopts.kaff_file),   kaff = K;   save(bopts.kaff_file,'kaff','-v7.3');   end
end


function v = vfill(v)
% Oracle verdicts: fill the fields the trace records.
df = struct('st',0,'why','oracle','rel_b',NaN,'psd_relmin',NaN,'t_mosek',NaN,'iter',NaN, ...
            'eta',NaN,'mosek_threads',NaN,'cert_viol',NaN,'rows','all','tight',false,'wall',0);
fn = fieldnames(df);
for i = 1:numel(fn),    if ~isfield(v,fn{i}),   v.(fn{i}) = df.(fn{i});     end,    end
end
