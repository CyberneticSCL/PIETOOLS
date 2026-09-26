function R = bl_bisect(dumpfile,opts)                                        % CC, 09/26/2026
% BL_BISECT  Certified bisection on the objective variable of an objective-form
% PIETOOLS SDP, one feasibility probe per step, solved by cuADMM or Mosek.
%
%   R = bl_bisect(dumpfile,struct('hi',g))   % defaults: cuadmm_settings(1).cuadmm.bisect
%   R = bl_bisect(dumpfile,opts)             % any field of opts overrides a default
%
% dumpfile is a <case>.mat written by bl_run (SeDuMi At/b/c/K, plus
% <case>_meta.mat).  The objective must be one unit entry in the free block, and
% that variable is what is pinned: gam^2, not gam, for the coercive H2
% executives (PIETOOLS_H2_norm_o_coercive.m:172 returns sqrt of it).
%
% THE RULE (agreed with the maintainer 2026-09-26; evidence: BASELINE_REPORT.md)
%  1. Verdict per probe gamma:
%       F  verified feasible.  The solver's X is REPAIRED onto the pinned affine
%          set (alternating projections, min-norm steps by lsqr), and the
%          repaired point meets the rows to rep_tol, has every PSD block's
%          smallest eigenvalue >= -psd_tol x the largest eigenvalue over all
%          blocks, and satisfies the UNPINNED dump's rows with x_j = gamma after
%          undoing bl_fix's scaling.  A residual gate alone is not a certificate:
%          gate slack becomes bound error.  psd_tol > 0 is forced: measured on
%          hinf_rd1 at 2 gamma*, two blocks keep two eigenvalues below 1e-8 x max
%          however much slack there is (no interior in those blocks), and 52 of
%          243 singular values of the row system are below 1e-10 x max, so an
%          exact repair leaves them at -5e-9 (1.3e-8 relative); on cuADMM's X
%          the projections stall at -1.1e-8 for 300 rounds.  Exactness would
%          need facial reduction (sospsimplify) or an operator-level margin.
%          CALIBRATION (hinf_rd1, cuADMM 20k iterations, 8 rounds), repaired
%          lambda_min/lambda_max at gamma/gamma* = 0.99, 0.999, 0.9999 |
%          1.0001, 1.001, 1.01, 1.2:  -3.6e-5, -4.3e-6, -6.9e-7 | -2.3e-7,
%          -2.9e-9, -3.6e-10, -4.4e-9.  psd_tol = 1e-7 rejects every infeasible
%          point measured (0.9999 by 7x); by the near-linear trend it admits an
%          infeasible gamma only within ~2.5e-5 of gamma* (extrapolated), 40x
%          below rtol.  The constant is case-specific: calibrate per class.
%          CC, 09/26/2026: THE RELATIVE TEST ALONE IS UNSOUND, measured.  SeDuMi
%          (numerr 1) at 0.9982 and 0.99686 gamma* of hinf_rd1, where Mosek's
%          Farkas certificates verify exactly, returned points with lambda_max
%          2920 / 3785 and ||X|| ~2.5e3 (500x a genuine point's): relative
%          lambda_min -3.5e-8 / -4.7e-8 passed, absolute -1.0e-4 / -1.8e-4.  So F
%          also requires lambda_min >= -psd_abs (1e-6; genuine certificates so
%          far are at <= 1.5e-7 absolute) and ||X|| <= norm_ratio (50) x the
%          seed's certified ||X||.  Both guards can only withhold an F.
%       I  verified infeasible.  y with A'y = 0 on the free block (after an
%          exact correction), A'y NSD on every PSD block, b'y > 0: no feasible X.
%       f  leans feasible: the running-min pinf reached tau (item 4) within
%          the cap.  cuADMM runs with tol = tau, so it stops early when its dual
%          side has also converged.
%       i  leans infeasible: it did not.  Mosek UNKNOWN is also i.
%     Leans only steer.  The bound reported is gamma_F.
%  2. Two brackets: certified [gamma_I, gamma_F]; estimate [gamma_i, gamma_f],
%     gamma_f = min over F,f and gamma_i = max over I,i.  opts.lo (e.g. a
%     closed-form gain) is an I point by theorem, never probed.
%  3. Steer on the estimate: the next gamma is the midpoint of [gamma_i,gamma_f]
%     (geometric once gamma_i > 0).
%  4. The lean threshold follows the bracket INTERVAL in residual units, not
%     the ratio.  Pinning moves the unit-normalised b along a chord, and cuADMM's
%     pinf divides by 1+||b|| = 2, so tau = c_lean*chord(gamma_i,gamma_f)/2.
%     ERROR BAND: the infeasible plateau is ~kappa*chord(gamma,gamma*)/2 with
%     kappa 0.26-0.32 (hinf_rd1/_hv; 0.87 hinfdu_rd1), so an infeasible probe
%     can lean f when it is within ~(c_lean/kappa) of the interval from
%     gamma*: ~1/6 of it at c_lean = 0.05.  Measured in the acceptance run
%     (hinf_rd1, 2026-09-26): 0.9982 gamma* read f in a 1.6% bracket, the
%     estimate slid below gamma*, and the endgame recovered (F at 1.0002 gamma*
%     after two failed certifications).  A wrong lean costs time, not validity.
%     A fixed 1e-3 threshold instead misreads 0.99 gamma* even in a 2%
%     bracket, an error that does not shrink with the bracket.
%  5. Cap: cap = kmin*(D0/D)^p, p = 0.7 (measured: iterations to a stable lean
%     ~ d^-0.7 on hinf_rd1/_hv), and never below 2x the iterations the
%     certified seed took to reach tau.  Wide bracket, few iterations.
%  6. Steering stops when the estimate is within rtol, when tau < tau_min, or
%     when even the seed needed more than kmax/2 iterations to reach tau (the
%     pin moves b by less than the solver resolves: est_rd1's failure mode).
%  7. Endgame: certify at gamma_f*(1+d), d = 0, delta, 2 delta, 4 delta, ...
%     (tol cert_tol, cert_iter iterations) until an F, or gamma_F is reached.
%     The first gamma_F is opts.hi, certified before any probe and doubled up
%     to hi_doublings times until it is F.  Endgame runs move only the
%     certified bracket: an unconverged certification run is not a lean.
%  Mosek probes decide F and I directly; the same verifiers check its output.
%
% NOT DONE: warm start (cuadmm_exe cannot, main.cu:66); a secant step on the
% infeasible plateaus (3 cases of evidence); steering by slope or Aitken (the
% slope's boundary value differs by executive).  Both are recorded per probe.
%
% PERFORMANCE.  Per case: one svec map.  Per probe only b.txt is rewritten (At
% does not depend on gamma); the checks are lsqr on As (O(nnz(At)) per
% iteration, no factorisation: As'As is rank-deficient on these programs and a
% shifted Cholesky amplified the residual 1e3-fold) plus one eig per PSD block.
%
% ACCEPTANCE (hinf_rd1, m = 243, gamma* = 0.182626, 2026-09-26):
%   mosek   certified [0.99932, 1.00000] gamma*, 12 verified certificates, 0.3 s
%   cuadmm  certified [0, 1.00020] gamma*, 15 probes (F2 f2 i11), 552 s solver,
%           245k iterations.  No exact Farkas certificate from cuADMM: its y
%           gives b'y > 0 but A'y leaves the NSD cone by ~1e-5 relative.
%
% OUTPUT R: gamma_F, gamma_I, gamma_f, gamma_i, the probe table R.probes (also
% <outdir>/<case>_bisect.tsv), R.flags (order violations, stops), R.t_solver,
% R.iters.
%
% CC, 09/26/2026 (harness for the unattended regime, bl_regime.m):
%   opts.mode  'bisect' (default, as above); 'probe': certification runs at the
%              gammas in opts.gammas, no steering; 'feas': one certification run
%              on a program with no objective (stability, well-posedness).
%   opts.solver 'sedumi' added (a second opinion where Mosek is UNKNOWN).
%   opts.pin_scale  s in bl_fix's pin row (s*e_j'x = s*gam); 'auto' = ||b_un||/hi.
%   opts.tag   suffix for the run directory and tsv names.
%   Kill switches, checked before EVERY solver call, seed loop included:
%   opts.stopfile exists, opts.max_wall (s since start), opts.deadline (posix s),
%   each against the call's estimated length; opts.run_timeout caps one cuADMM
%   call (a timed-out call has no CUADMM_TIMING line and reads 'i').
%   Each cuADMM run's X and y are kept as run_NNN_X.txt / run_NNN_y.txt, and
%   each probe is appended to <tag>_<solver>_probes.tsv as it finishes, so a
%   crash loses nothing already measured.

if nargin < 2, opts = struct(); end
cuadmm_path;
P = cuadmm_settings(1);  P = P.cuadmm.bisect;
fn = fieldnames(opts);
for i = 1:numel(fn), P.(fn{i}) = opts.(fn{i}); end
if isempty(P.delta), P.delta = P.rtol; end
if isempty(P.exe), P.exe = getenv('CUADMM_EXE'); end
if isempty(P.exe), P.exe = '/home/mpeet/solvers/cuADMM/build/cuadmm_exe'; end
if isempty(P.stopfile), P.stopfile = fullfile(cuadmm_outdir(),'STOP'); end
if strcmp(P.mode,'bisect') && ~(numel(P.hi) == 1 && P.hi > 0)
    error('bl_bisect:hi','opts.hi (a positive gamma expected to be feasible) is required');
end
if strcmp(P.mode,'probe') && isempty(P.gammas), error('bl_bisect:gammas','probe mode needs opts.gammas'); end
if isempty(P.probe_tol), P.probe_tol = P.cert_tol; end
if isempty(P.probe_cap), P.probe_cap = P.cert_iter; end
t0 = tic;

% ---- the case, pinned once: only b depends on gamma
[~,id] = fileparts(dumpfile);
tag = id;  if ~isempty(P.tag), tag = [id '_' P.tag]; end
if isempty(P.outdir), P.outdir = fullfile(cuadmm_outdir(),'bisect'); end
if ~exist(P.outdir,'dir'), mkdir(P.outdir); end
S  = load(dumpfile);
Mt = load(strrep(dumpfile,'.mat','_meta.mat'));
K  = S.K;  if ~isfield(K,'f') || isempty(K.f), K.f = 0; end
Ks = double(K.s(:)');
nvar = size(S.At,1);
C.b_un = full(S.b(:))*Mt.bscl;               % bl_fix's unnormalised b
C.B    = norm(C.b_un);
C.dir  = fullfile(P.outdir,[tag '_' P.solver]);
if strcmp(P.mode,'feas')                     % no objective, nothing to pin
    if nnz(S.c) > 0, error('bl_bisect:feas','feas mode needs c = 0; this dump has an objective'); end
    C.j = [];  C.s = 1;  C.At2 = S.At;
    dump2cuadmm(dumpfile,C.dir);             % the dump's b is already unit-norm
else
    g0 = P.hi;  if isempty(g0), g0 = max(P.gammas); end
    s  = P.pin_scale;  if ischar(s), s = C.B/g0; end   % 'auto'
    G0 = bl_fix(dumpfile,g0,C.dir,s);        % validates the objective; writes At once
    C.j = G0.obj_idx;  C.s = s;  C.At2 = G0.At;
end
C.At0  = S.At;  C.K = K;  C.Ks = Ks;  C.nvar = nvar;  C.P = P;  C.t0 = t0;
vec_len = K.f + sum(Ks.*(Ks+1)/2);
C.T    = svecmap(K.f,Ks,vec_len,nvar);
C.As   = C.T*C.At2;                          % vec_len x m: constraint rows on symmetric X
C.live = fullfile(P.outdir,sprintf('%s_%s_probes.tsv',tag,P.solver));

R.id = tag;  R.solver = P.solver;  R.mode = P.mode;  R.flags = {};  R.t_solver = 0;
R.iters = 0;  R.spi = 5e-3;  R.stopped = '';
R.probes = repmat(newpr(),0,1);
fid = fopen(C.live,'w');  fprintf(fid,'%s\n',strjoin(fieldnames(newpr())',sprintf('\t')));  fclose(fid);
cert = struct('F',inf,'I',P.lo);             % gamma_F, gamma_I
est  = struct('f',inf,'i',P.lo);             % gamma_f, gamma_i

% ---- probe / feas modes: certification runs only, no steering
if ~strcmp(P.mode,'bisect')
    gl = P.gammas;  if strcmp(P.mode,'feas'), gl = NaN; end
    for q = 1:numel(gl)
        [R,stop] = mustStop(R,C,P.probe_cap);  if stop, break; end
        pr = probe(C,gl(q),P.probe_tol,P.probe_cap,P.mode,true,numel(R.probes)+1);
        R  = logprobe(R,pr,C);
        if ~isnan(gl(q)), [cert,est,R] = update(cert,est,gl(q),pr.verdict,R,P); end
    end
    R = finish(R,cert,est,C);  return
end

% ---- 7. seed: a certified upper end before any probe
g = P.hi;  seedtr = [];
for d = 0:P.hi_doublings
    [R,stop] = mustStop(R,C,P.cert_iter);  if stop, R = finish(R,cert,est,C); return; end
    [pr,tr] = probe(C,g,P.cert_tol,P.cert_iter,'seed',true,numel(R.probes)+1);
    R = logprobe(R,pr,C);
    [cert,est,R] = update(cert,est,g,pr.verdict,R,P);   % a seed I or lean counts too
    if strcmp(pr.verdict,'F'), seedtr = tr;  C.normx_ref = pr.F_normx;  break; end
    g = 2*g;
end
if ~isfinite(cert.F)
    R.flags{end+1} = sprintf('no certified upper bound up to %g',g/2);
    R = finish(R,cert,est,C);  return
end

% ---- 3-6. steer on the estimate bracket
D0 = chord(C.s*est.i,C.s*est.f,C.B);
while (est.f-est.i) > P.rtol*est.f && numel(R.probes) < P.max_probes
    D   = chord(C.s*est.i,C.s*est.f,C.B);
    tau = P.c_lean*D/2;
    fp  = firstpass(seedtr,tau);
    if tau < P.tau_min || (strcmp(P.solver,'cuadmm') && fp > P.kmax/2)
        R.flags{end+1} = sprintf('resolution: tau %.2e, seed first passage %g',tau,fp);
        break
    end
    cap = min(P.kmax,max(P.kmin*(D0/D)^P.p,2*fp));
    [R,stop] = mustStop(R,C,cap);  if stop, break; end
    if est.i > 0, g = sqrt(est.i*est.f); else, g = (est.i+est.f)/2; end
    pr  = probe(C,g,tau,round(cap),'probe',false,numel(R.probes)+1);
    R   = logprobe(R,pr,C);
    [cert,est,R] = update(cert,est,g,pr.verdict,R,P);
end

% ---- 7. endgame: certify just above the estimate; moves the certified bracket only
d = 0;
while est.f*(1+d) < cert.F && numel(R.probes) < P.max_probes && isempty(R.stopped)
    [R,stop] = mustStop(R,C,P.cert_iter);  if stop, break; end
    g  = est.f*(1+d);
    pr = probe(C,g,P.cert_tol,P.cert_iter,'cert',true,numel(R.probes)+1);
    R  = logprobe(R,pr,C);
    if strcmp(pr.verdict,'F'), cert.F = g;  break; end
    if strcmp(pr.verdict,'I')
        cert.I = max(cert.I,g);
        R.flags{end+1} = sprintf('endgame: certified infeasible at %.8g, above gamma_f',g);
    end
    if d == 0, d = P.delta; else, d = 2*d; end
end
R = finish(R,cert,est,C);
end


% =========================================================================
function [pr,tr] = probe(C,g,tol,cap,phase,wantF,n)
% One feasibility SDP at gamma = g (g = NaN: the unpinned feasibility program).
% F and I only through the verifiers.
P  = C.P;
if isempty(C.j)
    bscl = C.B;  bn = C.b_un/bscl;
else
    b2 = [C.b_un; C.s*g];  bscl = norm(b2);  bn = b2/bscl;
end
pr = newpr();  pr.phase = phase;  pr.gam = g;  pr.tau = tol;  pr.cap = cap;
tr = [];
switch P.solver
case 'cuadmm'
    write_b(fullfile(C.dir,'b.txt'),bn);
    for f = {'X_opt.txt','y_opt.txt','S_opt.txt'}      % a stale certificate must not survive
        if exist(fullfile(C.dir,f{1}),'file'), delete(fullfile(C.dir,f{1})); end
    end
    lg  = fullfile(C.dir,sprintf('run_%03d_%s.log',n,phase));
    to  = min([P.run_timeout, P.max_wall - toc(C.t0), P.deadline - posixtime(datetime('now'))]);
    pre = '';  if isfinite(to), pre = sprintf('timeout %d ',max(1,floor(to))); end
    cmd = sprintf(['wsl -e bash -c "LD_LIBRARY_PATH=/usr/lib/wsl/lib %s''%s'' ''%s/'' ' ...
                   '%.6g %d 15 1 > ''%s'' 2>&1"'],pre,P.exe,wslpath(C.dir),tol,cap,wslpath(lg));
    system(cmd);
    txt = fileread(lg);
    tm  = regexp(txt,'CUADMM_TIMING init (\S+) solve (\S+)','tokens','once');
    if isempty(tm)                                    % cuadmm_exe exits 0 on CUDA failure;
        pr.note = 'cuADMM: no CUADMM_TIMING line (failed or timed out)';   % a timeout lands here too
        return
    end
    pr.t_s  = str2double(tm{2});
    pr.conv = ~isempty(strfind(txt,'Solver ended: converged'));
    tr = parse_trace(txt);
    if ~isempty(tr)
        rm = cummin(tr(:,2));
        pr.iters = tr(end,1);  pr.pinf_end = tr(end,2);  pr.pinf_min = rm(end);
        pr.dobj_end = tr(end,5);
        [pr.slope,pr.aitken] = diag_trace(tr(:,1),rm);
        tr = [tr(:,1) rm];
    end
    % lean on the PRIMAL residual (the retro-tested rule), not cuADMM's own
    % convergence flag, which also waits for dinf and the gap: measured on
    % hinf_rd1 at 1.0009 gamma*, pinf reached 2e-5 < tau = 3.3e-5 while the dual
    % side lagged to the cap, and the flag gave a wrong (conservative) 'i'.
    pr.verdict = 'i';  if pr.conv || pr.pinf_min <= tol, pr.verdict = 'f'; end
    % CC, 09/26/2026: no lean from the start-up transient.  Measured on ctrl_rd1
    % (bracket [0, 4700], tau 3.5e-2): ten probes 'converged' at iteration 6-7
    % with pinf 2.8e-2 < tau, five of them below Mosek's certified-infeasible 128.
    % A certification (F) can still come out of a short run; only the lean is voided.
    if ~(pr.iters >= P.lean_kmin), pr.verdict = 'i';  pr.note = 'stopped in the start-up transient: no lean'; end
    xf_ = fullfile(C.dir,'X_opt.txt');  yf_ = fullfile(C.dir,'y_opt.txt');
    if wantF || (pr.conv && pr.pinf_end <= P.try_F)
        xs = readvec(xf_);
        [ok,info] = certify_primal(C,xs,bn,g,bscl);
        pr = putF(pr,info);
        if ok, pr.verdict = 'F'; end
    end
    if ~strcmp(pr.verdict,'F') && ~pr.conv
        y = readvec(yf_);
        [ok,info] = certify_farkas(C,y,bn);
        pr.I_beta = info.beta;  pr.I_eps = info.eps;
        if ok, pr.verdict = 'I'; end
    end
    % keep this run's certificate candidates for offline re-analysis
    kx = fullfile(C.dir,sprintf('run_%03d_X.txt',n));  ky = fullfile(C.dir,sprintf('run_%03d_y.txt',n));
    if exist(xf_,'file'), movefile(xf_,kx); end
    if exist(yf_,'file'), movefile(yf_,ky); end
case 'mosek'
    prob = Sedumi2Mosek(C.At2',bn,sparse(C.nvar,1),C.K);
    [~,res] = mosekopt('minimize info echo(0)',prob);
    pr.t_s  = res.info.MSK_DINF_OPTIMIZER_TIME;
    pro = res.sol.itr.prosta;  pr.note = pro;
    if ~isempty(strfind(pro,'INFEASIBLE'))
        [ok,info] = certify_farkas(C,res.sol.itr.y,bn);
        pr.I_beta = info.beta;  pr.I_eps = info.eps;
        if ok, pr.verdict = 'I'; end
    elseif ~isempty(strfind(pro,'FEASIBLE'))
        x  = MosekSol2SedumiSol(C.K,res);
        [ok,info] = certify_primal(C,C.T*x(:),bn,g,bscl);
        pr = putF(pr,info);
        pr.verdict = 'f';  if ok, pr.verdict = 'F'; end
    end
case 'sedumi'
    t1 = tic;
    [x,y,info] = sedumi(C.At2',bn,sparse(C.nvar,1),C.K,struct('fid',0));
    pr.t_s = toc(t1);
    pr.note = sprintf('sedumi pinf=%d dinf=%d numerr=%d',info.pinf,info.dinf,info.numerr);
    if info.pinf == 1
        [ok,inf2] = certify_farkas(C,y,bn);
        pr.I_beta = inf2.beta;  pr.I_eps = inf2.eps;
        if ok, pr.verdict = 'I'; end
    elseif info.dinf == 0
        [ok,inf2] = certify_primal(C,C.T*x(:),bn,g,bscl);
        pr = putF(pr,inf2);
        pr.verdict = 'f';  if ok, pr.verdict = 'F'; end
    end
otherwise
    error('bl_bisect:solver','unknown solver %s',P.solver);
end
end


function [ok,info] = certify_primal(C,xs,bn,g,bscl)
% Alternating projections: affine (min-norm correction with gamma fixed; lsqr
% from zero converges to it and ignores the numerically null directions) then
% clip each block's negative eigenvalues to 0.  F only if a point AFTER an
% affine step passes.  Blocks are rebuilt by reshape of T'x, as cuimport does,
% so an svec error that T'T would hide shows up as lost positivity.
P = C.P;  K = C.K;  x = xs(:);  nb = norm(bn);  ok = false;  lmin = -inf;  rres = inf;
for r = 1:P.rep_rounds
    [dx,~] = lsqr(C.As',C.As'*x - bn,1e-14,P.lsqr_it);
    x  = x - dx;  xa = x;                      % xa: last affine-repaired point, what F is judged on
    rres = norm(C.As'*x - bn)/nb;
    xf = C.T'*x;
    [lmin,V,L,lmax] = blockeig(xf,K,C.Ks);
    labs = lmin;  lmin = lmin/max(lmax,realmin);   % relative to the largest eigenvalue over blocks
    if rres <= P.rep_tol && lmin >= -P.psd_tol, ok = true; break; end
    off = K.f;
    for k = 1:numel(C.Ks)
        N = C.Ks(k);
        Xk = V{k}*diag(max(L{k},0))*V{k}';
        xf(off+(1:N^2)) = Xk(:);
        off = off + N^2;
    end
    x = C.T*xf;
end
% independent of the pinned system: the unpinned rows and x_j = gamma, in the
% units the executive posed them (undo bl_fix's 1/bscl)
xt = (C.T'*xa)*bscl;                           % not x: after a failed round x is the clipped point
unpin = norm(C.At0'*xt - C.b_un)/max(norm(C.b_un),eps);
pin   = 0;  if ~isempty(C.j), pin = abs(xt(C.j) - g)/max(abs(g),eps); end
ok = ok && unpin <= 10*P.rep_tol && pin <= 10*P.rep_tol;
% scale guards (see header, item 1): a near-infeasible program's huge-norm point
% flatters the relative test, so bound the absolute negativity and the norm
ok = ok && labs >= -P.psd_abs;
if isfield(C,'normx_ref') && ~isempty(C.normx_ref), ok = ok && norm(xt) <= P.norm_ratio*C.normx_ref; end
% lmax and ||X|| are logged so a ratio flattered by a large lambda_max (a drifting
% recession direction) is visible against Mosek's; labs is the absolute margin,
% to set against eppos on the c = 0 classes
info = struct('lmin',lmin,'res',rres,'unpin',max(unpin,pin),'lmax',lmax, ...
              'labs',labs,'normx',norm(xt));
end


function [ok,info] = certify_farkas(C,y,bn)
% Farkas: At2*y zero on the free block and NSD on every PSD block, b'y > 0 =>
% no X in the cone meets At2'x = b.  Computed from the SeDuMi At (asymmetric A_k
% symmetrised here), not from the svec system the solver saw.  Both signs are
% tried: a certificate is valid whichever sign passes.
K = C.K;  Kf = K.f;  ok = false;  info = struct('beta',NaN,'eps',NaN);
best = -inf;
for s = [1 -1]
    yy = s*y(:);
    if Kf > 0                                  % exact free-block correction, twice
        F = C.At2(1:Kf,:);  FF = F*F';  FF = FF + 1e-14*max(1,max(abs(diag(FF))))*speye(Kf);
        for t = 1:2, yy = yy - F'*(FF\(F*yy)); end
        fres = norm(F*yy);
    else
        fres = 0;
    end
    z  = C.At2*yy;  ny = norm(yy);  beta = bn'*yy;
    lmax = -inf;  off = Kf;
    for k = 1:numel(C.Ks)
        N = C.Ks(k);  Zk = reshape(z(off+(1:N^2)),N,N);
        lmax = max(lmax,max(eig((Zk+Zk')/2)));
        off = off + N^2;
    end
    good = beta > 0 && lmax <= 0 && fres <= 1e-12*max(norm(z),eps);
    score = beta/ny - max(lmax,0)/ny;
    if good || score > best
        best = score;  info = struct('beta',beta/ny,'eps',max(lmax,0)/ny);
        ok = good;
    end
    if good, return; end
end
end


function [lmin,V,L,lmax] = blockeig(xf,K,Ks)
lmin = inf;  lmax = -inf;  off = K.f;  V = cell(1,numel(Ks));  L = V;
for k = 1:numel(Ks)
    N = Ks(k);  Xk = reshape(xf(off+(1:N^2)),N,N);
    [V{k},Dk] = eig((Xk+Xk')/2);  L{k} = diag(Dk);
    lmin = min(lmin,min(L{k}));  lmax = max(lmax,max(L{k}));
    off = off + N^2;
end
end


function T = svecmap(Kf,Ks,vec_len,nvar)
% SDPT3 svec, identical to cuimport.m and dump2cuadmm.m (local copies there).
r2 = sqrt(2)/2;  I = {};  J = {};  V = {};
if Kf > 0
    p = (1:Kf)';  I{end+1} = p;  J{end+1} = p;  V{end+1} = ones(Kf,1);
end
bs = Kf;  bv = Kf;
for k = 1:numel(Ks)
    N = Ks(k);
    idx = find(triu(true(N)));
    [ii,jj] = ind2sub([N N],idx);
    p = (1:numel(idx))';  dg = (ii==jj);  od = ~dg;
    I{end+1} = [bs+p;   bs+p(od)];                                   %#ok<AGROW>
    J{end+1} = [bv+idx; bv+sub2ind([N N],jj(od),ii(od))];            %#ok<AGROW>
    V{end+1} = [dg + od*r2; repmat(r2,nnz(od),1)];                   %#ok<AGROW>
    bs = bs + N*(N+1)/2;  bv = bv + N^2;
end
T = sparse(cat(1,I{:}),cat(1,J{:}),cat(1,V{:}),vec_len,nvar);
end


function D = chord(a,b,B)
% distance between the unit-normalised pinned b at gamma = a and gamma = b
D = 2*sin(abs(atan(a/B) - atan(b/B))/2);
end


function k = firstpass(tr,tau)
% iterations the certified seed needed for its running-min pinf to reach tau
if isempty(tr), k = 0; return; end
i = find(tr(:,2) <= tau,1);
if isempty(i), k = inf; else, k = tr(i,1); end
end


function tr = parse_trace(txt)
% rows ' iter | pinf dinf | pobj dobj relgap | t | sigma'
tok = regexp(txt,['^\s*(\d+)\s*\|\s*(\S+)\s+(\S+)\s*\|\s*(\S+)\s+(\S+)\s+(\S+)\s*\|'], ...
             'tokens','lineanchors');
tr = zeros(numel(tok),6);
for i = 1:numel(tok), tr(i,:) = str2double(tok{i}); end
tr = tr(all(isfinite(tr),2),:);
end


function [s,a] = diag_trace(k,rm)
% recorded, never used to steer: last-octave log-log slope of the running min,
% and the Aitken plateau from (k/4, k/2, k)
s = NaN;  a = NaN;  ke = k(end);
i2 = find(k <= ke/2,1,'last');  i4 = find(k <= ke/4,1,'last');
if ~isempty(i2) && k(i2) > 0, s = log(rm(end)/rm(i2))/log(ke/k(i2)); end
if ~isempty(i4) && ~isempty(i2)
    p1 = rm(i4);  p2 = rm(i2);  p3 = rm(end);  dn = p1 + p3 - 2*p2;
    if dn > 0, a = (p1*p3 - p2^2)/dn; end
end
end


function [cert,est,R] = update(cert,est,g,v,R,P)
switch v
    case 'F', cert.F = min(cert.F,g);  est.f = min(est.f,g);
    case 'f', est.f = min(est.f,g);
    case 'I', cert.I = max(cert.I,g);  est.i = max(est.i,g);
    case 'i', est.i = max(est.i,g);
end
% feasibility is monotone in gamma, so these orderings are forced
if cert.I >= cert.F
    R.flags{end+1} = sprintf('ORDER: certified I %.8g >= certified F %.8g -- verifier or pinning bug',cert.I,cert.F);
end
if strcmp(v,'F') && g < P.lo
    R.flags{end+1} = sprintf('ORDER: F at %.8g below the analytic bound %.8g',g,P.lo);
end
if est.i >= est.f                           % a lean contradicted: fall back to certified
    R.flags{end+1} = sprintf('lean conflict at %.8g; estimate reset to certified bracket',g);
    est.f = cert.F;  est.i = cert.I;
end
end


function R = logprobe(R,pr,C)
R.probes(end+1) = pr;
R.t_solver = R.t_solver + max(pr.t_s,0);
if isfinite(pr.iters), R.iters = R.iters + pr.iters; end
if isfinite(pr.iters) && pr.iters > 0 && isfinite(pr.t_s), R.spi = pr.t_s/pr.iters; end
% append to the live tsv now: a crash later must not lose this probe
fl = fieldnames(pr);  c = cell(1,numel(fl));
for q = 1:numel(fl), c{q} = cellval(pr.(fl{q})); end
fid = fopen(C.live,'a');  fprintf(fid,'%s\n',strjoin(c,sprintf('\t')));  fclose(fid);
fprintf('BIS %-5s gam=%-13.8g tau=%-9.2e cap=%-6d it=%-6g t=%7.2fs pinf=%-9.2e -> %s  (Flmin %.1e Fres %.1e | Ibeta %.1e Ieps %.1e) %s\n', ...
    pr.phase,pr.gam,pr.tau,pr.cap,pr.iters,pr.t_s,pr.pinf_end,pr.verdict, ...
    pr.F_lmin,pr.F_res,pr.I_beta,pr.I_eps,pr.note);
end


function pr = newpr()
% one probe record; field order is the tsv column order
pr = struct('phase','','gam',NaN,'tau',NaN,'cap',NaN,'iters',NaN,'t_s',NaN, ...
    'conv',false,'pinf_end',NaN,'pinf_min',NaN,'dobj_end',NaN,'verdict','i', ...
    'F_lmin',NaN,'F_res',NaN,'F_unpin',NaN,'F_lmax',NaN,'F_labs',NaN,'F_normx',NaN, ...
    'I_beta',NaN,'I_eps',NaN,'slope',NaN,'aitken',NaN,'note','');
end


function pr = putF(pr,info)
pr.F_lmin = info.lmin;  pr.F_res = info.res;  pr.F_unpin = info.unpin;
pr.F_lmax = info.lmax;  pr.F_labs = info.labs;  pr.F_normx = info.normx;
end


function [R,stop] = mustStop(R,C,cap)
% Kill switches, checked before every solver call: STOP file, the run's wall
% budget, and the regime deadline, each against the call's estimated length
% (cap x the last measured seconds per iteration, plus 10 s of overhead).
P = C.P;  stop = false;
est = 10;  if strcmp(P.solver,'cuadmm'), est = 10 + cap*R.spi; end
if exist(P.stopfile,'file'), R.stopped = 'STOP file';
elseif toc(C.t0) + est > P.max_wall, R.stopped = sprintf('max_wall %.0f s',P.max_wall);
elseif posixtime(datetime('now')) + est > P.deadline, R.stopped = 'deadline';
end
if ~isempty(R.stopped)
    stop = true;  R.flags{end+1} = ['STOPPED: ' R.stopped];
end
end


function s = cellval(x)
if ischar(x), s = x; elseif islogical(x), s = num2str(double(x)); else, s = num2str(x,'%.8g'); end
end


function R = finish(R,cert,est,C)
R.gamma_F = cert.F;  R.gamma_I = cert.I;  R.gamma_f = est.f;  R.gamma_i = est.i;
v = {R.probes.verdict};
R.counts = struct('F',nnz(strcmp(v,'F')),'f',nnz(strcmp(v,'f')), ...
                  'i',nnz(strcmp(v,'i')),'I',nnz(strcmp(v,'I')));
fprintf(['BISDONE %s [%s]  certified [%.8g, %.8g]  estimate [%.8g, %.8g]  ' ...
         'probes %d (F%d f%d i%d I%d)  solver %.1f s  iters %d\n'], ...
    R.id,R.solver,cert.I,cert.F,est.i,est.f,numel(v),R.counts.F,R.counts.f, ...
    R.counts.i,R.counts.I,R.t_solver,R.iters);
for i = 1:numel(R.flags), fprintf('  FLAG %s\n',R.flags{i}); end
tsv = fullfile(C.P.outdir,sprintf('%s_%s_bisect.tsv',R.id,R.solver));
fl  = fieldnames(R.probes);
fid = fopen(tsv,'w');  fprintf(fid,'%s\n',strjoin(fl',sprintf('\t')));
for i = 1:numel(R.probes)
    c = cell(1,numel(fl));
    for q = 1:numel(fl), c{q} = cellval(R.probes(i).(fl{q})); end
    fprintf(fid,'%s\n',strjoin(c,sprintf('\t')));
end
fclose(fid);
R.tsv = tsv;
end


function write_b(fn,bn)
% cuADMM's b.txt: 0-based COO, ascending rows (dump2cuadmm's store_coo format)
[r,~,v] = find(sparse(bn(:)));
fid = fopen(fn,'w');
fprintf(fid,'%d %d %.17g\n',[r(:)'-1; zeros(1,numel(r)); v(:)']);
fclose(fid);
end


function v = readvec(fn)
fid = fopen(fn,'r');
if fid < 0, error('bl_bisect:open','cannot open %s',fn); end
v = fscanf(fid,'%f');  fclose(fid);
end


function p = wslpath(p)
p = strrep(p,'\','/');
if numel(p) >= 2 && p(2) == ':', p = ['/mnt/' lower(p(1)) p(3:end)]; end
end
