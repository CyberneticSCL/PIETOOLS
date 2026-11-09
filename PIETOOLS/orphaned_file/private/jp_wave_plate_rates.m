function R = jp_wave_plate_rates(sys,par,opts)
% R = JP_WAVE_PLATE_RATES(SYS,PAR,OPTS) certified lower bound on the
% exponential PIE-to-PDE decay rate of the damped 2D wave or the clamped 2D
% plate of Jagt & Peet, arXiv:2508.14840 (Sec. 7.2.2-7.2.3), by a top-down
% search on k (default; monotone in k, see 'topdown' below) or bisection, in
% the LPI of Cor. 35 of that paper, written exactly as the Sec. 7.1 listing:
%
%   P = lpivar(T.dim,d),  P'*T = T'*P,  T'*P - ep^2*T'*T >= 0,
%   -(P'*A + A'*P + 2*k*P'*T) >= 0,
%
% optionally with the linear face generators added to both inequalities
% (OPTS.psatz = [3 4 5 6]; not in Cor. 35 as printed -- on this tree the
% printed listing certifies no positive rate for the heat example of the
% paper, measured in PIETOOLS_demos/sopvar_demos/README.md).
%
%   SYS    'wave' (PAR = kap, the exact rate is kap) | 'plate' (PAR = alp0)
%   OPTS   d (1 wave / 0 plate, as in the paper; [] = the lpivar default
%          degrees = the tensor Z_1 of (14) -- measured: the scalar d = 1 of
%          the listing is certified infeasible at k = 0 on the wave, Z_1 is
%          not), ep (0.1), psatz (0), mode ('topdown' | 'bisect'), fracs and
%          above (topdown: probes hi*(1-fracs) and (1+above)*hi; above = []
%          skips the latter), hi (wave: kap, the exact rate; plate: the value
%          in the paper), rtol (1e-4), maxsteps (24, bisect), deadline (posix
%          s, Inf), tight (true; heatNd_solve only), accept ('certified' | 'psd_clip' | 'mosek':
%          the MOSEK status, NOT certified; topdown only), scale (1; wave
%          only: state u2 = phi_t/scale)
% R      lo (largest certified k; NaN if none); hi (topdown: the top of the
%        probe grid; bisect: smallest k certified infeasible or uncertain
%        above lo); I_above (topdown: (1+above)*hi if certified infeasible,
%        else NaN); trace [k st t_mosek rel_b] (column 4 is eta of the
%        clipped point under 'psd_clip'); affine_err; the SDP shape, build and
%        solve times; the value in the paper, for comparison.
%
% VERDICTS are certified unless OPTS.accept = 'mosek' (heatNd_solve: +1 needs
% rel_b <= 1e-6, PSD blocks and x ~= 0; -1 a verified Farkas ray; else 0,
% 'uncertain', which bounds the steering but never the certified result).
% 'psd_clip' certifies with bl_bisect 'feas' and MOSEK (default tolerances)
% instead: +1 iff MOSEK reports a FEASIBLE status and its point, repaired onto
% the rows and clipped to PSD, has row-normwise backward error eta <= 1e-7, with
% lambda_min >= -1e-6 (unit-norm b) before clipping; -1 a verified Farkas ray.
% Rule: PIETOOLS_demos/cuadmm/bl_bisect.m; its calibration (certify_primal) is
% on other classes, and on the wave no infeasible k returned a FEASIBLE status,
% so the +1 side has no infeasible control on this class.
% On the same SDP (wave kap 3, scale 1, k = 3) heatNd_solve (tight) rejected
% its MOSEK point, rel_b 5.9e-6, and bl_bisect certified its own MOSEK point,
% eta 1.0e-10 (two solves, measured 10/06/2026).
% TABLE 2 (wave): jp_wave_plate_rates('wave',kap,struct('d',[],'psatz',[3 4 5 6],
% 'scale',kap,'accept','psd_clip')) certified k = kap exactly (eta 1.1e-10 to
% 2.3e-10) and 1.01*kap infeasible for kap = 1..7, 21-27 s per kap (15-20 s
% building, 3.5-4.9 s in the two MOSEK probes) (measured with this script
% 10/07/2026, before the eta print label, which was checked at kap = 3; paper
% 0.9999 .. 6.9645). Without the scaling, 'psd_clip' certifies k = kap for
% kap = 1..4 only.
% THE SDP IS AFFINE IN k (only the 2k*P'*T term depends on it): built at
% k = 1 and k = 2 and interpolated, At(k) = At1 + (k-1)*(At2 - At1); checked
% against a direct build at k = 0.5 before any solve.
% DEPENDENCY: heatNd_path / heatNd_sdp / heatNd_solve (MMP, 09/27/2026) in
% PIETOOLS_demos/sopvar_demos; for 'psd_clip', bl_bisect and cuadmm_path in
% PIETOOLS_demos/cuadmm.
%
% CC, 09/28/2026: initial coding, recreating the wave and plate examples of
% the paper (PIETOOLS_examples/Examples_Library/2D/PIETOOLS_PDE_Ex_2D_Wave_Eq_Damped.m
% and PIETOOLS_PDE_Ex_2D_Plate_Eq.m).
% CC, 09/28/2026: k-affine basis built at k = 1, 2 instead of 0, 1. The plate
% k = 0 build has fewer constraints (monomials reached only via 2k*P'T), so the
% affine check errored ("incompatible sizes") before any solve.
% CC, 09/28/2026: OPTS.above = [] skips the probe above hi (plate solves > 17 min).
% CC, 09/28/2026: the heatNd tools moved to PIETOOLS_demos/sopvar_demos; the
% addpath follows them. Header corrected to the code (top-down default,
% R.hi per mode, R.I_above, accept, scale).
% CC, 10/07/2026: OPTS.accept = 'psd_clip' (bl_bisect's certificate rule, with
% MOSEK). The heatNd_solve gate rejected k = kap at kap = 3, 4 (scale 1) and
% kap = 5, 6, 7 (scale = kap), where this rule certifies it (separate MOSEK
% solves); with scale = kap it reproduces Table 2 (kap = 1..7, exactly k = kap).
% The kap = 3 note under 'accept' below is the heatNd_solve gate only.

if nargin < 3, opts = struct(); end
isw = strcmp(sys,'wave');
df = struct('d',1-~isw*1,'ep',0.1,'psatz',0,'rtol',1e-4,'hi',[],'maxsteps',24, ...
            'deadline',Inf,'tight',true,'mode','topdown', ...
            'fracs',[0 1e-4 1e-3 1e-2 0.05 0.1 0.2 0.5 0.9],'above',0.01,'accept','certified','scale',1);
f = fieldnames(opts);  for i = 1:numel(f), df.(f{i}) = opts.(f{i}); end
o = df;  if ~isw && ~isfield(opts,'d'), o.d = 0; end
here = fileparts(mfilename('fullpath'));
% addpath(fullfile(here,'..','..','sopvar','Testfolder','sdopvar','claude_tests','heatNd')); % CC, 09/28/2026 (was)
addpath(fullfile(here,'..','sopvar_demos'));                                % CC, 09/28/2026
% the path: this checkout only.  A genpath of the checkout also picks up
% nested .claude/worktrees copies, and one of them shadowed poslpivar_2d with
% a version lacking psatz 3-6 (measured, 09/28/2026); heatNd_path strips them
heatNd_path();
ROOT = fileparts(fileparts(here));
for r = {'convert','poslpivar_2d','lpi_ineq_2d','lpivar_2d','lpiprogram','lpi_eq_2d', ...
         'PIETOOLS_PDE_Ex_2D_Wave_Eq_Damped','PIETOOLS_PDE_Ex_2D_Plate_Eq'}
    % every hit inside this checkout, and at most one that is not a class
    % method (convert is both @pde_struct/convert and @sys/convert: one owner each)
    w = unique(which(r{1},'-all'));
    plain = w(~contains(w,[filesep '@']));
    if isempty(w) || ~all(strncmpi(w,ROOT,numel(ROOT))) || numel(plain) > 1
        error('jp_wave_plate_rates:path','%s resolves as:\n  %s\n(ROOT %s)',r{1},strjoin(w,'\n  '),ROOT);
    end
end
if isw, paper = [0.9999 1.9998 2.9874 3.9973 4.9980 5.9906 6.9645];
        R.paper = NaN;  if par == round(par) && par >= 1 && par <= 7, R.paper = paper(par); end
else,   R.paper = NaN;  if abs(par-0.2) < 1e-12, R.paper = 3.6328; end
end

% ---- the PIE and the k-affine SDP
t0 = tic;
if isw && o.scale ~= 1
    % wave with the second state scaled, u = [phi; phi_t/s]: the diagonal
    % similarity S = diag(1,1/s) maps any Cor. 35 certificate P to S'PS and
    % leaves the rate unchanged, but the operator entries go from O(kap^2)
    % to O(kap) for s = kap (a conditioning test, not the form in the paper)
    pvar s1 s2;  s = o.scale;  kap = par;
    clear stateNameGenerator
    u = pde_var(2,[s1;s2],[0,1;0,1],[2;2]);
    PDE = [diff(u,'t')==[0,s;-kap^2/s,-2*kap]*u + [0,0;1/s,0]*(diff(u,s1,2)+diff(u,s2,2));
           subs(u,s1,0)==0;  subs(diff(u,s1),s1,1)==0;  subs(u,s2,0)==0;  subs(diff(u,s2),s2,1)==0];
elseif isw, PDE = PIETOOLS_PDE_Ex_2D_Wave_Eq_Damped(0,{sprintf('kap=%.17g;',par)});
else,       PDE = PIETOOLS_PDE_Ex_2D_Plate_Eq(0,{sprintf('alp0=%.17g;',par)}); end
PIE = convert(PDE,'pie');
% basis at k = 1, 2: the k = 0 build drops the constraints reached only via
% 2k*P'T (plate d = 0: 135 of 44,239), so it cannot be subtracted from k = 1.
% The k = 0 member of the family is that build plus those constraints as 0 = 0
% (measured: At equal to 9.2e-16, b equal)
% D0 = build(PIE,o,0);  D1 = build(PIE,o,1);  Dh = build(PIE,o,0.5);        % CC, 09/28/2026 (was)
D1 = build(PIE,o,1);  D2 = build(PIE,o,2);  Dh = build(PIE,o,0.5);          % CC, 09/28/2026
R.t_build = toc(t0);
% size before b: the old order errored on norm(D0.b-D1.b) for the plate
% if ~isequal(D0.K,D1.K) || norm(D0.b-D1.b) > 0 || ~isequal(size(D0.At),size(D1.At)) % CC, 09/28/2026 (was)
if ~isequal(D1.K,D2.K) || ~isequal(size(D1.At),size(D2.At)) || norm(D1.b-D2.b) > 0 % CC, 09/28/2026
    error('jp_wave_plate_rates:affine','the k = 1 and k = 2 SDPs differ in shape or b'); % CC, 09/28/2026
end
% dA = D1.At - D0.At;                                                       % CC, 09/28/2026 (was)
dA = D2.At - D1.At;                                                         % CC, 09/28/2026
D0 = setfield(D1,'At',D1.At - dA);  clear D2                                % CC, 09/28/2026 %#ok<SFLD>
aff = full(max(abs(nonzeros(D0.At + 0.5*dA - Dh.At))));  if isempty(aff), aff = 0; end
R.affine_err = aff;
if aff > 1e-10*max(1,full(max(abs(nonzeros(D0.At)))))
    error('jp_wave_plate_rates:affine','interpolated SDP at k = 0.5 differs from the build by %.2e',aff);
end
R.m = D0.m;  R.nx = D0.nx;  R.nnz = D0.nnz;  R.Ks = D0.K.s(:)';
fprintf('JP %s par=%g d=%d psatz=%s: m=%d nnz=%d Ks=%s, build %.1fs (affine err %.1e)\n', ...
    sys,par,o.d,mat2str(o.psatz),R.m,R.nnz,mat2str(R.Ks),R.t_build,aff);

% ---- 'topdown' (default): probe k = hi*(1 - f) for f in o.fracs, largest k
% first, and stop at the first CERTIFIED feasible point.  Sound because Cor. 35
% is monotone in k: its constraints give P'T >= ep^2 T'T >= 0, so
% -(P'A + A'P + 2*k1*P'T) = -(P'A + A'P + 2k P'T) + 2(k - k1)P'T >= 0 for
% k1 <= k with the same P.  Also probes k = (1 + o.above)*hi for a certified
% infeasible point.  Why not bisection: measured on the wave (kap = 1, tensor
% Z_1, faces), MOSEK returns UNKNOWN at every probed k in {0.8, 0.85, 0.9,
% 0.95, 0.99, 0.999} but a certified +1 at k = 1 and a certified -1 at
% k = 1.01; a bisection that steers down on 'uncertain' never tests k = 1
% and reported 0.7197.
sv = struct('tight',o.tight);
% o.accept 'mosek' (NOT certified; for comparing with the paper, which reports
% the MOSEK status): +1 whenever MOSEK returns PRIMAL_AND_DUAL_FEASIBLE/OPTIMAL,
% whatever rel_b.  Measured at kap = 3 (faces): no probe certified -- MOSEK
% OPTIMAL but rel_b 1.6e-6 > 1e-6 at k = 0.3 (eta 1e-8), UNKNOWN at the
% probes k = 1.5, 2.7, 2.97, 3.  'psd_clip' certifies k = 3 there (header).
if strcmp(o.accept,'mosek')
    hs = @heatNd_solve;
    solve = @(D,s) mosek_status(hs(D,s));
elseif strcmp(o.accept,'psd_clip')                                          % CC, 10/07/2026
    cuadmm_path();                  % adds MOSEK and SeDuMi; same checkout  % CC, 10/07/2026
    solve = @(D,s) psd_clip_solve(D);                                       % CC, 10/07/2026
else
    solve = @heatNd_solve;
end
if strcmp(o.mode,'topdown')
    if isw, hi = par; else, hi = o.hi; if isempty(hi), hi = R.paper; end, end
    R.trace = zeros(0,4);  R.lo = NaN;  R.I_above = NaN;
    at = @(k) setfield(D0,'At',D0.At + k*dA);                               %#ok<SFLD>
    % o.above = [] skips the probe: a plate solve exceeds 17 min (measured)
    % and the paper claims only a lower bound there, so k = hi goes first
    if ~isempty(o.above)                                                    % CC, 09/28/2026
    v = solve(at((1+o.above)*hi),sv);  R.trace(end+1,:) = [(1+o.above)*hi v.st v.t_mosek v.rel_b];
    if v.st == -1, R.I_above = (1+o.above)*hi; end
    fprintf('JP %s k=%.6f (above): st=%+d %s\n',sys,(1+o.above)*hi,v.st,v.why);
    end                                                                     % CC, 09/28/2026
    rlab = 'rel_b';  if strcmp(o.accept,'psd_clip'), rlab = 'eta'; end      % CC, 10/07/2026
    for fr = o.fracs
        if posixtime(datetime('now')) > o.deadline, break; end
        k = hi*(1-fr);
        v = solve(at(k),sv);  R.trace(end+1,:) = [k v.st v.t_mosek v.rel_b];
%       fprintf('JP %s k=%.6f (1-%g): st=%+d %s rel_b=%.1e t=%.0fs\n',sys,k,fr,v.st,v.prosta,v.rel_b,v.t_mosek); % CC, 10/07/2026 (was)
        fprintf('JP %s k=%.6f (1-%g): st=%+d %s %s=%.1e t=%.0fs\n',sys,k,fr,v.st,v.prosta,rlab,v.rel_b,v.t_mosek); % CC, 10/07/2026
        if v.st == 1, R.lo = k; break; end
    end
    R.hi = hi;  R.why = 'topdown';  R.t_total = toc(t0);
    fprintf('JP %s par=%g: certified rate >= %.6f; certified infeasible at %.6f; paper %.4f; %d solves, %.0fs\n', ...
        sys,par,R.lo,R.I_above,R.paper,size(R.trace,1),R.t_total);
    return
end
% ---- 'bisect': bisection on certified verdicts
R.trace = zeros(0,4);
at = @(k) setfield(D0,'At',D0.At + k*dA);                                   %#ok<SFLD>
v = heatNd_solve(at(0),sv);  R.trace(end+1,:) = [0 v.st v.t_mosek v.rel_b];
lo = NaN;  hi = Inf;
if v.st ~= 1
    R.lo = NaN;  R.hi = 0;  R.why = sprintf('k = 0 not certified (st %d: %s)',v.st,v.why);
    fprintf('JP %s: %s\n',sys,R.why);  return
end
lo = 0;
if isw, hi = par; else, hi = o.hi; if isempty(hi), hi = 8; end, end
if ~isw                                          % plate: grow hi while certified
    while true
        v = heatNd_solve(at(hi),sv);  R.trace(end+1,:) = [hi v.st v.t_mosek v.rel_b];
        if v.st == 1, lo = hi; hi = 2*hi; else, break; end
        if posixtime(datetime('now')) > o.deadline, break; end
    end
end
for s = 1:o.maxsteps
    if (hi-lo) <= o.rtol*max(lo,1) || posixtime(datetime('now')) > o.deadline, break; end
    k = (lo+hi)/2;
    v = heatNd_solve(at(k),sv);  R.trace(end+1,:) = [k v.st v.t_mosek v.rel_b];
    if v.st == 1, lo = k; else, hi = k; end      % uncertain steers down, never certifies
    fprintf('JP %s step %2d: k=%.6f st=%+d  [%.6f, %.6f]  t=%.0fs rel_b=%.1e\n',sys,s,k,v.st,lo,hi,v.t_mosek,v.rel_b);
end
R.lo = lo;  R.hi = hi;  R.why = 'bisection';
R.t_total = toc(t0);
fprintf('JP %s par=%g: certified rate >= %.6f (bracket [%.6f, %.6f]); paper %.4f; %d solves, %.0fs\n', ...
    sys,par,lo,lo,hi,R.paper,size(R.trace,1),R.t_total);
end


function D = build(PIE,o,k)
% the Sec. 7.1 listing at rate k, returned as the SDP sossolve would build
T = PIE.T;  A = PIE.A;
prog = lpiprogram(PIE.vars(:,1),PIE.vars(:,2),PIE.dom);
% d = [] uses the lpivar_2d default degrees, which (read from source) are the
% tensor basis Z_1 of (14) in the paper; a scalar d caps the joint degree
if isempty(o.d), [prog,P] = lpivar(prog,T.dim);
else,            [prog,P] = lpivar(prog,T.dim,o.d); end
prog = lpi_eq(prog,P'*T-T'*P);
io = struct('psatz',o.psatz);
prog = lpi_ineq(prog,T'*P-o.ep^2*(T'*T),io);
prog = lpi_ineq(prog,-(P'*A+A'*P+2*k*(P'*T)),io);
D = heatNd_sdp(prog);
end


% CC, 10/07/2026 (start): OPTS.accept = 'psd_clip'
function v = psd_clip_solve(D)
% one probe through bl_bisect 'feas' with MOSEK (its rule, not re-implemented
% here): the SDP goes to a temporary dump (unit-norm b, bscl in the _meta file,
% the format bl_bisect reads); the folder, including the cuADMM text copy of At
% that 'feas' mode writes (2.53 GB measured for the plate with faces at d = 0,
% 10/06/2026), is deleted on return, error or Ctrl+C (onCleanup).  bl_bisect
% gets its own stop file: its default (cuadmm_outdir()/STOP) is the harness
% kill switch, which would leave it with no probe
tmp = tempname;  mkdir(tmp);  cl = onCleanup(@() rmdir(tmp,'s'));  file = fullfile(tmp,'jp.mat');
S = struct('At',D.At,'b',D.b,'c',D.c,'K',D.K,'Ns',D.K.s(:)','Kf',D.K.f);
save(file,'-struct','S','-v7.3');  clear S
bscl = D.bscl;  save(fullfile(tmp,'jp_meta.mat'),'bscl');
R = bl_bisect(file,struct('solver','mosek','mode','feas','outdir',tmp,'tag','jp', ...
                          'stopfile',fullfile(tmp,'STOP')));
if isempty(R.probes), error('jp_wave_plate_rates:psd_clip','bl_bisect made no probe: %s',R.stopped); end
p = R.probes(end);
v.st = 0;  if p.verdict == 'F', v.st = 1; elseif p.verdict == 'I', v.st = -1; end
v.t_mosek = p.t_s;  v.rel_b = p.F_etapsd;  v.prosta = p.note;  v.solsta = p.verdict;
v.why = sprintf('psd_clip %s (eta %.1e, lambda_min %.1e; %s)',p.verdict,p.F_etapsd,p.F_labs,p.note);
end
% CC, 10/07/2026 (end)


function v = mosek_status(v)
% NOT a certificate: the MOSEK feasible status, whatever rel_b (o.accept 'mosek')
if strcmp(v.prosta,'PRIMAL_AND_DUAL_FEASIBLE') && strcmp(v.solsta,'OPTIMAL'), v.st = 1;
elseif v.st ~= -1, v.st = 0; end
v.why = ['MOSEK status (not certified): ' v.prosta '/' v.solsta];
end
