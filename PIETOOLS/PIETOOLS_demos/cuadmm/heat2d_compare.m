%HEAT2D_LOWRANK_MOSEK_COMPARE Compare LPI solution via:
% 1) Mosek interior-point method for SDPs
% 2) low-rank 
% 3) cuADMM 
% Uses the same PDE and LPI as PIETOOLS_demos/lowrank_2d_stability:
%
%   PIE  = 2D heat equation below
%   st   = 2D PIETOOLS settings below
%   LPI  = build_stab_2d_st2(PIE,st)
%
% PDE:
%   x_t = x_s1s1 + x_s2s2 + lam*x,    (s1,s2) in [0,1]^2, x in R^n
%   x = 0 on all four edges.
%
% Analytic stability threshold:
%   lam_star = 2*pi^2.
%
% Edit the variables in the parameters block below, then run this script.
% To assemble only, set run_mosek=false and run_lowrank=false, then comment
% out the cuADMM block.
%
% cuADMM:
%   The assembled LPI is exported in the text format used by
%   PIETOOLS_demos/cuadmm/dump2cuadmm.m.
% Note: the full Mosek reference for n=1, nd=4 is known to be a long solve.
%% parameters
n = 1;                          % Number of coupled heat-equation states.
lam = 2;                        % Reaction coefficient in x_t = x_s1s1 + x_s2s2 + lam*x.
nd = 4;                        % Polynomial degree offset for the negativity Gram.
dmult = [];                     % Optional multiplier degree data; [] uses PIETOOLS defaults.
run_mosek = true;               % Run the interior-point reference solve.
run_lowrank = true;             % Run the low-rank stability search.
gate = 1e-6;                    % Tolerance for accepting a certified stability result.
gate_k = 10;                    % Number of randomized/probe checks used by the gate.
gate_refmax = 1e-5;             % Maximum allowed reference residual used by the gate.
maxrank = 6;                    % Largest face/rank width the low-rank search may try.
rank = [];                      % Fixed face/rank width; [] lets the routine search up to maxrank.
seeds = [11 22 33];             % Random seeds for repeated low-rank discovery attempts.
lmit = 400;                     % Iteration limit for the low-rank nonlinear search.
lmtol = [];                     % Low-rank optimizer tolerance; [] uses the routine default.
cuadmm_tol = 1e-6;              % Target feasibility residual for the cuADMM probe.
cuadmm_maxiter = 50000;         % Maximum ADMM iterations for the cuADMM feasibility probe.
cuadmm_tag = 'heat2d';          % Label added to cuADMM output folders and probe tables.

%% path and solver defaults
pietools_root = '';            
solver = 'mosek';               % Interior-point solver name passed to the reference routine.
cuadmm_exe = '';                % Empty means use CUADMM_EXE/autodetect.
cuadmm_outdir = '';             % Empty means write under this script's working folder.
cuadmm_timeout = Inf;           % Wall-clock limit for one cuADMM executable call.
cuadmm_max_wall = Inf;          % Wall-clock limit for the whole cuADMM wrapper call.
cuadmm_pinf_norm = '2';         % Residual norm used by cuADMM stopping logic: '2' or 'inf'.
cuadmm_require_full_gpu = false;% If true, reject MIG slices for timing runs.
cuadmm_launcher = 'auto';       % cuADMM launch mode; 'auto' chooses Linux on SOL and WSL on Windows.
verbose = true;                 % Print solver progress and comparison diagnostics.

if isempty(pietools_root)
    env_pietools_root = getenv('PIETOOLS_ROOT');
    if ~isempty(env_pietools_root) && exist(fullfile(env_pietools_root, 'pietools_path_update.m'), 'file') == 2
        pietools_root = env_pietools_root;
    elseif exist('/scratch/mpeet/pietools/harness/PIETOOLS/pietools_path_update.m', 'file') == 2
        pietools_root = '/scratch/mpeet/pietools/harness/PIETOOLS';
    else
        pietools_root = '/Users/danilobraghini/Documents/GitHub/PIETOOLS/PIETOOLS';
    end
end
lr_root = fullfile(pietools_root, 'PIETOOLS_demos', 'lowrank_2d_stability');
cu_root = fullfile(pietools_root, 'PIETOOLS_demos', 'cuadmm');
if exist(fullfile(pietools_root, 'pietools_path_update.m'), 'file') ~= 2
    error('PIETOOLS root not found: %s', pietools_root);
end
if exist(fullfile(lr_root, 'build_stab_2d_st2.m'), 'file') ~= 2
    error('lowrank_2d_stability folder not found: %s', lr_root);
end
if exist(fullfile(cu_root, 'dump2cuadmm.m'), 'file') ~= 2
    error('cuADMM demo folder not found: %s', cu_root);
end

addpath(genpath(pietools_root));
addpath(lr_root);
addpath(cu_root);
sedumi_root = getenv('CUADMM_SEDUMI');
if ~isempty(sedumi_root) && exist(sedumi_root, 'dir')
    addpath(sedumi_root);
end
mosek_root = getenv('CUADMM_MOSEK');
if ~isempty(mosek_root) && exist(mosek_root, 'dir')
    addpath(mosek_root);
end

oldpwd = pwd;
cleanup = onCleanup(@() cd(oldpwd));
run_root = oldpwd;
if isempty(cuadmm_outdir)
    cuadmm_outdir = fullfile(run_root, 'heat2d_cuadmm_runs');
end
cuadmm_exe = pick_cuadmm_exe(cuadmm_exe);
cd(lr_root);
pielr_path();
%% PDE in PIE state-space 
lam_star = 2*pi^2;
fprintf('\n2D heat stability comparison\n');
fprintf('  PDE: x_t = x_s1s1 + x_s2s2 + lam*x on [0,1]^2, x in R^%d\n', n);
fprintf('  BC : x = 0 on all four edges\n');
fprintf('  lam = %.8g  (lam / (2*pi^2) = %.6g)\n', lam, lam/lam_star);
fprintf('  analytic stable iff lam < %.8g\n', lam_star);
fprintf('  LPI: build_stab_2d_st2 with Dup=%d\n', nd);
t0 = tic;
pvar s1 s2
clear stateNameGenerator
x = pde_var('state', n, [s1;s2], [0,1;0,1]);
sys = [diff(x,'t') == diff(x,s1,2) + diff(x,s2,2) + lam*x;
       subs(x,s1,0) == zeros(n,1);
       subs(x,s1,1) == zeros(n,1);
       subs(x,s2,0) == zeros(n,1);
       subs(x,s2,1) == zeros(n,1)];
PIE = convert(sys, 'pie');
t_pie = toc(t0);
%% SDP settings
st = settings_PIETOOLS_light_2D();
if ~isempty(dmult)
    st.LF_deg.d2{1,1} = dmult*ones(2,2);
end
dx = st.LF_deg.dx;
dy = st.LF_deg.dy;
d2 = st.LF_deg.d2;
st.eq_deg.dx = {nd+dx{1}; nd+dx{2}; nd+dx{3}};
st.eq_deg.dy = {nd+dy{1}, nd+dy{2}, nd+dy{3}};
st.eq_deg.d2 = cell(3,3);
for a = 1:3
    for b = 1:3
        st.eq_deg.d2{a,b} = nd+d2{a,b};
    end
end
st.LF_use_psatz = 0;
st.eq_use_psatz = [0;0];
st.eppos = 1e-2*ones(4,1);
st.epneg = 0;
st.use_sosineq = 0;
if ~isfield(st,'eq_opts') || ~isfield(st.eq_opts,'exclude')
    st.eq_opts = struct('psatz',0,'sep',zeros(1,6),'exclude',zeros(1,16));
end
if ~isfield(st,'eq_deg_psatz'),  st.eq_deg_psatz  = {st.eq_deg};  end
if ~isfield(st,'eq_opts_psatz'), st.eq_opts_psatz = {st.eq_opts}; end
%% Use low rank files created by CC
A = pielr_private('pielr_lpi', 'stability');
A.gate = struct('abs', gate, 'k', gate_k, ...
    'ref', NaN, 'refmax', gate_refmax);

PIEn = pielr_private('pielr_norm_pie', PIE);
t0 = tic;
[prog, H] = A.build(PIEn, st);
D = pielr_private('pielr_rawdata', prog);
t_setup = toc(t0);

shape = struct();
shape.m = numel(D.b);
shape.Kf = D.L.Kf;
shape.Ns = D.L.N;
shape.Ntot = D.L.Ntot;
shape.unknowns_full = sum(D.L.N .* (D.L.N + 1) / 2);
shape.nnzAt = nnz(D.At);

fprintf('\nAssembled LPI\n');
fprintf('  PIE build: %.2f s\n', t_pie);
fprintf('  LPI setup: %.2f s\n', t_setup);
fprintf('  m constraints: %d\n', shape.m);
fprintf('  free variables Kf: %d\n', shape.Kf);
fprintf('  Gram blocks Ns: [%s]\n', num2str(shape.Ns));
fprintf('  full Gram unknowns: %d\n', shape.unknowns_full);
fprintf('  nnz(At): %d\n', shape.nnzAt);

out = struct();
out.params = struct('n', n, 'lam', lam, 'Dup', nd, 'dmult', dmult, ...
    'run_mosek', run_mosek, 'run_lowrank', run_lowrank, ...
    'gate', gate, 'gate_k', gate_k, 'gate_refmax', gate_refmax, ...
    'maxrank', maxrank, 'rank', rank, 'seeds', seeds, 'lmit', lmit, ...
    'lmtol', lmtol, 'pietools_root', pietools_root, 'solver', solver, ...
    'cuadmm_exe', cuadmm_exe, 'cuadmm_outdir', cuadmm_outdir, ...
    'cuadmm_tol', cuadmm_tol, 'cuadmm_maxiter', cuadmm_maxiter, ...
    'cuadmm_timeout', cuadmm_timeout, 'cuadmm_max_wall', cuadmm_max_wall, ...
    'cuadmm_tag', cuadmm_tag, 'cuadmm_pinf_norm', cuadmm_pinf_norm, ...
    'cuadmm_require_full_gpu', cuadmm_require_full_gpu, ...
    'cuadmm_launcher', cuadmm_launcher, 'verbose', verbose);
out.lam_star = lam_star;
out.PIE = PIE;
out.settings = st;
out.prog = prog;
out.H = H;
out.D = D;
out.shape = shape;
out.t_pie = t_pie;
out.t_setup = t_setup;
out.mosek = [];
out.lowrank = [];
out.cuadmm = struct('ran', false, 'skip_reason', '');
%% Output reference interior-point solution
%ref.ok     % did the IPM solution pass the certificate gate?
%ref.rel    % main normalized operator residual, used for acceptance
%ref.rel_d  % older diagnostic residual with older denominator
%ref.t      % time in seconds for the IPM reference solve
ref_rel = NaN;
if run_mosek
    fprintf('\nMosek/reference solve on the exact same assembled LPI...\n');
    ref = pielr_private('pielr_ipm_ref', prog, H, D, A, verbose, solver);
    out.mosek = ref;
    if isempty(ref.err)
        ref_rel = ref.rel;
        fprintf('  solver: %s\n', ref.solver);
        fprintf('  ok: %d, rel: %.3e, rel_d: %.3e, time: %.2f s\n', ...
            ref.ok, ref.rel, ref.rel_d, ref.t);
    else
        fprintf('  reference solve failed: %s\n', ref.err);
    end
else
    fprintf('\nMosek/reference solve skipped.\n');
end
if run_lowrank
    fprintf('\nLow-rank solve...\n');
    lr_opts = struct();
    lr_opts.settings = st;
    lr_opts.gate = gate;
    lr_opts.gate_k = gate_k;
    lr_opts.gate_refmax = gate_refmax;
    lr_opts.maxrank = maxrank;
    lr_opts.rank = rank;
    lr_opts.seeds = seeds;
    lr_opts.lmit = lmit;
    lr_opts.lmtol = lmtol;
    lr_opts.verbose = verbose;
    if isfinite(ref_rel) && ref_rel > 0
        lr_opts.ref_rel = ref_rel;
    end
    cert = pielr_solve(PIE, 'stability', lr_opts);
    out.lowrank = cert;
    fprintf('  ok: %d, rel: %.3e, thresh: %.3e, score: %.3g\n', ...
        cert.ok, getf(cert,'rel'), getf(cert,'thresh'), getf(cert,'score'));
    fprintf('  setup: %.2f s, discovery: %.2f s\n', ...
        getf(cert,'t_setup'), getf(cert,'t_discover'));
    if isfield(cert, 'r')
        fprintf('accepted face widths r: [%s]\n', num2str(cert.r));
    end
else
    fprintf('\nLow-rank solve skipped.\n');
end

fprintf('\ncuADMM export and GPU solve for the same assembled LPI...\n');
out.cuadmm = export_cuadmm_problem(D, shape, n, lam, nd, dmult, cuadmm_outdir);
fprintf('  dump: %s\n', out.cuadmm.dumpfile);
fprintf('  text problem dir: %s\n', out.cuadmm.problem_dir);
fprintf('  export time: %.2f s\n', out.cuadmm.t_export);

if isempty(cuadmm_exe) || exist(cuadmm_exe, 'file') ~= 2
    error('cuADMM executable not found: %s', cuadmm_exe);
end
fprintf('  executable: %s\n', cuadmm_exe);
cu_opts = struct();
cu_opts.mode = 'feas';
cu_opts.solver = 'cuadmm';
cu_opts.exe = cuadmm_exe;
cu_opts.probe_tol = cuadmm_tol;
cu_opts.probe_cap = cuadmm_maxiter;
cu_opts.run_timeout = cuadmm_timeout;
cu_opts.max_wall = cuadmm_max_wall;
cu_opts.outdir = cuadmm_outdir;
cu_opts.tag = cuadmm_tag;
cu_opts.pinf_norm = cuadmm_pinf_norm;
cu_opts.require_full_gpu = cuadmm_require_full_gpu;
cu_opts.launcher = cuadmm_launcher;

t0 = tic;
R = bl_bisect(out.cuadmm.dumpfile, cu_opts);
out.cuadmm.wall = toc(t0);
out.cuadmm.R = R;
out.cuadmm.ran = true;
out.cuadmm.run_dir = fullfile(cuadmm_outdir, sprintf('%s_%s', R.id, R.solver));
out.cuadmm.solution_file = fullfile(out.cuadmm.run_dir, 'run_001_X.txt');
out.cuadmm.certificate_file = fullfile(out.cuadmm.run_dir, 'run_001_Frep.txt');
pr = R.probes(end);
fprintf('  verdict: %s, iters: %.0f, solver: %.2f s, wall: %.2f s\n', ...
    pr.verdict, pr.iters, R.t_solver, out.cuadmm.wall);
fprintf('  pinf_min: %.3e, pinf_end: %.3e, F_eta_psd: %.3e\n', ...
    pr.pinf_min, pr.pinf_end, getf(pr,'F_etapsd'));

fprintf('\nComparison summary\n');
fprintf('  lam/(2*pi^2): %.6g\n', lam/lam_star);
if ~isempty(out.mosek)
    fprintf('  Mosek:    ok=%d rel=%.3e t=%.2f s\n', ...
        out.mosek.ok, out.mosek.rel, out.mosek.t);
end
if ~isempty(out.lowrank)
    fprintf('  Low-rank: ok=%d rel=%.3e t=%.2f s\n', ...
        out.lowrank.ok, getf(out.lowrank,'rel'), ...
        getf(out.lowrank,'t_setup') + getf(out.lowrank,'t_discover'));
end
if isfield(out, 'cuadmm') && isstruct(out.cuadmm) && out.cuadmm.ran
    pr = out.cuadmm.R.probes(end);
    fprintf('  cuADMM:   verdict=%s pinf_min=%.3e it=%d solver=%.2f s wall=%.2f s\n', ...
        pr.verdict, pr.pinf_min, out.cuadmm.R.iters, ...
        out.cuadmm.R.t_solver, out.cuadmm.wall);
elseif isfield(out, 'cuadmm') && isstruct(out.cuadmm) && ...
        isfield(out.cuadmm, 'dumpfile') && ~isempty(out.cuadmm.dumpfile)
    fprintf('cuADMM:   exported, not solved (%s)\n', out.cuadmm.skip_reason);
end

function exe = pick_cuadmm_exe(requested)
if ~isempty(requested)
    exe = requested;
    return
end
env_exe = getenv('CUADMM_EXE');
if ~isempty(env_exe)
    exe = env_exe;
    return
end
candidates = { ...
    '/scratch/mpeet/pietools/solvers/cuADMM/build_sol/cuadmm_exe', ...
    '/scratch/mpeet/pietools/solvers/cuADMM/build/cuadmm_exe', ...
    '/home/mpeet/solvers/cuADMM/build/cuadmm_exe', ...
    '/Users/danilobraghini/Documents/GitHub/cuADMM/build/cuadmm_exe'};
exe = '';
for k = 1:numel(candidates)
    if exist(candidates{k}, 'file') == 2
        exe = candidates{k};
        return
    end
end
end


function cu = export_cuadmm_problem(D, shape, n, lam, Dup, dmult, cuadmm_outdir)
if ~exist(cuadmm_outdir, 'dir')
    mkdir(cuadmm_outdir);
end

if isfield(D, 'c') && nnz(D.c) > 0 && norm(full(D.c(:))) > 1e-12
    error('cuADMM feasibility export expects c = 0; norm(D.c) = %.3e', norm(full(D.c(:))));
end

case_id = sprintf('heat2d_n%d_lam%s_Dup%d', n, number_tag(lam), Dup);
dumpfile = fullfile(cuadmm_outdir, [case_id '.mat']);
metafile = fullfile(cuadmm_outdir, [case_id '_meta.mat']);

At = D.At;
bscl = norm(full(D.b(:)));
if bscl == 0
    bscl = 1;
end
b = sparse(D.b(:) / bscl);
c = sparse(size(At, 1), 1);
K = struct('f', D.L.Kf, 'l', 0, 'q', [], 's', D.L.N);
Ns = D.L.N;
Kf = D.L.Kf;
save(dumpfile, 'At', 'b', 'c', 'K', 'Ns', 'Kf', '-v7.3');

meta = struct();
meta.bscl = bscl;
meta.shape = shape;
meta.n = n;
meta.lam = lam;
meta.Dup = Dup;
meta.dmult = dmult;
save(metafile, '-struct', 'meta');

problem_dir = fullfile(cuadmm_outdir, [case_id '_cuadmm_text']);
t0 = tic;
info = dump2cuadmm(dumpfile, problem_dir);
t_export = toc(t0);

cu = struct();
cu.ran = false;
cu.skip_reason = '';
cu.dumpfile = dumpfile;
cu.metafile = metafile;
cu.problem_dir = problem_dir;
cu.t_export = t_export;
cu.info = info;
cu.bscl = bscl;
cu.solution_file = '';
cu.certificate_file = '';
end

function s = number_tag(x)
s = sprintf('%.8g', x);
s = strrep(s, '-', 'm');
s = strrep(s, '+', '');
s = strrep(s, '.', 'p');
s = strrep(s, 'e', 'em');
end

function v = getf(s, name)
if isfield(s, name)
    v = s.(name);
else
    v = NaN;
end
end
