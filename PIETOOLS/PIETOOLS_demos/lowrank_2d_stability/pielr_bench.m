function R = pielr_bench(opts)                                             % CC, 09/23/2026
% PIELR_BENCH  Run the low-rank certifier over the benchmark catalogue and
% score it against an interior-point reference measured in the same session.
%
%   R = pielr_bench                       % tier 1
%   R = pielr_bench(struct('tier',2))
%   R = pielr_bench(struct('lmit',4000,'tag','lmit4000'))
%
% WHAT THIS IS FOR.  To make a change to the solver MEASURABLE.  Run it once
% to bank a baseline, change one thing, run it again, and diff the two tables.
% Every number in the table is produced in the run that prints it; nothing is
% quoted from a previous session, a log file or a library comment.
%
% HOW THE COMPARISON IS MADE FAIR.  For each case the identical assembled
% program goes to both solvers, and BOTH answers are judged by the same gate
% (pielr_opcheck).  A solver's own feasibility flag never appears.  For the
% LPIs with an objective the score is the relative gap in gamma; for the
% feasibility LPIs it is whether each side's point passes the gate, plus the
% rank the low-rank side needed.
%
% WHAT THE COLUMNS MEAN
%   ipm_ok   the interior-point point passes the gate
%   lr_ok    the low-rank certificate passes the gate
%   rel      the low-rank operator residual (the gate quantity)
%   ipm_rel  the same quantity at the interior-point point: the residual a
%            full solve achieves on this program, and the reference the
%            threshold is now calibrated against
%   rel/ipm  THE SCORE.  1 would mean the low-rank certificate is as accurate
%            as the full solve on the same program.  This is the quality
%            measure; acceptance is rel <= max(opts.gate, opts.gate_k*ipm_rel)
%   r        per-block face widths
%   gam/ipm  gamma ratio, for the objective LPIs;  1.000 would be tight
%   t_lr     low-rank seconds (setup + discovery)
%   t_ipm    interior-point seconds
%
% A case where ipm_ok is FALSE is not a failure of the low-rank method: the
% LPI is simply infeasible at those settings, and the row is marked so.
%
% INPUT  opts  .tier (1..3, default 1)  .cases (name or cellstr filter)
%              .verbose (default false: one line per case)
%              .skip_ipm (default false)
%              .tag (a label stored in R for diffing two runs)
%              plus any pielr_solve option (lmit, maxrank, seeds, gate, ...)
% OUTPUT R     struct array, one entry per case, plus R(1).meta

if nargin<1, opts = struct(); end
def = struct('tier',1,'cases',[],'verbose',false,'skip_ipm',false,'tag','', ...
             'eppos',1e-2);                                            % CC, 09/23/2026
fn = fieldnames(def);
for k=1:numel(fn)
    if ~isfield(opts,fn{k}) || isempty(opts.(fn{k})), opts.(fn{k}) = def.(fn{k}); end
end
C = pielr_bench_cases(opts.tier);
if ~isempty(opts.cases)
    keep = false(1,numel(C));
    want = cellstr(opts.cases);
    for k=1:numel(C), keep(k) = any(strcmp(C(k).name,want)) || ...
            any(cellfun(@(w)contains(C(k).name,w),want)); end
    C = C(keep);
end
if isempty(C), error('pielr_bench:nocases','no cases matched'); end

solve_opts = rmfield_if(opts,{'tier','cases','skip_ipm','tag'});
solve_opts.verbose = opts.verbose;

fprintf('\npielr_bench: %d case(s), tier %d%s\n',numel(C),opts.tier, ...
    tern(isempty(opts.tag),'',['   tag ''' opts.tag '''']));
fprintf('matlab %s   threads %d   %s\n',version,maxNumCompThreads, ...
    char(datetime('now','Format','yyyy-MM-dd HH:mm:ss')));
hdr();

R = struct([]);
for k = 1:numel(C)
    c = C(k);
    r = struct('name',c.name,'family',c.family,'lpi',c.lpi,'tier',c.tier, ...
        'ok',false,'ipm_ok',false,'rel',NaN,'ipm_rel',NaN,'r',NaN, ...
        'gam',NaN,'ipm_gam',NaN,'ratio',NaN,'t_lr',NaN,'t_ipm',NaN, ...
        'score',NaN,'thresh',NaN, ...
        'N',[],'m',NaN,'Kf',NaN,'err','','note',c.note,'cert',[]);
    try
        % ONE evalc, here, rather than inside each loader.  The loaders are
        % invoked through function handles stored in a struct array, and an
        % evalc inside a function reached that way crashed MATLAB R2025b with
        % an access violation in its JIT, reproducibly, on the SECOND case of
        % any run (two separate case lists, two different second cases; the
        % crash trace names inEvalCmdWithLocalReturn).  A single evalc at a
        % fixed call site suppresses the same conversion output and does not.
        tmp = evalc('PIE = c.load();');   %#ok<NASGU>
        so = solve_opts;  so.settings = set_eppos(c,opts.eppos);

        % CC, 09/23/2026: the REFERENCE IS SOLVED FIRST, and its residual is
        % handed to pielr_solve as opts.ref_rel, because the acceptance
        % threshold is now relative to it: accept iff
        %   rel <= max(opts.gate, opts.gate_k * ref_rel).
        % A fixed 1e-6 is not a usable criterion across this catalogue -- the
        % interior-point residual spans 4.15e-10 to 9.18e-07 over the cases
        % that certify, and on the historical 2-D benchmark the
        % interior-point solution itself only reaches 3.497e-06, so 1e-6 there
        % asks the low-rank search to beat the full solve.
        if ~opts.skip_ipm
            A = pielr_private('pielr_lpi',c.lpi);
            PIEn = pielr_private('pielr_norm_pie',PIE);
            st = so.settings;   % the SAME settings the low-rank arm gets
            [prog,H] = A.build(PIEn,st);
            D = pielr_private('pielr_rawdata',prog);
            r.m = numel(D.b);
            ref = pielr_private('pielr_ipm_ref',prog,H,D,A,opts.verbose);
            r.ipm_ok = ref.ok;  r.ipm_rel = ref.rel;  r.t_ipm = ref.t;
            r.ipm_gam = ref.gam;
            if ~isempty(ref.err), r.err = ['ipm: ' ref.err]; end
            % Only a reference the gate itself accepts is worth calibrating
            % against: if the reference solve did not converge on this program
            % there is nothing to be relative to, and the absolute floor
            % stands.
            if isfinite(ref.rel) && ref.rel > 0
                so.ref_rel = ref.rel;
            end
        end

        t = tic;
        cert = pielr_solve(PIE,c.lpi,so);
        r.t_lr = toc(t);
        r.score  = getfd(cert,'score',NaN);
        r.thresh = getfd(cert,'thresh',NaN);
        r.ok = cert.ok;  r.rel = cert.rel;  r.N = cert.Ns;  r.Kf = cert.Kf;
        if cert.ok, r.r = cert.r; end
        if isfield(cert,'gam'), r.gam = cert.gam; end
        r.cert = cert;

        if ~isnan(r.gam) && ~isnan(r.ipm_gam) && r.ipm_gam>0
            r.ratio = r.gam/r.ipm_gam;
        end
    catch ME
        r.err = sprintf('%s @ %s:%d',ME.message,ME.stack(1).name,ME.stack(1).line);
    end
    row(r);
    if isempty(R), R = r; else, R(end+1) = r; end %#ok<AGROW>
end

R(1).meta = struct('tag',opts.tag,'tier',opts.tier,'opts',solve_opts, ...
    'matlab',version,'threads',maxNumCompThreads, ...
    'when',char(datetime('now','Format','yyyy-MM-dd HH:mm:ss')));
summary(R);
end

% =========================================================================
function hdr()
% 'ipm@1e-6', not 'ipm'.  The two columns answer DIFFERENT questions and must
% not be read as a head-to-head: the low-rank verdict is the relative gate
% (threshold derived from ipm_rel), while the reference is reported against the
% plain 1e-6, which asks "is 1e-6 reachable on this program at all".  Judging
% the reference by the relative gate would be vacuous -- its own residual is
% the reference, so it would pass at any k >= 1 by construction.  A row with
% lr = yes and ipm@1e-6 = no is therefore expected wherever 1e-6 is out of
% reach, and is NOT an anomaly; the comparison is the rel/ipm score.
fprintf('\n%-22s %-14s %-5s %-9s %-10s %-10s %-9s %-9s %-8s %-7s %s\n', ...
    'case','lpi','lr','ipm@1e-6','rel','ipm_rel','rel/ipm','r','gam/ipm','t_lr','t_ipm');
fprintf('%s\n',repmat('-',1,124));
end

function row(r)
g = '';
if ~isnan(r.ratio), g = sprintf('%.4f',r.ratio);
elseif ~isnan(r.gam), g = sprintf('%.4g',r.gam); end
s = '-';
if isfinite(r.score), s = sprintf('%.4g',r.score); end
fprintf('%-22s %-14s %-5s %-9s %-10.3e %-10.3e %-9s %-9s %-8s %-7.1f %.1f%s\n', ...
    r.name,r.lpi,yn(r.ok),yn(r.ipm_ok),r.rel,r.ipm_rel,s, ...
    tern(isnan(r.r(1)),'-',mat2str(r.r)),g,r.t_lr,r.t_ipm, ...
    tern(isempty(r.err),'',['   ! ' r.err]));
end

function summary(R)
n = numel(R);
ipm = [R.ipm_ok];  lr = [R.ok];
fprintf('%s\n',repmat('-',1,120));
sc = [R.score];  sc = sc(isfinite(sc) & [R.ok]);
if ~isempty(sc)
    fprintf(['score rel/ipm_rel over the %d accepted case(s): min %.4g  ' ...
             'median %.4g  max %.4g\n'],numel(sc),min(sc),median(sc),max(sc));
    fprintf(['  (1 would mean as accurate as the full solve; this is the ' ...
             'quality measure, the threshold is max(abs, k*ipm_rel))\n']);
end
fprintf('cases %d | reference reaches 1e-6 on %d | low-rank certifies %d', ...
    n,sum(ipm),sum(lr));
both = sum(ipm & lr);   missed = sum(ipm & ~lr);   extra = sum(~ipm & lr);
fprintf(' | agree %d | LOW-RANK MISSED %d',both,missed);
if extra>0
    fprintf('\n  %d case(s) certify where 1e-6 is out of reach for the reference too', extra);
    fprintf('\n  (expected: the threshold there is k*ipm_rel, not 1e-6 -- read the score)');
end
fprintf('\n');

% SOUNDNESS ASSERTION.  A face restriction can only LOSE feasible points, so
% the restricted optimum cannot lie below the full optimum: a gamma ratio under
% 1 is a contradiction, not tolerance, and means the gate is loose enough to be
% accepting points that do not satisfy the LPI.  MEASURED: at gate_k = 100 and
% 1000 this fires on gain-reacdiff (0.6143 and 0.1936), which is how k was
% bounded.  Any run that trips it is invalid, not merely surprising.
% Tolerance 1e-3, not 1e-6.  Both sides are approximate solutions of the same
% LPI, so at tight agreement the ratio crosses 1 on solver noise alone:
% MEASURED on Poincare tier 1 at lmit 8000, low-rank 0.426785 against the
% reference 0.4270904, a ratio of 0.99928 -- a 0.07% crossing that is not
% unsoundness.  The real violations this is for were 0.6143 and 0.1936.
gg = [R.ratio];  bad = find(isfinite(gg) & gg < 1-1e-3);
for i = bad
    fprintf(['*** UNSOUND: %s reports gamma %.4g, BELOW the reference %.4g ' ...
             '(ratio %.4f).\n    A face can only raise the optimum, so the ' ...
             'gate is accepting infeasible points.  Lower gate_k.\n'], ...
        R(i).name,R(i).gam,R(i).ipm_gam,R(i).ratio);
end
g = [R.ratio];  g = g(~isnan(g));
if ~isempty(g)
    fprintf('gamma ratio (low-rank / interior-point, 1.000 is tight): ');
    fprintf('min %.4f  median %.4f  max %.4f  over %d case(s)\n', ...
        min(g),median(g),max(g),numel(g));
end
e = {R.err};  e = e(~cellfun(@isempty,e));
for i=1:numel(e), fprintf('error: %s\n',e{i}); end
fprintf('\n');
end

function st = set_eppos(c,ep)                                              % CC, 09/23/2026
% Set the coercivity margin eppos for a STABILITY case.
%
% WHY THIS IS A BENCHMARK KNOB AND NOT A DETAIL.  For the stability LPI the
% negativity operator Dop = T'PA + A'PT is entirely LINEAR in P, so the only
% constant term in the equality Deop + Dop = 0 is the one carried by the
% eppos*I inside Pop.  The SDP's right-hand side is therefore PROPORTIONAL TO
% eppos, and the executives default it to 1e-4/1e-6.  That makes b tiny, the
% feasible Gram blocks tiny with it, and every error metric built on them
% confounded with the scale: MEASURED on the 2-D benchmark, max|Pop| = 1.35e-05
% against an eppos floor of 1e-6, i.e. the Lyapunov operator is barely more
% than its own epsilon.  Maintainer, 09/23/2026: the default "was foolishly
% set too low in stability tests, so all the solutions are small and confound
% the error metrics"; it should be nearer 1e-2, or 1.
%
% The l2gain LPI is NOT touched.  There Dop carries Dzwop and Czop as
% standalone constant blocks, so b is nonzero from the plant data whatever
% eppos is -- which is why PIETOOLS_Hinf_gain puts no eppos term on Rop at all,
% and why it can be exactly zero there.
st = c.settings;
if ischar(st) || isstring(st), st = lpisettings(char(st)); end
if isempty(ep) || ~startsWith(c.lpi,'stability'), return, end
if isfield(st,'eppos') && numel(st.eppos) == 4
    st.eppos = ep*ones(4,1);              % 2-D: [R, L2x, L2y, L2xy]
else
    st.eppos = ep;                        % 1-D: real-valued states
    st.eppos2 = ep;                       %      distributed states
end
if isfield(st,'settings_2d') && ~isempty(st.settings_2d)
    st.settings_2d.eppos = ep*ones(4,1);
end
end

function s = yn(t), if t, s='yes'; else, s='no'; end, end
function v = getfd(s,f,d)
if isstruct(s) && isfield(s,f) && ~isempty(s.(f)), v = s.(f); else, v = d; end
end
function s = tern(c,a,b), if c, s=a; else, s=b; end, end
function o = rmfield_if(o,f)
for i=1:numel(f), if isfield(o,f{i}), o = rmfield(o,f{i}); end, end
end
