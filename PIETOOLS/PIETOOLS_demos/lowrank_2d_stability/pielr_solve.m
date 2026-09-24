function cert = pielr_solve(PIE,lpi,opts)                                   % CC, 09/23/2026
% PIELR_SOLVE  Low-rank certification of a PIETOOLS LPI, for any adapted
% executive, in either spatial dimension.
%
%   cert = pielr_solve(PIE,'stability')
%   cert = pielr_solve(PIE,'l2gain',opts)
%
% The signature deliberately mirrors lpiscript(PIE,lpi,opts): same first two
% arguments, same LPI vocabulary, same settings argument (a struct, or a
% preset name that lpisettings understands).  Where lpiscript SOLVES the LPI
% with an interior-point method, this searches for a LOW-RANK certificate and
% verifies it at the operator level.
%
% TWO DIFFERENCES FROM lpiscript, both deliberate and both worth knowing:
%
%  1. 'stability' HERE MEANS THE DIRECT FORM.  lpiscript('stability') calls
%     PIETOOLS_PIE2PDEstability, the Q form.  This routine builds the LPI of
%     PIETOOLS_PDEstability / PIETOOLS_stability_2D -- a different and
%     stronger LPI, and the one this package's existing 1-D and 2-D drivers
%     have always built.  Named accordingly rather than silently.
%
%  2. lpiscript HAS NO 2-D DISPATCH AT ALL; every case routes to a 1-D
%     executive.  Here the dimension is taken from the PIE and the adapter
%     picks the matching transcription.
%
% WHAT IS CERTIFIED, AND HOW FAR IT GOES.  Per positivity block i the routine
% looks for a face V_i with orthonormal columns and a small S_i >= 0 with
% Gram block X_i = V_i S_i V_i'.  Restriction to a face can only LOSE feasible
% points, so any S_i >= 0 meeting the restricted system gives a feasible point
% of the original LPI: for a feasibility LPI that is a genuine certificate, and
% for an optimisation LPI (l2gain) a genuine UPPER bound on gamma.  A rank that
% fails is reported "not reached", never as a proven floor.
%
% ACCEPTANCE is pielr_opcheck: the operator residual of the executive's own
% constraint, measured against a normaliser the adapter names and justifies,
% AND every Gram block PSD.  See pielr_opcheck's header for why the
% normaliser is not the operator that must vanish.
%
% INPUT
%   PIE   pie_struct, pde_struct, or a struct carrying at least T and A
%   lpi   'stability' | 'stability-dual' | 'l2gain' | 'l2gain-dual'
%   opts  struct of options, all optional:
%     .settings  LPI settings struct, or a preset name for lpisettings.
%                Default 'light' in 1-D, set2d_deg(4,[]) in 2-D.
%     .rank      per-block rank (or scalar) to START the ladder at
%     .maxrank   ladder ceiling (default 6)
%     .gate      acceptance threshold on the operator residual (default 1e-6)
%     .seeds     BM random seeds per rank (default [11 22 33]).  NOTE the seed
%                INDEX also selects the initial magnitude, so the seeds are a
%                designed sweep over scale, not an i.i.d. ensemble.
%     .lmit      LM iterations per attempt (default 400, the shipped value)
%     .lmtol     LM residual target (default opts.gate)
%     .verbose   print progress (default true)
%
% OUTPUT (struct)
%   cert.ok        true iff a certificate passed the gate
%   cert.lpi       the LPI name;  cert.dim  1 or 2
%   cert.r         per-block face widths of the accepted certificate
%   cert.rank      per-block numerical rank of the lifted Gram blocks
%   cert.rel       operator residual against the adapter's normaliser  <- gate
%   cert.rel_d     the same residual against the operator that must vanish,
%                  i.e. what opcheck_2d used to report.  Kept so older 2-D
%                  numbers stay comparable and so a collapsing denominator is
%                  visible instead of silent.
%   cert.maxRes .maxNrm .maxDop   the maxima those ratios are built from.
%                  opcheck_2d computed these on every call and returned only
%                  the ratio, which is why no existing log can say whether two
%                  op_rel values differ in numerator or denominator.
%   cert.gam       gamma, for the optimisation LPIs (NaN otherwise)
%   cert.mineig .normQ  per Gram block, original units
%   cert.face .S   accepted face and its coefficients
%   cert.lay       the SeDuMi coordinate map actually used
%   cert.why       per-attempt log: rank, start, raw_rel, rel, exit reason,
%                  face np/rM/rel_eq_full.  bm_lm2 returns its exit reason and
%                  restrict_solve returns np, rM and the full-row residual;
%                  no existing caller captured either, which is why the
%                  package's own open questions about which exit fired and
%                  what rM is could not be answered from its logs.
%   cert.notes     fine print accumulated during the run
%
% Requires PIETOOLS and SeDuMi on the path: run pielr_path once per session.

if nargin<2 || isempty(lpi), lpi = 'stability'; end
if nargin<3, opts = struct(); end
def = struct('settings',[],'rank',[],'maxrank',6,'gate',1e-6, ...
             'seeds',[11 22 33],'lmit',400,'lmtol',[],'verbose',true, ...
             'refine',true,'w0',[],'ref_rel',NaN,'gate_k',10,'gate_refmax',[]); % CC, 09/23/2026
fn = fieldnames(def);
for k = 1:numel(fn)
    if ~isfield(opts,fn{k}) || isempty(opts.(fn{k})), opts.(fn{k}) = def.(fn{k}); end
end
% lmtol DEFAULTS TO gate/100, NOT gate (CC, 09/24/2026).  bm_lm2 compares tol
% against f = ||W r||, the WHITENED residual, while the gate is an operator
% residual; the two are not in the same units and their ratio is not 1.  Set
% equal to the gate, the 'tol' exit stops the search while the operator
% residual is still above threshold, and no amount of budget helps because
% the run never reaches maxit.
%
% MEASURED, tier 2, the three 1-D stability cases that missed at EVERY budget
% (8000, 32000 and 128000 all gave the same answer in the same 17 s, exiting
% tol:13 / tol:11 / tol:11):
%   case          lmtol 1e-6            lmtol 1e-8
%   rd1d-lam0.5   no  (tol:13)          yes rel 4.9327e-08  rank [1 1 1]
%   rd1d-lam0.9   no  (tol:11)          yes rel 2.9318e-08  rank [2 2 2]
%   rd1d-n2       no  (tol:11)          yes rel 2.6312e-08  rank [2 2 2]
% The residuals reached are 20-40x BELOW the gate, so the certificates were
% always in reach and the exit was stopping two orders short.  A factor 100
% is chosen to sit well inside that margin; the right long-term fix is for
% the exit to test the quantity the gate actually uses.
if isempty(opts.lmtol), opts.lmtol = opts.gate/100; end                    % CC, 09/24/2026
vb = opts.verbose;

A = pielr_lpi(lpi);
% The acceptance threshold is relative to a REFERENCE residual on the same
% program when one is supplied (opts.ref_rel, e.g. from pielr_ipm_ref):
%   accept iff  rel <= max(opts.gate, opts.gate_k * opts.ref_rel)
% and the reported score is rel/ref_rel.  Measured on the historical 2-D
% benchmark, the interior-point solution itself only reaches rel = 3.497e-06,
% so a fixed 1e-6 there demands more precision than the program admits from
% ANY method.  With no reference the behaviour is the old absolute gate.  See
% pielr_opcheck's header for why the absolute floor is kept alongside the
% ratio rather than replaced by it.
A.gate = struct('abs',opts.gate,'k',opts.gate_k,'ref',opts.ref_rel);   % CC, 09/23/2026
% refmax is passed ONLY when the caller set it, so pielr_opcheck's own default
% is the single definition.  It was duplicated here at first and the two
% drifted within the hour: opcheck's default was tightened 1e-4 -> 1e-5 while
% this line kept passing 1e-4, which silently won, and the negative control it
% was meant to reject went on certifying at rel = 1.020e-04.
if ~isempty(opts.gate_refmax)                                          % CC, 09/23/2026
    A.gate.refmax = opts.gate_refmax;                                  % CC, 09/23/2026
end                                                                    % CC, 09/23/2026

% ---- normalise the PIE, and settle the dimension ------------------------
PIE = pielr_norm_pie(PIE);
dim = pielr_pie_dim(PIE);
if isempty(opts.settings)
    if dim==2, opts.settings = set2d_deg(4,[]); else, opts.settings = lpisettings('light'); end
elseif ischar(opts.settings) || isstring(opts.settings)
    opts.settings = lpisettings(char(opts.settings));
end

notes = {};
cert = struct('ok',false,'lpi',A.name,'dim',dim,'r',NaN,'rank',NaN, ...
              'rel',NaN,'rel_d',NaN,'gam',NaN);

% ---- assemble ------------------------------------------------------------
if vb, fprintf('pielr_solve: assembling the %d-D ''%s'' LPI ...\n',dim,A.name); end
tS = tic;
[prog,H] = A.build(PIE,opts.settings);
D = pielr_rawdata(prog);
L = D.L;
t_setup = toc(tS);
if vb
    fprintf('  Gram blocks N = [%s]  free = %d  m = %d  unknowns = %d  [setup %.1fs]\n', ...
        num2str(L.N),L.Kf,numel(D.b),sum(L.N.*(L.N+1)/2),t_setup);
    if ~isempty(D.hineq)
        fprintf('  %d inequality row(s) held out of the equality system\n',numel(D.hineq));
    end
    if D.has_obj
        fprintf('  objective present (%d nonzero coefficient(s))\n',nnz(D.c));
    end
end

% ---- discovery -----------------------------------------------------------
% An LPI with an objective is certified by BISECTION on the objective
% variable: nothing in the search optimises anything (bm_lm2 minimises the
% equality residual, in which c does not appear, and restrict_solve's
% determined solve has no degree of freedom left), so run straight it returns
% whatever gamma the face happens to carry.  MEASURED on
% Ex_Transport_Eq_with_Disturbance at 'light': 1.42095 unpinned against
% 0.516333 from the interior-point solve of the same program.  See
% pielr_bisect_obj, including the cost objection.
tD = tic;
if D.has_obj && ~(isfield(opts,'no_bisect') && opts.no_bisect)
    [V,R,q,why,dnotes,P] = pielr_bisect_obj(prog,H,D,A,opts,vb);
    t_disc = toc(tD);
    notes = [notes dnotes];
    cert = finish(cert,notes,L,t_setup,t_disc,why,V,R,q,P,opts,vb,prog,H,D,A);
    return
end

% ---- package for Burer-Monteiro -----------------------------------------
P = bm_setup(D.At,D.b,L.N,L.Kf,1,L);      % pre=1: BM needs the whitened metric
if P.bout > 1e-10
    % b outside the numerical range of A: the whitened residual is blind to
    % that component by construction, so ||F|| -> 0 is attainable while the
    % raw residual keeps an irreducible floor.  bm_setup has always computed
    % this and nothing has ever read it.
    notes{end+1} = sprintf( ...
        'b lies %.3e outside range(A); the whitened objective cannot see that component', ...
        P.bout);
    if vb, fprintf('  WARNING: %s\n',notes{end}); end
end

% ---- discover a face, certify it ----------------------------------------
[V,R,q,why,dnotes] = pielr_discover(prog,H,D,P,A,opts,vb);
t_disc = toc(tD);
notes = [notes dnotes];
cert = finish(cert,notes,L,t_setup,t_disc,why,V,R,q,P,opts,vb,prog,H,D,A);
end

% =========================================================================
function cert = finish(cert,notes,L,t_setup,t_disc,why,V,R,q,P,opts,vb,prog,H,D,A)
% Fill and report.  Shared by the feasibility and the bisection paths so the
% two cannot drift in what they record.
% Shrink the accepted face before reporting: the rank the ladder stopped at
% is the first rung that worked, not the narrowest face that holds, and the
% rank is the headline number of every low-rank claim.
if opts.refine && ~isempty(R) && R.ok
    [V,R,q,rnotes] = pielr_refine(prog,H,D,P,A,V,R,q,vb);
    notes = [notes rnotes];
end
% A bisection repeats the same note once per certifying trial; the reader
% wants to know the fact, not the count.
if ~isempty(notes), notes = unique(notes,'stable'); end
cert.lay = L;   cert.Ns = L.N;   cert.Kf = L.Kf;
cert.t_setup = t_setup;   cert.t_discover = t_disc;
cert.why = why;
cert.unknowns_full = sum(L.N.*(L.N+1)/2);

% EVERY FIELD IS SET ON BOTH PATHS.  The failure path used to return before
% score/thresh/ref_rel existed, so a caller reading cert.score on a run that
% did not certify got "Unrecognized field name" rather than NaN -- which is
% how a tier sweep over 18 problems silently lost the five cases that failed,
% i.e. exactly the cases the sweep existed to characterise.
cert.score   = NaN;   cert.thresh  = NaN;   cert.rel_i    = NaN;
cert.ref_rel = opts.ref_rel;  cert.gate_k = opts.gate_k;  cert.gate_abs = opts.gate;
cert.maxRes  = NaN;   cert.maxNrm  = NaN;   cert.maxDop   = NaN;
cert.maxRes_i = NaN;  cert.maxNrm_i = NaN;  cert.worst_eq = NaN;
cert.mineig  = NaN;   cert.normQ   = NaN;   cert.psd      = false;
cert.nrm_why = '';    cert.face    = {};    cert.S        = {};
cert.unknowns_face = NaN;
if isempty(R)
    cert.notes = [notes {sprintf(['no certificate at rank <= %d over %d seed(s)/rank; ' ...
        'reported as NOT REACHED, never as a proven floor'],opts.maxrank,numel(opts.seeds))}];
    if vb, pielr_report(cert); end
    return
end

cert.ok      = R.ok;
cert.rel     = R.rel;       cert.rel_d  = R.rel_d;
cert.score   = R.score;     cert.thresh = R.thresh;   % CC, 09/23/2026
cert.rel_i   = R.rel_i;     cert.maxRes_i = R.maxRes_i;  % per-equality  % CC, 09/23/2026
cert.maxNrm_i = R.maxNrm_i; cert.worst_eq = R.worst;                  % CC, 09/23/2026
cert.ref_rel = R.gate_ref;  cert.gate_k = R.gate_k;                    % CC, 09/23/2026
cert.gate_abs = R.gate_abs;                                  % CC, 09/23/2026
cert.maxRes  = R.maxRes;    cert.maxNrm = R.maxNrm;   cert.maxDop = R.maxDop;
cert.nrm_why = R.nrm_why;
cert.mineig  = R.mineig;    cert.normQ  = R.normQ;    cert.rank = R.rank;
cert.psd     = R.psd;
if isfield(R,'aux') && isfield(R.aux,'gam'), cert.gam = R.aux.gam; end
cert.face = V;
cert.r    = cellfun(@(v)size(v,2),V);
cert.unknowns_face = sum(cert.r.*(cert.r+1)/2);
cert.S    = pielr_face_coeffs(V,q,P);
cert.notes = notes;
if vb, pielr_report(cert); end
end
