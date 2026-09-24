function cert = pielr_certify(PIE,opts)
% PIELR_CERTIFY  Low-rank certification of 2-D PIE stability.
%
%   cert = pielr_certify(PIE)
%   cert = pielr_certify(PIE,opts)
%   cert = pielr_certify(PIE,prev_cert)   % a previous output IS a valid input:
%   cert = pielr_certify(PIE,faceN)       % certification only, no discovery
%
% Builds the direct-form 2-D stability LPI for the PIE (T,A) -- the same LPI as
% executives/2D/PIETOOLS_stability_2D, psatz 0 -- then finds a LOW-RANK
% certificate: per positivity block i, a face V_i (orthonormal columns) and a
% small S_i >= 0 with Gram block X_i = V_i*S_i*V_i'.  The certificate is
% verified AT THE OPERATOR LEVEL, never by an equality residual or a solver
% flag: on this LPI family norm(b) ~ 5e-6, so X = 0 satisfies the equality rows
% to ~1e-9 relative while certifying nothing, and a least-squares point was
% measured at op-residual 1.7e-11 with mineig < 0.  Acceptance is therefore
%     (operator residual over all 36 opvar2d cells < opts.gate)  AND
%     (every Gram block PSD).
% Every reported rank is an UPPER BOUND exhibited by a verified certificate; a
% rank that fails is reported as "not reached", never as a proven floor.
%
% INPUT
%   PIE   What convert() returns (pie_struct or plain struct) for a 2-D PDE.
%         Only PIE.T and PIE.A are used (autonomous stability).  You can also
%         hand in RAW OPERATORS without calling convert: assemble
%             PIE.T = Top;  PIE.A = Aop;          % opvar2d objects
%             PIE.vars = [Top.var1, Top.var2];    % optional: defaulted from T
%             PIE.dom  = Top.I;                   % optional: defaulted from T
%   opts  struct, all fields optional:
%     .settings  2D settings struct (default set2d_deg(4,[]))
%     .rank      per-block rank (or scalar) to START discovery at.  Default:
%                the ESTIMATE 2*n2 (n2 = number of 2-D PDE states), the
%                measured law on this LPI family; ranks below the start are
%                recovered by the refine shrink, never swept blind
%     .face      warm-start face: struct with .V (1xB cell, orthonormal cols)
%                and optionally .S (1xB cell, ORIGINAL-b units; see UNITS).
%                Skips discovery entirely.  Build one with pielr_tensor, load
%                the shipped face_n1.mat, or -- simplest -- pass a PREVIOUS
%                pielr_certify OUTPUT directly:
%                    opts.face = cert;                % cert from an earlier call
%                its cert.face / cert.S are picked up automatically.
%     .route     'bm' | 'mintrace' | 'auto' (default 'auto': Burer-Monteiro if
%                its helpers shipped, else min-trace SeDuMi + truncation)
%     .maxrank   discovery search ceiling (default 6, lifted to rank+2 when
%                the rank was estimated)
%     .gate      operator-residual acceptance gate (default 1e-6)
%     .refine    per-block face shrinking after acceptance (default: true when
%                discovery ran, false when opts.face was supplied)
%     .seeds     BM random seeds per rank (default [11 22 33]).  NOTE the seed
%                INDEX also picks the initial magnitude from [1e-1 1 1e1 1e2
%                1e4 1e6], so the 4th..6th scales are only reached by passing
%                four or more seeds.
%     .w0        {Y_1,...,Y_B} per-block factor to SEED the search with       % CC, 09/21/2026
%                (X_i ~ Y_i*Y_i'), e.g. the factor returned by an external     % CC, 09/21/2026
%                low-rank solver.  Tried FIRST at every rank rung, ahead of    % CC, 09/21/2026
%                the random seeds.  Unlike .face, which freezes a subspace     % CC, 09/21/2026
%                and skips discovery, this one only starts the descent -- the  % CC, 09/21/2026
%                columns stay free to move.  Columns beyond the current rank   % CC, 09/21/2026
%                are filled with small random values, never zeros (see fitw).  % CC, 09/21/2026
%     .lmit      LM iterations per attempt (default 400, the shipped value)  % CC, 09/21/2026
%                MEASURED 1-D, 22 lam x 5 degrees: 400 is a HARD cut, not a  % CC, 09/21/2026
%                convergence test -- bm_lm2 exits only when its damping loop  % CC, 09/21/2026
%                stalls or ||F|| < 1e-16 (never), so whether the residual     % CC, 09/21/2026
%                crosses the gate inside exactly 400 steps is near-arbitrary  % CC, 09/21/2026
%                in the operating point, which shows up as a RAGGED reach     % CC, 09/21/2026
%                (verdict flipping on 0.01 steps in lam).  At lmit = 4000 the % CC, 09/21/2026
%                1-D reach went 0.41 -> 0.999 of lam* at 'heavy' and the gap  % CC, 09/21/2026
%                to the best known relaxation closed to zero at every degree  % CC, 09/21/2026
%                tested -- and where 400 FAILS, 4000 is usually FASTER too,   % CC, 09/21/2026
%                since it accepts at low rank early instead of grinding the   % CC, 09/21/2026
%                rank ladder (heavy/0.50: 39 attempts/23.4 s fail -> 5        % CC, 09/21/2026
%                attempts/6.4 s certify).  Where 400 already works, 4000      % CC, 09/21/2026
%                costs ~4.5x, because there is no early exit.  The real fix   % CC, 09/21/2026
%                is a stagnation test in bm_lm2; this knob exposes the cut.   % CC, 09/21/2026
%     .verbose   print progress and the final report (default true)
%
% OUTPUT (struct)
%   cert.ok             true iff a certificate passed the gate
%   cert.r              per-block face widths of the accepted certificate
%   cert.rank           per-block NUMERICAL rank of the lifted Gram blocks
%   cert.op_rel         operator residual max|Qop+Deop| / max|Qop|
%   cert.mineig         per-block minimum eigenvalue (original units)
%   cert.face           {V_i} accepted face
%   cert.S              {S_i} face coefficients, ORIGINAL-b units
%   cert.Ns             Gram block sizes;  cert.Kf free variables
%   cert.unknowns_full  sum Ns(Ns+1)/2;  cert.unknowns_face  sum r(r+1)/2
%   cert.t_setup        seconds: PIE -> LPI assembly + SDP data extraction
%   cert.t_discover     seconds: face search (0 when opts.face given)
%   cert.t_certify      seconds: accepted solve/lift + verification + refine
%   cert.route          'face-replicate' | 'face-sdp' | 'bm' | 'mintrace'
%   cert.part           per-block group partition (for pielr_tensor), or {}
%   cert.notes          honest fine print accumulated during the run
%
% UNITS.  The SDP data is solved with b normalised to unit norm; opcheck takes
% Gram vectors in NORMALISED-b units and converts internally.  Everything this
% function RETURNS (cert.S, cert.mineig) is in ORIGINAL units: X_i = V_i S_i V_i'
% is the Gram block of the unnormalised program.  The conversion is owned here
% -- callers never touch nb0.  Getting this wrong reports op_rel = 1 for a
% perfectly good certificate.
%
% Requires PIETOOLS and SeDuMi on the path: run pielr_path once per session.
%
% CC, 09/19/2026: a previous pielr_certify output (or a bare face struct) is
%                 now accepted BOTH as opts.face AND directly as the second
%                 argument, so certificates chain between calls with no
%                 intermediate struct: cert2 = pielr_certify(PIE2,cert).

% ---- options ---------------------------------------------------------------
if nargin<2, opts = struct(); end
% The second argument may be a PREVIOUS CERTIFICATE or a bare FACE instead of  % CC, 09/19/2026
% an options struct: cert2 = pielr_certify(PIE2,cert) or (PIE2,faceN).  A     % CC, 09/19/2026
% cert is recognised by its ok+face fields, a face by a .V cell; neither      % CC, 09/19/2026
% collides with any option name.  Both wrap into opts.face and skip discovery. % CC, 09/19/2026
if isstruct(opts) && isfield(opts,'ok') && isfield(opts,'face')               % CC, 09/19/2026
    opts = struct('face',opts);                                               % CC, 09/19/2026
elseif isstruct(opts) && isfield(opts,'V') && iscell(opts.V) ...              % CC, 09/19/2026
        && ~isfield(opts,'face')                                              % CC, 09/19/2026
    opts = struct('face',opts);                                               % CC, 09/19/2026
end                                                                           % CC, 09/19/2026
def = struct('settings',[],'rank',[],'face',[],'route','auto','maxrank',6, ...
             'gate',1e-6,'refine',[],'seeds',[11 22 33],'verbose',true, ...
             'w0',[],'lmit',400,'lmtol',[]);                                % CC, 09/22/2026
fn = fieldnames(def);
for k = 1:numel(fn)
    if ~isfield(opts,fn{k}) || isempty(opts.(fn{k})), opts.(fn{k}) = def.(fn{k}); end
end
if isempty(which('poslpivar_2d'))
    error('pielr_certify:path','PIETOOLS is not on the path; run pielr_path first.');
end
if isempty(opts.settings), opts.settings = set2d_deg(4,[]); end
if isempty(opts.refine)
    opts.refine = isempty(opts.face);   % shrink only when we searched ourselves
end
if isempty(opts.lmtol)                  % raw-residual target for bm_lm2      % CC, 09/22/2026
    opts.lmtol = opts.gate;             % see the bm_lm2 header for the ratio % CC, 09/22/2026
end
okg = @(R) (R.rel < opts.gate) && R.psd;   % THE acceptance gate, used everywhere
vb  = opts.verbose;
notes = {};

% ---- normalise the PIE input ----------------------------------------------
PIE = norm_pie(PIE);

% Start discovery at the ESTIMATED rank instead of r = 1 (maintainer request: % CC, 09/19/2026
% the r < 2*n2 sweeps always failed and wasted seeds*LM time).  Measured law: % CC, 09/19/2026
% per-block minimum rank = 2*n2 (n2 = number of 2-D PDE states) at n2 = 1..6  % CC, 09/19/2026
% on the replicated family, [2 2] on Heat_Eq_with_ODE, [3 3] (one above) on   % CC, 09/19/2026
% the anisotropic variant.  Overshooting is repaired for free by the refine   % CC, 09/19/2026
% shrink; undershoot escalates, so maxrank is lifted to leave headroom.       % CC, 09/19/2026
if isempty(opts.rank) && isempty(opts.face)                                   % CC, 09/19/2026
    n2est = PIE.T.dim(4,1);                                                   % CC, 09/19/2026
    opts.rank = max(1,2*n2est);                                               % CC, 09/19/2026
    opts.maxrank = max(opts.maxrank,opts.rank+2);                             % CC, 09/19/2026
    notes{end+1} = sprintf(['discovery started at estimated rank %d = 2 x ' ...% CC, 09/19/2026
        '%d L2 states (measured law); ranks below were not searched blind ' ...% CC, 09/19/2026
        '-- the refine shrink recovers them when they exist'], ...             % CC, 09/19/2026
        opts.rank,n2est);                                                      % CC, 09/19/2026
    if vb, fprintf('  starting discovery at estimated rank %d (2 x %d L2 states)\n', ... % CC, 09/19/2026
                   opts.rank,n2est); end                                       % CC, 09/19/2026
end                                                                            % CC, 09/19/2026

% ---- SETUP: LPI assembly + SDP data ----------------------------------------
tS = tic;
if vb, fprintf('pielr_certify: assembling the 2-D stability LPI ...\n'); end
[prog,H] = build_stab_2d_st2(PIE,opts.settings);
[Atf,bf,Ns,Kf] = raw_data(prog);
% pre=0: the BM preconditioner needs an eig of a dense m x m Gram; it is built
% inside the BM route only, and its cost is charged to discovery
P = bm_setup(Atf,bf,Ns,Kf,0);
B = numel(Ns);
[part,pnote] = local_partition(H,opts.settings,Ns);
if ~isempty(pnote), notes{end+1} = pnote; end
cert = struct();
cert.t_setup = toc(tS);
if vb
    fprintf('  blocks {%s}  Ns=[%s]  m=%d  Gram unknowns=%d  [setup %.1fs]\n', ...
        strjoin(H.tags,','),num2str(Ns),P.m,sum(Ns.*(Ns+1)/2),cert.t_setup);
end

% ---- CERTIFY / DISCOVER -----------------------------------------------------
t_discover = 0;  t_certify = 0; %#ok<NASGU> % t_certify set on every branch below
V = [];  R = [];  q = [];  route = '';
if ~isempty(opts.face)
    % ------- a face was supplied: certification only -------------------------
    face = opts.face;
    if isfield(face,'face') && ~isfield(face,'V')                          % CC, 09/19/2026
        % a previous pielr_certify OUTPUT was passed: remap cert.face/.S   % CC, 09/19/2026
        assert(~isfield(face,'ok') || face.ok, ...                         % CC, 09/19/2026
            'pielr_certify: the supplied certificate has ok=false; it carries no usable face.'); % CC, 09/19/2026
        fV = face.face;                                                    % CC, 09/19/2026
        if isfield(face,'S'), fS = face.S; else, fS = []; end              % CC, 09/19/2026
        face = struct();  face.V = fV;  face.S = fS;                       % CC, 09/19/2026
    end                                                                    % CC, 09/19/2026
    assert(iscell(face.V) && numel(face.V)==B, ...
        'pielr_certify: face has %d blocks, program has %d',numel(face.V),B);
    for i = 1:B
        assert(size(face.V{i},1)==Ns(i), ...
            ['pielr_certify: face block %d has %d rows but this program''s Gram ' ...
             'block is %d x %d -- wrong face (or wrong settings/state count)'], ...
            i,size(face.V{i},1),Ns(i),Ns(i));
    end
    tC = tic;
    if isfield(face,'S') && ~isempty(face.S)
        % direct lift X_i = V_i S_i V_i' (ORIGINAL units) -- no SDP at all.
        % Measured: the replicated heat certificate verifies this way at
        % n = 1..6 at the n=1 residual, where the tiny-SDP re-solve started
        % missing the gate at n=4.
        q = zeros(P.Ntot,1);
        for i = 1:B
            Xi = face.V{i}*face.S{i}*face.V{i}';   Xi = (Xi+Xi')/2;
            q(P.rows{i}) = Xi(:)/P.nb0;            % opcheck: normalised-b units
        end
        R = opcheck_2d(prog,H,P,q);
        if okg(R)
            V = face.V;  route = 'face-replicate';
            notes{end+1} = 'certified by direct lift of the supplied face coefficients; no SDP solved';
        elseif vb
            fprintf('  direct lift missed the gate (op rel %.3g); re-solving on the face\n',R.rel);
        end
    end
    if isempty(route)
        [R,q,info] = restrict_solve(prog,H,P,Atf,bf,face.V); %#ok<ASGLU>
        if okg(R)
            V = face.V;  route = 'face-sdp';
        else
            % keep R for the report; cert.ok will be false
            V = face.V;  route = 'face-sdp';
            notes{end+1} = sprintf(['supplied face does NOT certify here ' ...
                '(op rel %.3g, mineig %.3g) -- not reached on this face; ' ...
                'not a lower bound'],R.rel,min(R.mineig));
        end
    end
    t_certify = toc(tC);
else
    % ------- cold discovery ---------------------------------------------------
    route = opts.route;
    privdir = fullfile(fileparts(mfilename('fullpath')),'private');
    bmfiles = {'bm_setup','bm_proj','bm_dr','bm_resid','bm_report','bm_lm2'};
    hasbm = all(cellfun(@(f)exist(fullfile(privdir,[f '.m']),'file')==2,bmfiles));
    if strcmp(route,'auto'), if hasbm, route = 'bm'; else, route = 'mintrace'; end, end
    if strcmp(route,'bm') && ~hasbm
        notes{end+1} = 'BM helpers missing; fell back to mintrace';
        route = 'mintrace';
    end
    tD = tic;
    switch route
        case 'bm'
            [V,R,q,t_cert_core,dnotes] = disc_bm(prog,H,Atf,bf,Ns,Kf,opts,okg,vb);
        case 'mintrace'
            [V,R,q,t_cert_core,dnotes] = disc_mintrace(prog,H,P,Atf,bf,Ns,Kf,opts,okg,vb);
        otherwise
            error('pielr_certify:route','unknown route ''%s''',route);
    end
    notes = [notes,dnotes];
    t_discover = toc(tD) - t_cert_core;
    t_certify  = t_cert_core;
end

% ---- optional per-block refine (charged to certification) -------------------
if opts.refine && ~isempty(R) && okg(R)
    tR = tic;
    [V,R,q,rnotes] = refine_blocks(prog,H,P,Atf,bf,V,R,q,okg,vb);
    notes = [notes,rnotes];
    t_certify = t_certify + toc(tR);
end

% ---- package the result -----------------------------------------------------
cert.ok = ~isempty(R) && okg(R);
if ~isempty(R)
    cert.op_rel = R.rel;   cert.mineig = R.mineig;   cert.rank = R.rank;
    cert.psd = R.psd;
else
    cert.op_rel = NaN;  cert.mineig = NaN;  cert.rank = NaN;  cert.psd = false;
end
if ~isempty(V)
    cert.r = cellfun(@(v)size(v,2),V);
    cert.face = V;
    cert.unknowns_face = sum(cert.r.*(cert.r+1)/2);
else
    cert.r = NaN;  cert.face = {};  cert.unknowns_face = NaN;
end
% face coefficients in ORIGINAL units, so the face can be re-used / tensored
cert.S = cell(1,B);
if ~isempty(V) && ~isempty(q)
    for i = 1:B
        Xi = reshape(q(P.rows{i}),Ns(i),Ns(i))*P.nb0;
        Si = V{i}'*((Xi+Xi')/2)*V{i};
        cert.S{i} = (Si+Si')/2;
    end
end
cert.Ns = Ns;  cert.Kf = Kf;  cert.nb0 = P.nb0;
cert.unknowns_full = sum(Ns.*(Ns+1)/2);
cert.t_discover = t_discover;   cert.t_certify = t_certify;
cert.route = route;   cert.part = part;   cert.gate = opts.gate;
if ~cert.ok
    notes{end+1} = sprintf(['NO certificate at rank <= %d passed the gate: ' ...
        'reported as NOT REACHED by this search, never as a proven floor'],opts.maxrank);
end
notes{end+1} = 'all ranks are upper bounds exhibited by verified certificates';
cert.notes = notes;

% ---- one-screen report ------------------------------------------------------
if vb
    fprintf('\n================ pielr_certify report ================\n');
    fprintf('dims [n0 nx ny n2] = [%s]   blocks {%s}  Ns = [%s]\n', ...
        num2str(PIE.T.dim(:,1)'),strjoin(H.tags,','),num2str(Ns));
    fprintf('route  : %s\n',cert.route);
    if cert.ok
        fprintf('ranks  : r = %s of %s   (Gram unknowns %d -> %d, %.0fx fewer)\n', ...
            mat2str(cert.r),mat2str(Ns),cert.unknowns_full,cert.unknowns_face, ...
            cert.unknowns_full/max(cert.unknowns_face,1));
        fprintf('verdict: CERTIFIES   op rel = %.4g (gate %.1g)   PSD = %d   mineig = %s\n', ...
            cert.op_rel,cert.gate,cert.psd,mat2str(cert.mineig,3));
    else
        fprintf('verdict: NOT CERTIFIED (op rel = %.4g, psd = %d) -- see notes\n', ...
            cert.op_rel,cert.psd);
    end
    fprintf('times  : setup %.1fs | discovery %.1fs | certification %.1fs\n', ...
        cert.t_setup,cert.t_discover,cert.t_certify);
    for k = 1:numel(cert.notes), fprintf('note   : %s\n',cert.notes{k}); end
    fprintf('======================================================\n');
end
end

% =============================================================================
function PIE = norm_pie(PIE)
% Accept a pie_struct or a plain struct; require T and A; default vars/dom
% from the T operator so raw opvar2d objects can be handed in directly.
try
    Top = PIE.T;   Aop = PIE.A;
catch
    error('pielr_certify:input','PIE must carry T and A (opvar2d objects).');
end
if ~isa(Top,'opvar2d') || ~isa(Aop,'opvar2d')
    error('pielr_certify:input', ...
        'PIE.T and PIE.A must be opvar2d objects; this package is 2-D only.');
end
% field assignment, not struct('T',Top,...): with an opvar2d argument that
% call dispatches to opvar2d's own struct method and errors
S = struct();
S.T = Top;   S.A = Aop;
try
    S.vars = PIE.vars;
catch
    S.vars = [];
end
try
    S.dom = PIE.dom;
catch
    S.dom = [];
end
if isempty(S.vars), S.vars = [Top.var1, Top.var2]; end
if isempty(S.dom),  S.dom  = Top.I; end
PIE = S;
end

% =============================================================================
function [part,note] = local_partition(H,s2d,Ns)
% Per-block Gram row partition into poslpivar_2d monomial groups, CHECKED
% against the actual Gram side lengths.  Needed only by pielr_tensor (state-
% ladder replication); {} plus a note when the layout does not tile that way.
part = {};  note = '';
try
    eqo = get_eq_opts_2D(H.Qop,s2d.eq_opts,1e-12);
    RL = nb_2d(s2d.LF_deg,H.Top.dim(:,1)',getf(s2d.LF_opts,'sep',0), ...
               getf(s2d.LF_opts,'exclude',zeros(1,16)));
    RE = nb_2d(s2d.eq_deg,H.Qop.dim(:,1)',getf(eqo,'sep',0), ...
               getf(eqo,'exclude',zeros(1,16)));
    RR = {RL,RE};
    pt = cell(1,numel(Ns));
    for i = 1:numel(Ns)
        Rj = RR{i};   gi = find(Rj.inc & Rj.nZ>0 & Rj.Zdim>0);
        pt{i} = Rj.nZ(gi);
        % sum|Z| = Ns holds only when every active group carries state-dim 1
        % (the n = 1 member of a replicated family); anything else cannot be
        % state-tensored, so refuse the partition rather than mislead
        if sum(pt{i})~=Ns(i)
            note = sprintf(['group partition does not tile the Gram ' ...
                '(block %d: sum|Z|=%d vs Ns=%d); pielr_tensor unavailable ' ...
                'for this certificate'],i,sum(pt{i}),Ns(i));
            return
        end
    end
    part = pt;
catch ME
    note = sprintf('group partition unavailable (%s)',ME.message);
end
end

% =============================================================================
function [V,R,q,t_cert_core,notes] = disc_bm(prog,H,Atf,bf,Ns,Kf,opts,okg,vb)
% Burer-Monteiro discovery + face certification, the measured pipeline:
% BM explores unconfined (its residual is only approximate), restrict_solve
% converts the SUBSPACE it found into an exact certificate (sound one-way:
% S_i >= 0 gives X_i >= 0, so a restriction can only lose feasibility).
notes = {};  V = [];  R = [];  q = [];  t_cert_core = 0;
B = numel(Ns);
if vb, fprintf('  BM discovery: building the preconditioner (eig of %d x %d) ...\n', ...
        size(Atf,2),size(Atf,2)); end
P = bm_setup(Atf,bf,Ns,Kf,1);          % pre=1: BM needs the whitened metric
% rank ladder: opts.rank (if any) first, then the common ranks above it
rvl = {};
if ~isempty(opts.rank)
    rv0 = opts.rank(:)';
    if isscalar(rv0), rv0 = rv0*ones(1,B); end
    assert(numel(rv0)==B,'opts.rank must be scalar or one entry per block');
    rvl{end+1} = min(rv0,Ns);
    rnext = max(rv0)+1;
else
    rnext = 1;
end
for r = rnext:opts.maxrank, rvl{end+1} = min(r,Ns); end %#ok<AGROW>
if vb
    fprintf('  %-11s | %-8s | %-11s | %-11s | %-11s | %s\n', ...
        'rank','start','BM raw_rel','BM op rel','face op rel','verdict');
end
prevw = [];  prevrv = [];
for ri = 1:numel(rvl)
    rv = rvl{ri};
    bestraw = inf;  bestw = [];
    % CC, 09/21/2026: a supplied .w0 is tried FIRST at every rung.  It must
    % lead: on the heat family seed 1 certifies on the first attempt, so an
    % arm appended after the seeds would never execute and the option would
    % silently measure nothing.
    nw0 = double(~isempty(opts.w0));                                        % CC, 09/21/2026
    nst = nw0 + numel(opts.seeds) + ~isempty(prevw);                        % CC, 09/21/2026
    for si = 1:nst
        if nw0 && si==1                                                     % CC, 09/21/2026
            w = fitw(opts.w0,rv,Ns);   sname = 'w0';                        % CC, 09/21/2026
        elseif si <= nw0 + numel(opts.seeds)                                % CC, 09/21/2026
            sj = si - nw0;                                                  % CC, 09/21/2026
            rng(opts.seeds(sj),'twister');                                  % CC, 09/21/2026
            scl = [1e-1 1 1e1 1e2 1e4 1e6];  sc = scl(mod(sj-1,numel(scl))+1); % CC, 09/21/2026
            q0 = zeros(P.Ntot,1);
            for i = 1:B
                Mr = randn(Ns(i));  Mr = sc*(Mr+Mr')/sqrt(2*Ns(i));
                q0(P.rows{i}) = Mr(:);
            end
            q0 = bm_proj(P,q0,rv,'affine');
            [Vd,~] = bm_dr(P,rv,q0,300);
            w = [];
            for i = 1:B
                Yi = zeros(Ns(i),rv(i));
                if ~isempty(Vd) && ~isempty(Vd{i})
                    kk = min(size(Vd{i},2),rv(i));  Yi(:,1:kk) = Vd{i}(:,1:kk);
                end
                w = [w;Yi(:)]; %#ok<AGROW>
            end
            sname = sprintf('seed%d',opts.seeds(sj));                       % CC, 09/21/2026
        else
            % warm start: a rank-r' solution with r' < r is also a rank-r
            % solution, so padding keeps the sweep monotone in r
            w = padw(prevw,prevrv,rv,Ns);   sname = 'prev';
        end
%       w = bm_lm2(P,rv,w,400,1e-16);                                       % CC, 09/21/2026 (was)
        % tol at the gate value, not 1e-16: op_rel ~ 0.35-0.40*raw_rel here,  % CC, 09/22/2026
        % so this exits with ~2.5x margin instead of never exiting at all.    % CC, 09/22/2026
        w = bm_lm2(P,rv,w,opts.lmit,opts.lmtol);                            % CC, 09/22/2026
        R1 = bm_report(w,P,rv);
        if R1.raw_rel < bestraw, bestraw = R1.raw_rel;  bestw = w; end
        [~,~,qk] = bm_resid(w,P,rv);
        Rk = opcheck_2d(prog,H,P,qk);
        % certify the SUBSPACE BM found -- the point of the find/certify split
        Vc = cell(1,B);
        for i = 1:B
            Yi = reshape(wblk(w,Ns,rv,i),Ns(i),rv(i));
            Vi = orth(Yi);
            if isempty(Vi)                 % Y_i == 0: keep a 1-dim face, which
                Vi = zeros(Ns(i),1); Vi(1) = 1;   % still contains X_i = 0
            end
            Vc{i} = Vi;
        end
        Rf = [];  qf = [];  frel = NaN;  tcc = 0;
        try
            tc0 = tic;
            [Rf,qf] = restrict_solve(prog,H,P,Atf,bf,Vc);
            tcc = toc(tc0);   frel = Rf.rel;
        catch ME
            if vb, fprintf('      (restrict_solve: %s)\n',ME.message(1:min(70,end))); end
        end
        okf = ~isempty(Rf) && okg(Rf);
        if vb
            fprintf('  %-11s | %-8s | %-11.4g | %-11.4g | %-11.4g | %s\n', ...
                mat2str(rv),sname,R1.raw_rel,Rk.rel,frel, ...
                tern(okf,'FACE OK',tern(okg(Rk),'BM OK','fail')));
        end
        if okf
            V = Vc;  R = Rf;  q = qf;  t_cert_core = tcc;  return
        elseif okg(Rk)
            % the BM point itself passes the operator gate
            V = Vc;  R = Rk;  q = qk;  t_cert_core = 0;
            notes{end+1} = 'accepted the BM point directly (operator gate passed without a face re-solve)'; %#ok<AGROW>
            return
        end
    end
    prevw = bestw;  prevrv = rv;
end
notes{end+1} = sprintf('BM found no certificate at rank <= %d over %d seeds/rank', ...
                       opts.maxrank,numel(opts.seeds));
end

% =============================================================================
function [V,R,q,t_cert_core,notes] = disc_mintrace(prog,H,P,Atf,bf,Ns,Kf,opts,okg,vb)
% Min-trace discovery: one full SeDuMi solve of the LPI with a trace objective
% (trace minimisation is the truncation source -- the analytic centre of the
% plain feasibility program is near-full-rank), then per-block truncation to
% the leading eigenvectors and face certification via restrict_solve.
notes = {};  V = [];  R = [];  q = [];  t_cert_core = 0;
B = numel(Ns);
c = zeros(P.Ntot,1);
off = Kf;
for i = 1:B
    Ei = eye(Ns(i));   c(off+(1:Ns(i)^2)) = Ei(:);   off = off + Ns(i)^2;
end
K = struct('f',Kf,'s',Ns);
pars.fid = 0;
if vb, fprintf('  mintrace discovery: full SeDuMi solve (m=%d, blocks [%s]) ...\n', ...
        P.m,num2str(Ns)); end
t0 = tic;
% solve with b at UNIT NORM: ||b|| ~ 5e-6 on this family, and SeDuMi's
% termination measures are relative to the data scale, so the raw b loses
% digits exactly where the certificate needs them (measured; the same reason
% bm_setup normalises).  The program is linear, so x = nb0 * y rescales back.
[y,~,infm] = sedumi(Atf,full(bf(:))/P.nb0,sparse(c),K,pars);
if vb, fprintf('  mintrace solve %.1fs  (SeDuMi numerr %d)\n',toc(t0),getf(infm,'numerr',-1)); end
x = full(y)*P.nb0;
% per-block eigendecomposition of the min-trace point (ORIGINAL units)
W = cell(1,B);  ev = cell(1,B);
off = Kf;
for i = 1:B
    Xi = reshape(x(off+(1:Ns(i)^2)),Ns(i),Ns(i));   Xi = (Xi+Xi')/2;
    [Wi,Di] = eig(Xi);   [di,ord] = sort(diag(Di),'descend');
    W{i} = Wi(:,ord);  ev{i} = di;   off = off + Ns(i)^2;
    if vb
        nr = sum(di > max(di(1),eps)*1e-8);
        fprintf('  blk %d: mintrace numerical rank %d of %d (top eigs %s)\n', ...
            i,nr,Ns(i),mat2str(di(1:min(4,end))',2));
    end
end
% rank ladder on the truncated faces
rvl = {};
if ~isempty(opts.rank)
    rv0 = opts.rank(:)';  if isscalar(rv0), rv0 = rv0*ones(1,B); end
    rvl{end+1} = min(rv0,Ns);   rnext = max(rv0)+1;
else
    rnext = 1;
end
for r = rnext:opts.maxrank, rvl{end+1} = min(r,Ns); end %#ok<AGROW>
for ri = 1:numel(rvl)
    rv = rvl{ri};
    Vc = cell(1,B);
    for i = 1:B, Vc{i} = W{i}(:,1:rv(i)); end
    Rf = [];  qf = [];  tcc = 0;
    try
        tc0 = tic;
        [Rf,qf] = restrict_solve(prog,H,P,Atf,bf,Vc);
        tcc = toc(tc0);
    catch ME
        if vb, fprintf('      (restrict_solve: %s)\n',ME.message(1:min(70,end))); end
    end
    if vb && ~isempty(Rf)
        fprintf('  truncate r=%s : op rel %.4g  mineig %.3g  %s\n', ...
            mat2str(rv),Rf.rel,min(Rf.mineig),tern(okg(Rf),'CERTIFIES','fail'));
    end
    if ~isempty(Rf) && okg(Rf)
        V = Vc;  R = Rf;  q = qf;  t_cert_core = tcc;  return
    end
end
notes{end+1} = sprintf('mintrace truncation found no certificate at rank <= %d',opts.maxrank);
end

% =============================================================================
function [V,R,q,notes] = refine_blocks(prog,H,P,Atf,bf,V,R,q,okg,vb)
% Shrink the certified face one direction at a time, per block.  Directions are
% ordered by the certified S_i's own spectrum (measured weight, not guessed);
% each trial is a tiny SDP.  A rejected shrink is floor EVIDENCE, not a proof.
notes = {};
B = numel(P.Ns);
improved = true;
while improved
    improved = false;
    for i = 1:B
        if size(V{i},2) <= 1, continue; end
        Vt = V;   Vt{i} = shrink1(V{i},q,P,i);
        try
            [Rt,qt] = restrict_solve(prog,H,P,Atf,bf,Vt);
            if okg(Rt)
                V = Vt;  R = Rt;  q = qt;  improved = true;
                if vb, fprintf('  refine blk %d -> r=%d  OK    op rel %.3g\n', ...
                        i,size(Vt{i},2),Rt.rel); end
            else
                if vb, fprintf('  refine blk %d -> r=%d  FAILS op rel %.3g mineig %.2g  <== resists\n', ...
                        i,size(Vt{i},2),Rt.rel,min(Rt.mineig)); end
            end
        catch ME
            if vb, fprintf('  refine blk %d ERROR %s\n',i,ME.message(1:min(70,end))); end
        end
    end
end
end

% =============================================================================
function Vs = shrink1(Vi,q,P,i)
% drop the direction carrying the least weight in the CERTIFIED solution:
% diagonalise S_i = V_i' X_i V_i and keep its top eigenvectors
N = P.Ns(i);
X = reshape(q(P.rows{i}),N,N);   X = (X+X')/2;
Si = Vi'*X*Vi;   Si = (Si+Si')/2;
[W,D] = eig(Si);   [~,p] = sort(diag(D),'descend');   W = W(:,p);
Vs = Vi*W(:,1:end-1);
[Vs,~] = qr(Vs,0);               % re-orthonormalise
end

function v = wblk(w,Ns,rv,i)
k = 0;
for j = 1:i-1, k = k + Ns(j)*rv(j); end
v = w(k+(1:Ns(i)*rv(i)));
end

function w = fitw(Y0,rv,Ns)                                                 % CC, 09/21/2026
% Fit a supplied per-block factor {Y_1..Y_B} to the current rank profile rv.
% Keep the leading columns; fill any surplus with SMALL RANDOM values scaled
% to the factor's own entries.  Never zeros: a zero column of Y has an
% identically zero Jacobian block (J_i = 2*W*Ssym(:,rows_i)*kron(Y_i,I), whose
% column block for column c is Y_i(:,c) (x) I), so LM can never move it and the
% extra rank would be inert -- the same fixed point the 'prev' pad sits on.
w = [];
for i = 1:numel(Ns)
    Yi = zeros(Ns(i),rv(i));
    Y  = Y0{i};
    kk = min(size(Y,2),rv(i));
    if kk>0, Yi(:,1:kk) = Y(:,1:kk); end
    if rv(i) > kk
        s = norm(Y,'fro')/max(sqrt(numel(Y)),1);       % rms entry of the seed
        if ~isfinite(s) || s<=0, s = 1; end
        Yi(:,kk+1:end) = 1e-3*s*randn(Ns(i),rv(i)-kk);
    end
    w = [w;Yi(:)]; %#ok<AGROW>
end
end

function w = padw(wold,rold,rnew,Ns)
w = [];  k = 0;
for i = 1:numel(Ns)
    Yo = reshape(wold(k+(1:Ns(i)*rold(i))),Ns(i),rold(i));  k = k+Ns(i)*rold(i);
    Yn = zeros(Ns(i),rnew(i));  cc = min(rold(i),rnew(i));  Yn(:,1:cc) = Yo(:,1:cc);
    w = [w;Yn(:)]; %#ok<AGROW>
end
end

function s = tern(c,a,b)
if c, s = a; else, s = b; end
end

function v = getf(s,f,d)
if isstruct(s) && isfield(s,f) && ~isempty(s.(f)), v = s.(f); else, v = d; end
end
