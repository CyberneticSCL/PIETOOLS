function cert = pielr_certify_pos(Ptgt,opts)
% PIELR_CERTIFY_POS  Low-rank certificate that a GIVEN 2-D PI operator is >= 0.
%
%   cert = pielr_certify_pos(Ptgt)
%   cert = pielr_certify_pos(Ptgt,opts)
%   cert = pielr_certify_pos(Ptgt,prev_cert)   % a previous output IS a valid
%   cert = pielr_certify_pos(Ptgt,face)        % input: certification only
%
% The companion of pielr_certify for BARE POSITIVITY: instead of building the
% stability LPI from (T,A), the program here is
%     find  Q >= 0  such that  Ptgt = Z' Q Z     (every opvar2d cell matched)
% i.e. one Gram block, equated to the operator you supply.  Weighted integral
% inequalities are the canonical use (see demo_poincare_2d.m); note the
% structural difference from stability problems: here the certificate is
% FORCED to realise every cell of Ptgt, including a multiplier cell
% R22{1,1} = a(s1,s2) if the target has one -- and then the Gram sub-block
% carrying the multiplier has rank >= the number of squares a needs (measured:
% exactly that number, 1/2/3 on a designed ladder).  A target whose multiplier
% cell is nonnegative on the box but NOT a global sum of squares is INFEASIBLE
% at these settings (psatz 0) at any operating point; that is a property of
% the formulation, not a bug (theory.pdf, Section 7).
%
% INPUT
%   Ptgt  opvar2d, square (dim(:,1) == dim(:,2)), self-adjoint as an operator
%         (Ptgt' == Ptgt); spatial variables and domain are read off Ptgt.
%   opts  struct, all fields optional -- or a previous cert / a bare face,
%         accepted directly exactly as in pielr_certify:
%     .spec     [n1 n2 n3] Gram degree spec (default [2 2 1]): n1 = individual
%               degree of the multiplier basis Z2oo in each variable (full
%               tensor, (n1+1)^2 monomials); n2/n3 = primary/dummy individual
%               degrees of the four double-integral bases.  This default is a
%               STARTING spec, not derived from the target: if the equality is
%               infeasible at it, raise the spec (the §7 escalation of
%               GUIDE.md applies here too).
%     .deg      full poslpivar_2d degree struct {dx,dy,d2}; overrides .spec.
%     .rank     rank to try first in discovery.
%     .face     warm-start face (struct with .V and optionally .S, or a
%               previous cert); skips discovery.
%     .maxrank  discovery search ceiling (default 6).
%     .gate     operator-residual acceptance gate (default 1e-6).
%     .seeds    BM random seeds per rank (default [11 22 33]).
%     .refine   within-face shrink after acceptance (default: true after
%               discovery, false when a face was supplied).
%     .ref      also run a full-rank SeDuMi reference solve and use its
%               truncation as an extra BM start (default false; it is the
%               expensive diagnostic, and its failure does not gate anything).
%     .drit / .lmit   BM inner iteration budgets (defaults 200 / 150).
%     .verbose  print progress and the report (default true).
%
% OUTPUT: the same cert struct as pielr_certify (ok, r, rank, op_rel, mineig,
% face, S, Ns, unknowns_full/face, t_setup/t_discover/t_certify, route,
% notes), with cert.part = {} -- tensor replication across states does not
% apply to a bare-positivity certificate.
%
% ACCEPTANCE is opcheck_pos2d.ok only:  max|Ptgt - Peop(q)| / max|Ptgt| < gate
% over ALL 36 opvar2d leaf cells AND the Gram PSD.  Never an equality residual,
% never a solver flag.  Every reported rank is an UPPER BOUND exhibited by a
% verified certificate; a failed rank is "not reached on the faces tried" --
% except that when the target HAS a multiplier cell, rank(Q) >= the SOS length
% of that cell is an analytic LOWER bound you can compute by hand.
%
% Component structure is derived from the TARGET by PIETOOLS' own
% get_eq_opts_2D: Gram components whose target cells are identically zero are
% excluded, so the basis matches the operator instead of carrying zero-forced
% rows.  psatz is held at 0 throughout.
%
% Requires PIETOOLS and SeDuMi on the path: run pielr_path once per session.
%
% CC, 09/19/2026: initial version -- the ft_rank.m pipeline (BM find ->
%                 restrict_solve_pos certify -> within-face shrink), packaged
%                 with the pielr_certify option surface and cert struct.

% ---- options (mirrors pielr_certify, incl. direct cert/face second arg) ----
if nargin<2, opts = struct(); end
if isstruct(opts) && isfield(opts,'ok') && isfield(opts,'face')
    opts = struct('face',opts);
elseif isstruct(opts) && isfield(opts,'V') && iscell(opts.V) && ~isfield(opts,'face')
    opts = struct('face',opts);
end
def = struct('spec',[2 2 1],'deg',[],'rank',[],'face',[],'maxrank',6, ...
             'gate',1e-6,'seeds',[11 22 33],'refine',[],'ref',false, ...
             'drit',200,'lmit',150,'verbose',true);
fn = fieldnames(def);
for k = 1:numel(fn)
    if ~isfield(opts,fn{k}) || isempty(opts.(fn{k})), opts.(fn{k}) = def.(fn{k}); end
end
if isempty(which('poslpivar_2d'))
    error('pielr_certify_pos:path','PIETOOLS is not on the path; run pielr_path first.');
end
if ~isa(Ptgt,'opvar2d')
    error('pielr_certify_pos:input', ...
        'the target must be an opvar2d object; this package is 2-D only.');
end
dimc = Ptgt.dim;
assert(all(dimc(:,1)==dimc(:,2)), ...
    'pielr_certify_pos: the target must be square (dim(:,1) == dim(:,2)).');
okg = @(R) (R.rel < opts.gate) && R.psd;
vb  = opts.verbose;
notes = {};
if isempty(opts.rank) && isempty(opts.face)                                  % CC, 09/19/2026
    % rank 1 never certified in any measured 2-D bare-positivity case (the    % CC, 09/19/2026
    % measured floors were 3/4-5/5-6); start at 2 and escalate, refine shrinks % CC, 09/19/2026
    opts.rank = 2;                                                            % CC, 09/19/2026
    notes{end+1} = 'discovery started at rank 2 (measured: rank 1 never certifies here)'; % CC, 09/19/2026
end                                                                           % CC, 09/19/2026

% ---- SETUP: the equality program --------------------------------------------
tS = tic;
if vb, fprintf('pielr_certify_pos: assembling the bare-positivity LPI ...\n'); end
prog = lpiprogram(Ptgt.var1,Ptgt.var2,Ptgt.I);
if ~isempty(opts.deg)
    dd = opts.deg;
else
    dd = spec2deg(opts.spec);
end
gopts = struct();  gopts.psatz = 0;
gopts.exclude = zeros(1,16);  gopts.sep = zeros(1,6);
% derive the component structure from the target, exactly as the shipped 2-D
% executive derives its negativity-Gram structure from Qop
gopts = get_eq_opts_2D(Ptgt,gopts,1e-12);
[prog,Peop] = poslpivar_2d(prog,dimc(:,1),dd,gopts);
prog = lpi_eq(prog,Ptgt-Peop,'symmetric');
[Atf,bf,Ns,Kf] = raw_data(prog);
assert(isscalar(Ns), ...
    'pielr_certify_pos: expected one Gram block, got %d.',numel(Ns));
P = bm_setup(Atf,bf,Ns,Kf,isempty(opts.face));  % BM preconditioner only if searching
cert = struct();
cert.t_setup = toc(tS);
if vb
    fprintf('  Gram N=%d  m=%d  unknowns=%d  [setup %.1fs]\n', ...
            Ns(1),P.m,sum(Ns.*(Ns+1)/2),cert.t_setup);
end

% ---- CERTIFY / DISCOVER -----------------------------------------------------
t_discover = 0;  t_certify = 0;
V = [];  R = [];  q = [];  route = '';
if ~isempty(opts.face)
    face = opts.face;
    if isfield(face,'face') && ~isfield(face,'V')
        assert(~isfield(face,'ok') || face.ok, ...
            'pielr_certify_pos: the supplied certificate has ok=false; it carries no usable face.');
        fV = face.face;
        if isfield(face,'S'), fS = face.S; else, fS = []; end
        face = struct();  face.V = fV;  face.S = fS;
    end
    assert(iscell(face.V) && isscalar(face.V), ...
        'pielr_certify_pos: face must carry exactly one block.');
    assert(size(face.V{1},1)==Ns(1), ...
        'pielr_certify_pos: face has %d rows but the Gram is %d x %d.', ...
        size(face.V{1},1),Ns(1),Ns(1));
    tC = tic;
    if isfield(face,'S') && ~isempty(face.S)
        % direct lift, no SDP (S in ORIGINAL units; opcheck wants /nb0)
        q = zeros(P.Ntot,1);
        Xi = face.V{1}*face.S{1}*face.V{1}';   Xi = (Xi+Xi')/2;
        q(P.rows{1}) = Xi(:)/P.nb0;
        R = opcheck_pos2d(prog,Ptgt,Peop,P,q);
        if okg(R)
            V = face.V;  route = 'face-replicate';
            notes{end+1} = 'certified by direct lift of the supplied face coefficients; no SDP solved';
        elseif vb
            fprintf('  direct lift missed the gate (op rel %.3g); re-solving on the face\n',R.rel);
        end
    end
    if isempty(V)
        [R,q] = restrict_solve_pos(prog,Ptgt,Peop,P,Atf,bf,face.V);
        if okg(R), V = face.V;  route = 'face-sdp'; end
    end
    t_certify = toc(tC);
    if isempty(V)
        notes{end+1} = sprintf('the supplied face did not certify (op rel %.3g)',R.rel);
    end
else
    % ------- discovery: the measured ft_rank pipeline -------------------------
    tD = tic;
    q0 = [];
    if opts.ref
        q0 = sed_ref_pos(Atf,bf,Ns,Kf,P);
        if ~isempty(q0)
            R0 = opcheck_pos2d(prog,Ptgt,Peop,P,q0);
            if vb, fprintf('  reference solve: op rel %.3g psd %d\n',R0.rel,R0.psd); end
            if R0.rel >= opts.gate, q0 = []; end   % init only, never a gate
        end
    end
    r0 = 1;  if ~isempty(opts.rank), r0 = opts.rank(1); end               % CC, 09/19/2026
    % start AT the requested/estimated rank and escalate; ranks below are    % CC, 09/19/2026
    % recovered by the refine shrink when they exist, never swept blind     % CC, 09/19/2026
    rlist = r0:min(opts.maxrank,Ns(1));                                    % CC, 09/19/2026
    prevw = [];  prevr = [];
    for r = rlist
        rv = min(r,Ns);
        bestraw = inf;  bestw = [];
        nst = numel(opts.seeds) + double(~isempty(prevw)) + double(~isempty(q0));
        for si = 1:nst
            if si <= numel(opts.seeds)
                rng(opts.seeds(si),'twister');
                scl = [1e-1 1 1e1];  sc = scl(mod(si-1,3)+1);
                qi = zeros(P.Ntot,1);
                Mr = randn(Ns(1));  Mr = sc*(Mr+Mr')/sqrt(2*Ns(1));
                qi(P.rows{1}) = Mr(:);
                qi = bm_proj(P,qi,rv,'affine');
                sname = sprintf('seed%d',opts.seeds(si));
            elseif si == numel(opts.seeds)+1 && ~isempty(q0)
                qi = q0;  sname = 'sedtrunc';
            else
                sname = 'prev';
            end
            if strcmp(sname,'prev')
                w = padw(prevw,prevr,rv,Ns);
            else
                [Vd,~] = bm_dr(P,rv,qi,opts.drit);
                Yi = zeros(Ns(1),rv(1));
                if ~isempty(Vd) && ~isempty(Vd{1})
                    kk = min(size(Vd{1},2),rv(1));  Yi(:,1:kk) = Vd{1}(:,1:kk);
                end
                w = Yi(:);
            end
            w = bm_lm2(P,rv,w,opts.lmit,1e-16);
            Rb = bm_report(w,P,rv);
            if Rb.raw_rel < bestraw, bestraw = Rb.raw_rel;  bestw = w; end
            [~,~,qk] = bm_resid(w,P,rv);
            Rk = opcheck_pos2d(prog,Ptgt,Peop,P,qk);
            Yw = reshape(w,Ns(1),rv(1));
            Vi = orth(Yw);
            if isempty(Vi), Vi = zeros(Ns(1),1);  Vi(1) = 1; end
            Rf = [];  qf = [];
            try
                [Rf,qf] = restrict_solve_pos(prog,Ptgt,Peop,P,Atf,bf,{Vi});
            catch ME
                if vb, fprintf('    (restrict at r=%d: %s)\n',r,ME.message(1:min(60,end))); end
            end
            okf = ~isempty(Rf) && okg(Rf);
            if vb
                if isempty(Rf), frel = NaN; else, frel = Rf.rel; end
                fprintf('  r=%d %-9s BM %.3g -> face %.3g %s\n',r,sname, ...
                        Rk.rel,frel,tern(okf,'FACE OK',tern(okg(Rk),'BM OK','')));
            end
            if okf
                V = {Vi};  R = Rf;  q = qf;  route = 'bm';  break
            elseif okg(Rk)
                V = {orth(Yw)};  R = Rk;  q = qk;  route = 'bm';  break
            end
        end
        prevw = bestw;  prevr = rv;
        if ~isempty(V), break; end
    end
    t_discover = toc(tD);
    if isempty(V)
        notes{end+1} = sprintf(['no rank <= %d reached on the faces tried ' ...
            '(NOT a lower bound; raise .maxrank/.seeds, or the target may be ' ...
            'infeasible at this .spec)'],min(opts.maxrank,Ns(1)));
    end
end

% ---- REFINE: within-face shrink ---------------------------------------------
doref = opts.refine;  if isempty(doref), doref = isempty(opts.face); end
if ~isempty(V) && doref
    tR = tic;
    while size(V{1},2) > 1
        Vs = shrink1(V{1},q,P);
        try
            [Rs,qs] = restrict_solve_pos(prog,Ptgt,Peop,P,Atf,bf,{Vs});
        catch
            break
        end
        if okg(Rs), V = {Vs};  R = Rs;  q = qs;  else, break; end
    end
    t_certify = t_certify + toc(tR);
end

% ---- package + report -------------------------------------------------------
cert.ok = ~isempty(V);
cert.Ns = Ns;  cert.Kf = Kf;  cert.nb0 = P.nb0;
cert.unknowns_full = sum(Ns.*(Ns+1)/2);
cert.t_discover = t_discover;  cert.t_certify = t_certify;
cert.route = route;  cert.part = {};
notes{end+1} = 'tensor replication across states does not apply to bare-positivity certificates';
notes{end+1} = 'all ranks are upper bounds exhibited by verified certificates';
if cert.ok
    cert.r = size(V{1},2);
    cert.face = V;
    Xi = reshape(q(P.rows{1}),Ns(1),Ns(1));  Xi = (Xi+Xi')/2;
    cert.S = {V{1}'*Xi*V{1}*P.nb0};          % ORIGINAL units, like pielr_certify
    cert.rank = R.rank;  cert.op_rel = R.rel;  cert.mineig = R.mineig;
    cert.unknowns_face = cert.r*(cert.r+1)/2;
else
    cert.r = NaN;  cert.face = {};  cert.S = {};
    cert.rank = NaN;  cert.mineig = NaN;  cert.unknowns_face = NaN;
    if ~isempty(R), cert.op_rel = R.rel; else, cert.op_rel = NaN; end
end
cert.notes = notes;
if vb
    fprintf('============ pielr_certify_pos report ============\n');
    fprintf('target dims [%s]   Gram N = %d\n',num2str(dimc(:,1)'),Ns(1));
    if cert.ok
        fprintf('verdict: CERTIFIES   op rel = %.4g (gate %g)   PSD = %d\n', ...
                cert.op_rel,opts.gate,R.psd);
        fprintf('rank   : r = %d of %d   (%d unknowns -> %d)\n', ...
                cert.r,Ns(1),cert.unknowns_full,cert.unknowns_face);
    else
        fprintf('verdict: NOT CERTIFIED (best op rel %.4g)\n',cert.op_rel);
    end
    fprintf('times  : setup %.1fs | discovery %.1fs | certification %.1fs\n', ...
            cert.t_setup,cert.t_discover,cert.t_certify);
    for k = 1:numel(notes), fprintf('note   : %s\n',notes{k}); end
    fprintf('==================================================\n');
end
end

% ---- local helpers -----------------------------------------------------------
function dd = spec2deg(spec)
% the [n1 n2 n3] convenience spec -> full poslpivar_2d degree struct: full
% tensor bases, mirroring the 1-D {n1,[n2 n3 n2+n3],[n2 n3 n2+n3]} pattern in
% each direction (the spec the measured rank ladder used)
n1 = spec(1);  n2 = spec(2);  n3 = spec(3);
dd.dx = {n1;[n2;n3;n2+n3];[n2;n3;n2+n3]};
dd.dy = {n1,[n2,n3,n2+n3],[n2,n3,n2+n3]};
o  = [0,  n1;  n1,  2*n1];
ii = [0,      n2,      n3,      n2+n3;
      n2,     2*n2,    n2+n3,   2*n2+n3;
      n3,     n2+n3,   2*n3,    n2+2*n3;
      n2+n3,  2*n2+n3, n2+2*n3, 2*(n2+n3)];
so = [0, n2, n3, n2+n3; n1, n1+n2, n1+n3, n1+n2+n3];
dd.d2 = {o,   so,  so;
         so', ii,  ii;
         so', ii,  ii};
end
function w = padw(wold,rold,rnew,Ns)
Yo = reshape(wold,Ns(1),rold(1));
Yn = zeros(Ns(1),rnew(1));  c = min(rold(1),rnew(1));  Yn(:,1:c) = Yo(:,1:c);
w = Yn(:);
end
function Vs = shrink1(Vi,q,P)
N = P.Ns(1);
X = reshape(q(P.rows{1}),N,N);  X = (X+X')/2;
Si = Vi'*X*Vi;  Si = (Si+Si')/2;
[W,D] = eig(Si);  [~,p] = sort(diag(D),'descend');  W = W(:,p);
Vs = Vi*W(:,1:end-1);
[Vs,~] = qr(Vs,0);
end
function s = tern(c,a,b)
if c, s = a; else, s = b; end
end
