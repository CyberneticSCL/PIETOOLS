function [prog,meta] = heatNd_lpi(pie,d,k,ep,opts)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% [PROG,META] = HEATND_LPI(PIE,D,K,EP,OPTS) the UNSOLVED stability LPI of
% Jagt & Peet, arXiv:2508.14840v4, Cor. 35, for a PIE from HEATND_PIE, on
% the container path (copvar/cdopvar; any number of spatial variables):
%
%   find P = c' Z_d (indefinite) and R, Q >= 0 with
%   (15a)  P*T - ep^2 T*T - R = 0           (so P*T = T*P = ep^2 T*T + R),
%   (15b)  P*A + A*P + 2k P*T + Q = 0,
%
% whence the PDE is exponentially PIE-to-PDE stable with rate K (Thm. 34,
% which needs ep > 0 and K >= 0; both are checked).
%
% KAPPA. A = A0 + r T (HEATND_PIE), so (15b) is
%   P*A0 + A0*P + 2 kappa P*T + Q = 0,   kappa = r + K,
% and (15a) does not involve r: the program at (r, K) IS the program at
% (0, r + K). It is built that way (from pie.A0), so r only shifts K and
% the base is reusable across r. MEASURED: the former build from pie.A at
% (r, k) gave bitwise the SDP at (0, r + k) in 1-D, 2-D and 3-D (review,
% 09/27/2026), and this build reproduces the former one's dumps bitwise
% (1-D, 2-D, incl. r = 12); TEST_HEATND_LPI (7) checks the literal Cor. 35
% encoding against this one (equal to <= 3e-17 relative, 0 at r = 12).
%
% ENCODING
% - (15a), OPTS.eq15a: 'split' (default) is the Sec. 7.1 listing's form,
%   P*T - T*P = 0 and P*T - ep^2 T*T - R = 0, each on one cell per adjoint
%   pair ('symmetric'; exact for the anti-self-adjoint P*T - T*P as well),
%   which together are (15a). 'full' imposes (15a) on every coefficient of
%   every cell. MEASURED: same row rank, and the stacked system has the same
%   rank, so the same affine set (N = 1, 2 here; N = 3 by review). With the
%   default 'Tspan' bases split has 20-30% fewer rows (1-D DD 64 -> 51, 2-D
%   1712 -> 1242); NOT in general: with 'paper' dp = 0 at d = 1 it has MORE
%   (review: 2-D 277 -> 301, 3-D 3593 -> 4295). Verdicts agree
%   (test_heatNd_lpi);
% - (15b) is X1 + kappa X2 + Q = 0, 'symmetric', X1 = P*A0 + (P*A0)*, X2 =
%   P*T + (P*T)*: both self-adjoint for every decision value, and X2 = 2 P*T
%   on (15a). X1, X2 do not depend on kappa, so a new kappa is one 'plus'
%   and one 'lpi_eq_cdopvar' on a stored base (OPTS.base), and the SDP is
%   affine in kappa (HEATND_BISECT uses that).
%
% BASES (paper eqs. 13-14; mu(d) = d^2+4d+3 per direction)
% - P, OPTS.Pbasis: 'paper' (default) exactly Z_d, per direction
%     multiplier s^a, a <= d                   (d+1),
%     lower/upper s^a th^c, JOINT a + c <= d   ((d+1)(d+2)/2 each),
%   Kronecker over directions, one decision variable each: mu(d)^N.
%   'tensor': 'lpivar_cdopvar' degree d (a, c <= d separately), the stock
%   'lpivar' convention of the Sec. 7.1 listing: ((d+1)+2(d+1)^2)^N.
%   OPTS.Pmult = false drops the multiplier cells (paper Sec. 6.2 remark);
%   MEASURED certified infeasible at 0.5 k* in 2-D, so the default keeps
%   them.
% - R, Q, OPTS.RQ with degree OPTS.dp (scalar or [dR dQ]):
%   'Tspan'   (default) per direction the smallest set {th^a s^c : a <= 1,
%             c <= 1, a + c <= cj} whose span contains T_k (cj = 2 for DD,
%             1 for DN/ND, read off T_k's kernels), each cap raised by dp;
%             T is then in the span of R's basis, so T*T is a Gram form in
%             it. (P = T itself is in the P space only for d >= 2 in a DD
%             direction, d >= 1 in DN/ND: the DD kernel s th has joint
%             degree 2.) This is Z_{d'} of Cor. 35 with the per-direction
%             joint cap adapted to the BC;
%   'paper'   Z_{d'} of Cor. 35 exactly, d' = dp (joint cap a + c <= d'):
%             Gram mu(d')^N before pruning;
%   'tensor'  copquadvar scalar degree dp (a, c <= dp, no joint cap);
%   'balance' the listing's rule: lpi_ineq sizes the positive variable by
%             'degbalance' of the operator cancelled; container analogue
%             'eq_opts_sopvar' + '@sdopvar/degbalance' (dp unused).
%   'custom'  OPTS.Rdeg and OPTS.Qdeg, copquadvar degree specs used as      % MMP, 10/08/2026
%             given (dp unused): R and Q sized from their targets by the    % MMP, 10/08/2026
%             lift/weight rules of sopvar_lift_notes Sec. 7 (HEATND_TAILOR). % MMP, 10/08/2026
% - OPTS.prune (default true): only the basis blocks that can reach the
%   support of the cancelled operator ('eq_opts_sopvar'; lossless by its
%   argument); 'balance' always prunes, as lpi_ineq excludes zero blocks.
% - OPTS.psatz, terms added at the SAME degree (stock eq_deg_psatz = eq_deg;
%   -1 destroys certificates, memory psatz-generators-2d), to R and/or Q
%   (OPTS.psatzon 'RQ' default, or 'Q'):
%   'linear'  (default) the 2N normalised face generators (th_i-a_i)/L_i,
%             (b_i-th_i)/L_i (OPTS.faces, N x 2 logical, default all) via
%             poscopvar psatz = 2i+1 / 2i+2 (copquadvar face codes, MMP
%             09/27/2026; formerly the copy HEATND_POSW, bitwise the same
%             SDPs): stock poslpivar_2d psatz = [3 4 5 6] (CC,
%             09/23/2026) extended to N-D. MEASURED 2-D (review): dropping
%             the Neumann face or the s2 Dirichlet face gives no
%             certificate (MOSEK UNKNOWN, rel_b 3.5e-3 to 7e-3, no valid
%             ray) at k = 0, 0.1, 0.5 k*; psatz on Q only ('psatzon' 'Q')
%             is certified INFEASIBLE there (Farkas rel. violation <= 9e-13),
%             so psatz on R is needed for feasibility;
%   'product' poscopvar psatz = 1, g = prod_i (th_i-a_i)(b_i-th_i): stock
%             lpi_ineq opts.psatz = 1 (complete in 1-D; certified
%             infeasible at 0.5 k* in 2-D here, as for stock);
%   'none'    Cor. 35 as printed (no Psatz term). MEASURED: certifies no
%             rate in 1-D or 2-D at practical degree, as the stock listing.
%   'tensor'  the Markov-Lukacs tensor set: g_i = (th_i-a_i)(b_i-th_i)     % MMP, 10/08/2026
%             alone in direction i and the product of all g_i, each at     % MMP, 10/08/2026
%             'int' lowered by OPTS.psatz_offset in the directions it acts % MMP, 10/08/2026
%             in (poscopvar_direct codes 2N+2+i and 1); in 1-D with offset % MMP, 10/08/2026
%             1 it is the 'product' pair. MEASURED 2-D DD x DN, d = 0,     % MMP, 10/08/2026
%             R (1,1), Q (2,2) uncapped: faces at w-1 certify 12.336731 at % MMP, 10/08/2026
%             nx 339209, the tensor set at w-1 12.326178 at nx 395785 in   % MMP, 10/08/2026
%             twice the time, the tensor set at w nx 762889; it loses to   % MMP, 10/08/2026
%             'linear' with psatz_offset 1 (README Sec. 12).               % MMP, 10/08/2026
%   The 2N face terms are summed by one 'plus_batch' (10/01/2026); summed
%   one 'plus' at a time, each re-merged the decision list: MEASURED 6.0 s
%   of the 56.5 s Q stage in 3-D (review).
% - OPTS.psatz_offset (default 0): the Psatz terms at 'int' reduced by this, % MMP, 10/08/2026
%   floored at 0; a product term that would need 'int' < 0 is omitted. 1    % MMP, 10/08/2026
%   with 'product' is the Markov-Lukacs pair S0 + g S1, deg S1 = deg S0 - 2 % MMP, 10/08/2026
%   (sopvar_lift_notes Sec. 3). Only for 'int'/'mult' specs, not 'subset'.  % MMP, 10/08/2026
% - OPTS.preset sets defaults, individual fields still override:
%   'bench' (default) Tspan dp 0, linear;  'heavy' Tspan dp 1, linear;
%   'listing' balance, none (the Sec. 7.1 listing as printed);
%   'listingL' balance, linear (the listing with the linear generators).
%
% INPUT
% - pie:  struct from HEATND_PIE (needs pie.A0);   d: degree of P;
% - k:    decay rate in (15b), k >= 0 (Thm. 34); the program is built at
%         kappa = pie.r + k; [] returns the program with (15a) only;
% - ep:   epsilon of (15a), ep > 0 (Thm. 34; paper: 0.1);
% - opts: fields above, plus 'base' (a META.base of an earlier call with
%         the same pie.T/A0, d, ep and options, any r: only (15b) at the
%         new kappa is added).
%
% OUTPUT
% - prog: SOSTOOLS/LPI program with 'eq' expressions only (no 'ineq'), so
%         the pre-solve SDP is complete (HEATND_SDP);
% - meta: N, r, d, k, kappa (= r + k, the only rate the SDP sees), ep, opts
%         (resolved), kstar, lambda1 (= kappa*); t (stage wall times, s);
%         mem (MB: MATLAB in use after the build, process PEAK working set
%         and PEAK private bytes - lifetime peaks, per build only in a fresh
%         process or in increasing size order - and current private bytes);
%         nP (decision variables of P), ndec (numel(prog.decvartable));
%         sdp (m rows, Kf, Ks Gram sizes, nnzAt, nx); base.
%
% Cost: dominated by the R and Q declarations (copquadvar / int_semisep)
% and the compositions P*T, P*A, whose rows scale with the decision count;
% nothing is dense in the decision count (CLAUDE.md sec. 2). MEASURED 3-D
% (bench, d = 0): base 55-63 s plus ~10 s per finalization, 3.8 GB, 2.7M
% decision variables, m = 31770.
%
% Initial coding MMP, 09/27/2026
% MMP, 09/27/2026 (review fixes): build (15b) at kappa = r + k from pie.A0
%   (r only shifts k); reject ep <= 0 and k < 0 (Thm. 34); qualify the
%   split row-count and 'P = T' header claims; correct the ablation wording.
% MMP, 09/27/2026 (final reviews): comment only - the RQ_deg 'Tspan' note
%   no longer implies P = T is admissible at every d (as the header).
% MMP, 09/27/2026 (library face codes): 'linear' faces now come from
%   poscopvar psatz = 2i+1 / 2i+2, the option copquadvar gained today, not
%   from the copy heatNd_posw (retired). MEASURED: every SDP of the
%   switch-over set (N = 1, 2, 3-D bench d = 0) bit-identical before/after.
% MMP, 10/01/2026: The program is lpiprogram for every N; it no longer refuses
%   N > 2, so the hand-built copy for N > 2 is commented out. Same program
%   (lpiprogram builds exactly what the copy built).
% MMP, 10/01/2026: The 'linear' face terms are summed by one plus_batch
%   instead of Z = Z + Zg per face: each '+' remapped the growing sum
%   onto the union of the decision lists (O(N^2) block rows in the number
%   of terms). Same program; see @cdopvar/plus_batch.
% MMP, 10/08/2026: OPTS.psatz 'tensor' (the Markov-Lukacs tensor set through
%   poscopvar_direct's single-direction quadratic codes), to test it against
%   the faces on the 2-D benchmark before it becomes the N-D default.
% MMP, 10/08/2026: OPTS.RQ 'custom' (Rdeg, Qdeg given) and OPTS.psatz_offset,
%   so that HEATND_TAILOR can size R and Q from their targets by the
%   lift/weight rules and put the product term one degree below the plain
%   term. Defaults reproduce every earlier program: 'custom' is never a
%   preset and offset 0 is the former behaviour. Both enter the base key.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if nargin<5 || isempty(opts),   opts = struct();    end
if ~(isscalar(ep) && ep>0)                  % Thm. 34 / Cor. 35: epsilon > 0
    error('heatNd_lpi:ep','EP must be > 0 (Thm. 34); ep = 0 makes b = 0.')
end
if ~isempty(k) && ~(isscalar(k) && k>=0)    % Thm. 34 / Cor. 35: k >= 0
    error('heatNd_lpi:k','K must be >= 0 (Thm. 34); for rate kappa use pie with r = 0.')
end
if ~isfield(pie,'A0'),  error('heatNd_lpi:A0','PIE has no A0; rebuild it with HEATND_PIE.'),  end
if ~isfield(opts,'preset') || isempty(opts.preset),  opts.preset = 'bench';  end
switch opts.preset
    case 'bench',       pr = struct('RQ','Tspan','dp',0,'psatz','linear');
    case 'heavy',       pr = struct('RQ','Tspan','dp',1,'psatz','linear');
    case 'listing',     pr = struct('RQ','balance','dp',0,'psatz','none');
    case 'listingL',    pr = struct('RQ','balance','dp',0,'psatz','linear');
    otherwise,          error('heatNd_lpi:preset','Unknown preset ''%s''.',opts.preset)
end
df = struct('Pbasis','paper','Pmult',true,'RQ',pr.RQ,'dp',pr.dp,'psatz',pr.psatz, ...
            'psatzon','RQ','faces',[],'prune',true,'eq15a','split','base',[], ...
            'Rdeg',[],'Qdeg',[],'psatz_offset',0);                          % MMP, 10/08/2026
fn = fieldnames(df);
for i = 1:numel(fn),    if ~isfield(opts,fn{i}),    opts.(fn{i}) = df.(fn{i});  end,    end
dp = opts.dp;   if isscalar(dp),    dp = [dp dp];   end
N = pie.N;  vars = pie.vars;    dom = pie.dom;
t = struct('ops',0,'P',0,'compose',0,'R',0,'Q',0,'eq15a',0,'eq15b',0);

if isempty(opts.base)
    % % % Program in N variables. lpiprogram refuses N > 2 (lpiprogram.m,
    % 'more than 2 spatial variables'); for N > 2 build what it would return.
    % (10/01/2026: it no longer refuses; one call for every N.)             % MMP, 10/01/2026
    t0 = tic;
%   if N<=2                                                                 % MMP, 10/01/2026 (was)
        prog = lpiprogram(polynomial(vars(:)),[],dom);
%   else                                                                    % MMP, 10/01/2026 (was)
%       prog = sosprogram(polynomial([]),dpvar(zeros(0,1)));                % MMP, 10/01/2026 (was)
%       prog.vartable = [prog.vartable; polynomial(vars(:)); polynomial(strcat(vars(:),'_dum'))]; % MMP, 10/01/2026 (was)
%       prog.dom = dom;                                                     % MMP, 10/01/2026 (was)
%   end                                                                     % MMP, 10/01/2026 (was)
    T = pie.T;  A = pie.A0; t.ops = toc(t0);    % r enters via kappa only
    % % % P = c' Z_d.
    t0 = tic;
    switch opts.Pbasis
        case 'paper',   [prog,P] = lpivar_Zd(prog,vars,dom,d,opts.Pmult);
        case 'tensor'
            if ~opts.Pmult, error('heatNd_lpi:Pmult','Pmult = false needs Pbasis ''paper''.'),  end
            [prog,P] = lpivar_cdopvar(prog,1,{vars},dom,d);
        otherwise,      error('heatNd_lpi:Pbasis','Pbasis is ''paper'' or ''tensor''.')
    end
    nP = numel(prog.decvartable);   t.P = toc(t0);
    % % % Compositions, once.
    t0 = tic;
    PT = P'*T;  PA = P'*A;  TT = T'*T;
    X1 = PA + PA';  X2 = PT + PT';
    E1 = PT - ep^2*TT;                          % (15a) is E1 - R = 0
    t.compose = toc(t0);
    % % % R >= 0 sized on E1, Q >= 0 on X1 + X2 (support of X1 + k X2, k ~= 0).
%   dR = RQ_deg(pie,opts.RQ,dp(1));     dQ = RQ_deg(pie,opts.RQ,dp(2));     % MMP, 10/08/2026 (was)
    if strcmp(opts.RQ,'custom')     % specs given, sized from the targets   % MMP, 10/08/2026
        if isempty(opts.Rdeg) || isempty(opts.Qdeg)                         % MMP, 10/08/2026
            error('heatNd_lpi:custom','RQ ''custom'' needs OPTS.Rdeg and OPTS.Qdeg.') % MMP, 10/08/2026
        end                                                                 % MMP, 10/08/2026
        dR = opts.Rdeg;     dQ = opts.Qdeg;                                 % MMP, 10/08/2026
    else                                                                    % MMP, 10/08/2026
        dR = RQ_deg(pie,opts.RQ,dp(1));     dQ = RQ_deg(pie,opts.RQ,dp(2)); % MMP, 10/08/2026
    end                                                                     % MMP, 10/08/2026
    t0 = tic;   [prog,R] = posvar(prog,E1,vars,dom,opts,dR,contains(opts.psatzon,'R'));    t.R = toc(t0);
    t0 = tic;   [prog,Q] = posvar(prog,X1+X2,vars,dom,opts,dQ,contains(opts.psatzon,'Q')); t.Q = toc(t0);
    % % % (15a), in the form OPTS.eq15a.
    t0 = tic;
    switch opts.eq15a
        case 'full'     % every coefficient of every cell
            prog = lpi_eq_cdopvar(prog,E1 - R);
        case 'split'    % Sec. 7.1 listing form: P*T = T*P, then the
            % self-adjoint part; one cell per adjoint pair each ('symmetric'
            % is exact for the anti-self-adjoint PT - PT' too). Same affine
            % set as 'full', about half the rows.
            prog = lpi_eq_cdopvar(prog,PT - PT','symmetric');
            prog = lpi_eq_cdopvar(prog,E1 - R,'symmetric');
        otherwise
            error('heatNd_lpi:eq15a','eq15a is ''split'' or ''full''.')
    end
    t.eq15a = toc(t0);
    % Key without r: the base is r-independent (built from A0).
    base = struct('prog',prog,'X1',X1,'X2',X2,'Q',Q,'P',P,'R',R,'E1',E1,'nP',nP,'t',t, ...
                  'key',{{N,pie.bc,pie.dom,d,ep,opts.Pbasis,opts.Pmult,opts.RQ,dp,opts.psatz,opts.psatzon,opts.faces,opts.prune,opts.eq15a, ...
                          opts.Rdeg,opts.Qdeg,opts.psatz_offset}});         % MMP, 10/08/2026
else
    base = opts.base;   prog = base.prog;   nP = base.nP;   t = base.t;
    if ~isequal(base.key,{N,pie.bc,pie.dom,d,ep,opts.Pbasis,opts.Pmult,opts.RQ,dp,opts.psatz,opts.psatzon,opts.faces,opts.prune,opts.eq15a, ...
                          opts.Rdeg,opts.Qdeg,opts.psatz_offset})           % MMP, 10/08/2026
        error('heatNd_lpi:base','OPTS.base was built for different inputs.')
    end
end
% % % (15b) at kappa = r + k (X1 from A0).
kappa = [];
if ~isempty(k)
    t0 = tic;
    kappa = pie.r + k;
    prog = lpi_eq_cdopvar(prog,base.X1 + kappa*base.X2 + base.Q,'symmetric');
    t.eq15b = toc(t0);
end

opts.base = [];
meta = struct('N',N,'r',pie.r,'d',d,'k',k,'kappa',kappa,'ep',ep,'opts',opts,'kstar',pie.exact.kstar, ...
              'lambda1',pie.exact.lambda1,'t',t,'mem',memnow(),'nP',nP,'ndec',numel(prog.decvartable));
meta.sdp = shape(prog);
meta.base = base;
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [prog,Z] = posvar(prog,X,vars,dom,opts,deg,usepsatz)
% Z = Z* >= 0 on L2[vars] to cancel the 1x1 container X: base Gram term
% plus (if USEPSATZ) the psatz terms of OPTS.psatz, at one degree and basis.
% DEG: copquadvar degree spec, or [] for 'balance' (sized on X here).
N = numel(vars);
po = struct();
if strcmp(opts.RQ,'balance')
    [eo,~] = eq_opts_sopvar(X.C{1,1});
    deg = {degbalance(X.C{1,1},eo)};            % one struct per kept block
    po.include = {eo.include};
    if isfield(eo,'sep'),   po.sep = eo.sep;    end
elseif opts.prune
    [eo,~] = eq_opts_sopvar(X.C{1,1});
    po.include = {eo.include};
    if isfield(eo,'sep'),   po.sep = eo.sep;    end
end
[prog,Z] = poscopvar(prog,1,vars,dom,deg,po);
if ~usepsatz,   return,     end
% The Psatz terms at 'int' reduced by OPTS.psatz_offset (0: as before).     % MMP, 10/08/2026
% Empty when a product term would need 'int' < 0: omitted, as the           % MMP, 10/08/2026
% Markov-Lukacs pair has no S1 at w = 0.                                    % MMP, 10/08/2026
degp = psatz_deg(deg,opts.psatz_offset);                                    % MMP, 10/08/2026
switch opts.psatz
    case 'none'
    case 'product'
        if isempty(degp),   return,     end                                 % MMP, 10/08/2026
        pp = po;    pp.psatz = 1;
%       [prog,Zg] = poscopvar(prog,1,vars,dom,deg,pp);                      % MMP, 10/08/2026 (was)
        [prog,Zg] = poscopvar(prog,1,vars,dom,degp,pp);                     % MMP, 10/08/2026
        Z = Z + Zg;
    case 'linear'
        if isempty(degp),   return,     end     % no face term below w = 0  % MMP, 10/08/2026
        % OPTS.faces: N x 2 logical, face (i,1) at a_i, (i,2) at b_i.
        fc = opts.faces;    if isempty(fc),     fc = true(N,2);     end
        % Face (i,e): poscopvar psatz 2i+1 (e=1) / 2i+2 (e=2), weight       % MMP, 09/27/2026
        % (th_i-a_i)/L_i / (b_i-th_i)/L_i at theta; i indexes the SORTED    % MMP, 09/27/2026
        % registry, which vars = s1..sN is (heatNd_pie: N <= 9), as dom     % MMP, 09/27/2026
        % row order assumes too.                                            % MMP, 09/27/2026
        % The terms are summed once at the end ('plus_batch'): each '+'     % MMP, 10/01/2026
        % remapped the growing sum onto the union of the decision lists.    % MMP, 10/01/2026
        terms = {Z};                                                        % MMP, 10/01/2026
        for i = 1:N
%           si = polynomial(vars(i));   L = dom(i,2)-dom(i,1);              % MMP, 09/27/2026 (was)
%           gi = {(si-dom(i,1))/L, (dom(i,2)-si)/L};                        % MMP, 09/27/2026 (was)
            for e = find(fc(i,:))
%               pg = po;    pg.gfun = gi{e};                                % MMP, 09/27/2026 (was)
%               [prog,Zg] = heatNd_posw(prog,1,vars,dom,deg,pg);            % MMP, 09/27/2026 (was)
                pg = po;    pg.psatz = 2*i+e;                               % MMP, 09/27/2026
%               [prog,Zg] = poscopvar(prog,1,vars,dom,deg,pg);              % MMP, 09/27/2026 % MMP, 10/08/2026 (was)
                [prog,Zg] = poscopvar(prog,1,vars,dom,degp,pg);             % MMP, 10/08/2026
%               Z = Z + Zg;                                                 % MMP, 10/01/2026 (was)
                terms{end+1} = Zg;                                          % MMP, 10/01/2026
            end
        end
        Z = plus_batch(terms{:});                                           % MMP, 10/01/2026
    case 'tensor'                                                           % MMP, 10/08/2026
        % Markov-Lukacs tensor set (sopvar_lift_notes Sec. 9): g_i =        % MMP, 10/08/2026
        % (th_i-a_i)(b_i-th_i) alone in direction i at 'int' lowered by the % MMP, 10/08/2026
        % offset in that direction only, and the product of all g_i at      % MMP, 10/08/2026
        % 'int' lowered in every direction: poscopvar_direct codes 2N+2+i   % MMP, 10/08/2026
        % and 1 in one call. A 'subset' spec takes offset 0 only.           % MMP, 10/08/2026
        pt = po;    pt.psatz = [2*N+2+(1:N), 1];                            % MMP, 10/08/2026
        pt.psatz_offset = opts.psatz_offset*ones(1,N+1);                    % MMP, 10/08/2026
        [prog,Zg] = poscopvar_direct(prog,1,vars,dom,deg,pt);               % MMP, 10/08/2026
        Z = Z + Zg;                                                         % MMP, 10/08/2026
    otherwise
%       error('heatNd_lpi:psatz','psatz is ''none'', ''product'' or ''linear''.') % MMP, 10/08/2026 (was)
        error('heatNd_lpi:psatz','psatz is ''none'', ''product'', ''linear'' or ''tensor''.') % MMP, 10/08/2026
end
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function degp = psatz_deg(deg,off)                                          % MMP, 10/08/2026
% The degree spec of a Psatz term: DEG with 'int' lowered by OFF, floored   % MMP, 10/08/2026
% at 0; [] when some 'int' would go below 0 (no such term). A 'subset'      % MMP, 10/08/2026
% spec (Tspan) is refused, since its caps would have to move with 'int'.    % MMP, 10/08/2026
degp = deg;                                                                 % MMP, 10/08/2026
if isempty(off) || off==0,  return,     end                                 % MMP, 10/08/2026
if iscell(deg)                                                              % MMP, 10/08/2026
    for i = 1:numel(deg)                                                    % MMP, 10/08/2026
        degp{i} = psatz_deg(deg{i},off);                                    % MMP, 10/08/2026
        if isempty(degp{i}),    degp = [];  return,     end                 % MMP, 10/08/2026
    end                                                                     % MMP, 10/08/2026
    return                                                                  % MMP, 10/08/2026
end                                                                         % MMP, 10/08/2026
if isnumeric(deg),  deg = struct('int',deg);    degp = deg;     end         % MMP, 10/08/2026
if isfield(deg,'subset') && ~isempty(deg.subset)                            % MMP, 10/08/2026
    error('heatNd_lpi:psatz_offset','psatz_offset is not defined for a ''subset'' spec.') % MMP, 10/08/2026
end                                                                         % MMP, 10/08/2026
w = deg.int;    if ~isfield(deg,'int') || isempty(w),   w = 1;  end         % MMP, 10/08/2026
if any(w-off<0),    degp = [];  return,     end                             % MMP, 10/08/2026
degp.int = w-off;                                                           % MMP, 10/08/2026
if ~isfield(degp,'mult') || isempty(degp.mult),    degp.mult = w;   end     % MMP, 10/08/2026
end                                                                         % MMP, 10/08/2026


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function deg = RQ_deg(pie,RQ,dp)
% copquadvar degree spec of the R, Q basis on L2[s1..sN] ([] for 'balance',
% sized per operator in POSVAR). copquadvar builds the exponent grid over
% [th_1..th_N, s_1..s_N] ('int' caps the th_k, 'mult' the s_k, s_k zeroed in
% a multiplier direction), so a Kronecker product of per-direction sets
%   {th^a s^c : a <= ct_k, c <= cs_k, a + c <= cj_k}
% is exactly the 'subset' cap: subset b gets sum_k of the cap of b's part
% in direction k (0, ct_k, cs_k or cj_k) - every multi-direction subset cap
% is then a sum of per-direction ones and removes nothing further.
%   'paper'   ct = cs = cj = d'  (Z_{d'}, eq. 13, mu(d')^N before pruning);
%   'tensor'  scalar d' (a, c <= d', no joint cap);
%   'Tspan'   the smallest such set whose span contains T_k (kernel
%             monomials of T_k's lower/upper cells: DD th s -> cj = 2, DN/ND
%             -> cj = 1; ct = cs = 1), raised by d' in every cap. T then
%             lies in the span of R's basis, so T*T is a Gram form there.
%             P = T itself is admissible only for d >= 2 in a DD direction,
%             d >= 1 in DN/ND (the DD kernel s th has joint degree 2).
N = pie.N;
switch RQ
    case 'balance',     deg = [];   return
    case 'tensor',      deg = dp;   return
    case 'paper',       ct = dp*ones(1,N);  cs = ct;    cj = ct;
    case 'Tspan'
        ct = zeros(1,N);    cs = ct;    cj = ct;
        for k = 1:N
            [a,c] = find(pie.op1(k).C{2} | pie.op1(k).C{3});   % C(a,c): s_out^(a-1) s_in^(c-1)
            ct(k) = max(a)-1;   cs(k) = max(c)-1;   cj(k) = max(a+c)-2;
        end
        ct = ct+dp;     cs = cs+dp;     cj = cj+dp;
    otherwise
        error('heatNd_lpi:RQ','RQ is ''balance'', ''paper'', ''tensor'' or ''Tspan''.')
end
n = 2*N;    sc = zeros(1,2^n);
for i = 1:2^n
    b = bitget(i-1,1:n)>0;      bt = b(1:N);    bs = b(N+1:n);
    sc(i) = sum(ct.*(bt & ~bs) + cs.*(bs & ~bt) + cj.*(bt & bs));
end
deg = struct('int',ct,'mult',cs,'subset',sc);
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [prog,Pop] = lpivar_Zd(prog,vars,dom,d,Pmult)
% P = c' Z_d with the paper's joint-degree basis (eq. 13-14); the layout of
% 'lpivar_cdopvar' (ZL = ZR = {(0:d)'} per variable, gamma over sorted
% vars with direction 1 fastest, coefficient index row + NL*(col-1), ZL/ZR
% Kronecker first variable slowest), with the joint cap a + c <= d in an
% integral direction and right degree 0 in a multiplier direction
% (canonical form, sopvar.m), so the constructor folds nothing.
% One triplet per decision variable, one sparse per gamma cell.
N = numel(vars);
ZL = repmat({(0:d)'},1,N);      ZR = ZL;
n1 = d+1;   NL = n1^N;  nC = NL*NL;
st = n1.^(N-1:-1:0);                        % strides, first variable slowest
[a1,c1] = ndgrid(0:d,0:d);
okI = a1(:)+c1(:)<=d;
pairI = [a1(okI),c1(okI)];      pairM = [(0:d)',zeros(n1,1)];
pos = cell([3*ones(1,N),1]);
for g = 1:3^N
    c = cell(1,N);  [c{:}] = ind2sub([3*ones(1,N),1],g);    gam = [c{:}];
    L = 0;  Rr = 0;                         % 0-based row / column indices
    for kk = 1:N
        if gam(kk)==1,  pr = pairM;     else,   pr = pairI;     end
        L  = reshape(L(:)  + st(kk)*pr(:,1)',[],1);
        Rr = reshape(Rr(:) + st(kk)*pr(:,2)',[],1);
    end
    pos{g} = L + NL*Rr + 1;
    if ~Pmult && any(gam==1),   pos{g} = zeros(0,1);    end     % integral-only P
end
nvar = sum(cellfun(@numel,pos(:)));
[prog,dv] = lpidecvar(prog,[nvar,1]);
Zd = reshape(dv.dvarname,[],1);
A = cell(size(pos));    B = cell(size(pos));    k0 = 1;
for g = 1:numel(pos)
    ng = numel(pos{g});
    A{g} = sparse(nC,1);
    B{g} = sparse(k0+(0:ng-1)',pos{g},ones(ng,1),nvar,nC);
    k0 = k0+ng;
end
dm = struct('in',dom,'out',dom);
blk = sdopvar(struct('A',{A},'B',{B}),struct('in',{vars},'out',{vars}),Zd,ZL,ZR,dm,[1 1]);
Pop = cdopvar({blk});
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function s = shape(prog)
% Pre-solve SDP shape from the program itself (what processvars in
% sossolve.m builds for an eq-only program): rows, free and PSD columns,
% nonzeros. O(number of expressions).
s = struct('m',0,'Kf',prog.var.idx{1}-1,'Ks',[],'nnzAt',0,'nx',0);
for i = 1:prog.expr.num
    s.m = s.m + numel(prog.expr.b{i});  s.nnzAt = s.nnzAt + nnz(prog.expr.At{i});
end
for i = 1:prog.var.num
    sz = prog.var.idx{i+1}-prog.var.idx{i};
    if strcmp(prog.var.type{i},'sos'),  s.Ks(end+1) = round(sqrt(sz));
    else,                               s.Kf = s.Kf + sz;
    end
end
s.nx = s.Kf + sum(s.Ks.^2);
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function m = memnow()
% MB: MATLAB's own use (memory), and the process peaks from .NET (Windows).
m = struct('matlab',NaN,'peak_ws',NaN,'peak_private',NaN,'private',NaN);
try,    u = memory;     m.matlab = u.MemUsedMATLAB/2^20;    catch,  end
try
    p = System.Diagnostics.Process.GetCurrentProcess();     p.Refresh();
    m.peak_ws = double(p.PeakWorkingSet64)/2^20;
    m.peak_private = double(p.PeakPagedMemorySize64)/2^20;
    m.private = double(p.PrivateMemorySize64)/2^20;
catch
end
end
