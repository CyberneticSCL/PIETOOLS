function v = heatNd_solve(in,opts)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% V = HEATND_SOLVE(IN,OPTS) solves a heatNd stability SDP with MOSEK and
% returns a CERTIFIED verdict - never the solver status alone:
%
%   st = +1  prosta PRIMAL_AND_DUAL_FEASIBLE, solsta OPTIMAL, AND the
%            returned primal satisfies rel_b = ||At'x - b||/||b|| <= btol
%            (1e-6) AND every PSD block has min eig >= -psdtol (1e-8) times
%            the largest eigenvalue over all blocks, AND x ~= 0;
%   st = -1  solsta PRIMAL_INFEASIBLE_CER AND the Farkas ray is verified:
%            scaled to b'y = -1 (||b|| = 1, HEATND_SDP), the ABSOLUTE
%            violation of At*y in K*, cert_viol, is <= ctol (1e-8). Then no
%            feasible x has ||x_f|| + sum_j tr(X_j) < cert_radius =
%            1/max(||w_f||, max_j neg_j) (w = At*y, neg_j = the negative
%            part of W_j's least eigenvalue): -1 = <w,x> >= -||w_f|| ||x_f||
%            - sum_j neg_j tr(X_j). MEASURED (review, 09/27/2026): recorded
%            rays have cert_viol <= 1.9e-9, radius >= 6.6e8, against tr(X)
%            <= 1.6e4 at certified +1 points. A relative test (viol over the
%            ray's own slack scale, used before) passes a crafted ray of
%            violation 1 on a feasible SDP. An unverified ray is 0. Or the
%            program is DEFICIENT: an equality
%            row with no variable and b_i ~= 0 (sossolve's presolve rejects
%            it before any solver runs; 0 = b_i is itself a certificate, but
%            it says the DEGREE is too low, not that k is too large);
%   st =  0  anything else ('uncertain': reported, never coerced).
% MEASURED: the 1e-6 gate sits at MOSEK's default accuracy floor on these
% SDPs, so near-threshold verdicts flip under 1e-15 changes of the data
% (review, 09/27/2026); HEATND_BISECT retries uncertain verdicts. With
% OPTS.tight, 2-D rel_b is 1e-10 to 1.2e-8 (MEASURED: 16 solves at 8
% points near kappa-hat, bench and heavy, all rows and independent rows)
% for about 5 more iterations (+35% time). The +1 gate can accept an
% infeasible SDP within ~1e-8 of the threshold (MEASURED, review: 2-D heavy
% d = 0 +1 at lambda1 + 1e-8, rel_b 4.1e-7); HEATND_BISECT rejects any +1
% above lambda1.
%
% IN is an unsolved program (HEATND_LPI), a dump struct (HEATND_SDP), or
% the path of a dump file saved by HEATND_SDP. For a FILE, the route that
% produced the dump's reference verdict (info.ref.rows / info.ref.tight)
% and its stored row set (info.keep) are the defaults, so heatNd_solve(file)
% repeats the documented route; OPTS.rows / OPTS.tight override.
% OPTS (defaults): btol 1e-6, psdtol 1e-8, ctol 1e-8, route 'dump'
%   (Sedumi2Mosek + mosekopt on exactly the SDP sossolve would build:
%   HEATND_SDP) or 'lpisolve' (sossolve, Mosek; verdict from its
%   pinf/numerr + cx_resid: for cross-checking the dump route), keepx
%   (false), param (mosekopt parameter struct, default none: Mosek
%   defaults, as sossolve), tight (false: true sets MSK_DPAR_INTPNT_CO_TOL_
%   PFEAS/DFEAS/REL_GAP = 1e-10 on top of param), keep (row indices from
%   HEATND_LINDEP: Mosek sees only these rows - a presolve of linearly
%   dependent equalities; rel_b is still measured on ALL rows), rows
%   ('keep' uses OPTS.keep or the file's info.keep, 'all' ignores them;
%   default: 'keep' if OPTS.keep is given, else the file's reference route).
%
% OUTPUT struct V: st, why (text), prosta, solsta, rel_b, eta (normwise
%   backward error ||r||_inf / (max_i |At_i|'|x| + ||b||_inf), diagnostic only),
%   psd_min, psd_relmin, trivial, deficient, cert_viol (Farkas ray scaled to b'y =
%   -1: max violation of At*y in K*, the GATED quantity; NaN unless
%   infeasible), cert_radius (above), cert_rel (cert_viol / max(1, max
%   |slack|), diagnostic only), rows ('all' |
%   'keep'), nrows (rows given to Mosek), tight, t_mosek
%   (MSK_DINF_OPTIMIZER_TIME), iter, wall (s, including conversion), m, nx,
%   Ks, threads (MATLAB maxNumCompThreads), mosek_threads (Mosek's own count,
%   MSK_IINF_INTPNT_NUM_THREADS: Mosek does NOT inherit maxNumCompThreads -
%   measured 24 with MATLAB at 8), mem_peak (MB,
%   process peak private bytes - lifetime peak), x (if keepx).
%
% Initial coding MMP, 09/27/2026
% MMP, 09/27/2026 (review fixes): gate -1 on the measured Farkas violation
%   (ctol); 'tight' tolerances; a dump file's stored rows and reference
%   route are the defaults, so the documented route is reproducible.
% MMP, 09/27/2026 (final reviews): -1 gated on the ABSOLUTE Farkas
%   violation (cert_viol; the relative gate of the previous entry was not
%   sound), cert_radius reported; a given OPTS.keep now takes precedence
%   over a file's reference row set, as the header states.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if nargin<2 || isempty(opts),   opts = struct();    end
given = opts;                                   % caller's own fields
df = struct('btol',1e-6,'psdtol',1e-8,'ctol',1e-8,'route','dump','keepx',false,'param',struct(), ...
            'tight',false,'keep',[],'rows','');
fn = fieldnames(df);
for i = 1:numel(fn),    if ~isfield(opts,fn{i}),    opts.(fn{i}) = df.(fn{i});  end,    end
v = struct('st',0,'why','','prosta','','solsta','','rel_b',NaN,'eta',NaN,'psd_min',NaN, ...
           'psd_relmin',NaN,'trivial',true,'deficient',false,'cert_viol',NaN,'cert_radius',NaN,'cert_rel',NaN, ...
           'rows','all','nrows',NaN,'tight',false, ...
           't_mosek',NaN,'iter',NaN,'wall',NaN,'m',NaN,'nx',NaN,'Ks',[], ...
           'threads',maxNumCompThreads,'mosek_threads',NaN,'mem_peak',NaN,'x',[]);
tw = tic;
if strcmp(opts.route,'lpisolve')
    v = via_lpisolve(in,opts,v);
    v.wall = toc(tw);   v.mem_peak = peakMB();
    return
end
if ischar(in) || isstring(in)                   % a HEATND_SDP dump file
    S = load(in);   D = S;  D.m = size(S.At,2);     D.nx = size(S.At,1);
    D.deficient = S.info.deficient;
    % The file's documented route is the default (dumps before 09/27/2026
    % review fixes carry neither field: all rows, Mosek defaults).
    if isempty(opts.keep) && isfield(S.info,'keep'),    opts.keep = S.info.keep;    end
    if isfield(S.info,'ref') && isstruct(S.info.ref)
        gk = isfield(given,'keep') && ~isempty(given.keep);     % a given row set wins (header)
        if ~isfield(given,'rows') && ~gk && isfield(S.info.ref,'rows'),  opts.rows = S.info.ref.rows;  end
        if ~isfield(given,'tight') && isfield(S.info.ref,'tight'),  opts.tight = S.info.ref.tight;  end
    end
elseif isstruct(in) && isfield(in,'At'),    D = in;
else,                                       D = heatNd_sdp(in);
end
if isempty(opts.rows),  if isempty(opts.keep),  opts.rows = 'all';  else,   opts.rows = 'keep';     end,    end
switch opts.rows
    case 'all',     opts.keep = [];
    case 'keep',    if isempty(opts.keep),  error('heatNd_solve:keep','rows ''keep'' but no row set (OPTS.keep / info.keep).'),  end
    otherwise,      error('heatNd_solve:rows','OPTS.rows is ''keep'' or ''all''.')
end
if opts.tight                                   % tighter IPM stopping tolerances
    opts.param.MSK_DPAR_INTPNT_CO_TOL_PFEAS = 1e-10;
    opts.param.MSK_DPAR_INTPNT_CO_TOL_DFEAS = 1e-10;
    opts.param.MSK_DPAR_INTPNT_CO_TOL_REL_GAP = 1e-10;
end
v.rows = opts.rows;     v.tight = opts.tight;
v.m = D.m;  v.nx = D.nx;    v.Ks = D.K.s;   v.deficient = D.deficient;
if D.deficient
    v.st = -1;  v.why = 'deficient: 0 = b_i ~= 0 (sossolve presolve; degree too low)';
    v.wall = toc(tw);   v.mem_peak = peakMB();
    return
end
if isempty(opts.keep)
    prob = Sedumi2Mosek(D.At',full(D.b),full(D.c),D.K);        % as sossolve.m
    v.nrows = D.m;
else                            % independent rows only; verdict on ALL rows
    prob = Sedumi2Mosek(D.At(:,opts.keep)',full(D.b(opts.keep)),full(D.c),D.K);
    v.nrows = numel(opts.keep);
end
[~,res] = mosekopt('minimize info echo(0)',prob,opts.param);
v.prosta = res.sol.itr.prosta;  v.solsta = res.sol.itr.solsta;
v.t_mosek = res.info.MSK_DINF_OPTIMIZER_TIME;   v.iter = res.info.MSK_IINF_INTPNT_ITER;
if isfield(res.info,'MSK_IINF_INTPNT_NUM_THREADS'),  v.mosek_threads = res.info.MSK_IINF_INTPNT_NUM_THREADS;  end
% Primal in SeDuMi cone order: free part, then each block from its lower
% triangle (column-major, as MosekSol2SedumiSol).
Kf = D.K.f;     x = zeros(D.nx,1);  x(1:Kf) = res.sol.itr.xx(1:Kf);
off = Kf;   idx = 0;    emin = inf;     emax = -inf;
for j = 1:numel(D.K.s)
    s = D.K.s(j);   I = find(tril(true(s)));
    X = zeros(s);   X(I) = res.sol.itr.barx(idx+(1:numel(I)));   idx = idx+numel(I);
    X = X + X' - diag(diag(X));
    x(off+(1:s^2)) = X(:);  off = off+s^2;
    ev = eig((X+X')/2);     emin = min(emin,ev(1));     emax = max(emax,ev(end));
end
rr = D.At'*x - D.b;
v.rel_b = norm(rr)/max(norm(D.b),eps);
% Normwise backward error (diagnostic only, never in the verdict). The
% row-wise form is useless here: rows that vanish at x give 0/0 noise.
v.eta = full(norm(rr,inf)/(max(abs(D.At)'*abs(x)) + norm(D.b,inf)));
v.psd_min = emin;   v.psd_relmin = emin/max(emax,eps);
v.trivial = norm(x)<=1e-12 || abs(v.rel_b-1)<=1e-6;
if strcmp(v.solsta,'PRIMAL_INFEASIBLE_CER')
    [v.cert_viol,v.cert_rel,v.cert_radius] = farkas_viol(D,res,opts.keep);
    if v.cert_viol<=opts.ctol                   % absolute, at b'y = -1 (header)
        v.st = -1;  v.why = sprintf('Mosek infeasibility certificate (ray violation %.1e, radius %.1e)', ...
                                    v.cert_viol,v.cert_radius);
    else
        v.why = sprintf('infeasibility certificate not verified (ray violation %.1e > %.0e, radius %.1e)', ...
                        v.cert_viol,opts.ctol,v.cert_radius);
    end
elseif strcmp(v.prosta,'PRIMAL_AND_DUAL_FEASIBLE') && strcmp(v.solsta,'OPTIMAL')
    if v.rel_b<=opts.btol && v.psd_relmin>=-opts.psdtol && ~v.trivial
        v.st = +1;  v.why = 'optimal, residuals certified';
    else
        v.why = sprintf('optimal but not certified (rel_b %.1e, psd_relmin %.1e, trivial %d)', ...
                        v.rel_b,v.psd_relmin,v.trivial);
    end
else
    v.why = sprintf('uncertain: %s / %s',v.prosta,v.solsta);
end
if opts.keepx,  v.x = x;    end
v.wall = toc(tw);   v.mem_peak = peakMB();
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [viol,rel,rad] = farkas_viol(D,res,keep)
% Farkas: x in K, At'x = b is infeasible if some y has At*y in K* (free
% part 0, blocks PSD) and b'y < 0. Mosek's y is signed per its own form;
% take the sign with b'y < 0, scale to b'y = -1, and return the largest
% violation of At*y in K*: VIOL = max(|free part|, -min eig of each block),
% absolute since ||b|| = 1 (HEATND_SDP normalizes b); RAD = 1/max(||w_f||_2,
% max_j neg_j), the radius in ||x_f|| + sum tr(X_j) inside which no
% feasible x exists (header); REL = viol / max(1, largest |free entry| or
% |eig|), diagnostic only.
y = res.sol.itr.y;
if ~isempty(keep),  yf = zeros(D.m,1);  yf(keep) = y;  y = yf;  end   % dropped rows: 0
by = D.b'*y;
if by==0,   viol = Inf;     rel = Inf;  rad = 0;    return,     end
y = -y/by;                                      % now b'y = -1
w = D.At*y;     Kf = D.K.f;
viol = max([0; abs(w(1:Kf))]);  wmax = viol;    neg = 0;
off = Kf;
for j = 1:numel(D.K.s)
    s = D.K.s(j);   W = reshape(w(off+(1:s^2)),s,s);    off = off+s^2;
    ev = eig((W+W')/2);
    viol = max(viol,-ev(1));    wmax = max(wmax,max(abs(ev)));  neg = max(neg,-ev(1));
end
viol = full(viol);  rel = viol/max(1,full(wmax));
rad = 1/full(max(norm(w(1:Kf)),neg));           % Inf for an exact ray
end


function v = via_lpisolve(prog,opts,v)
% Cross-check route: sossolve with Mosek. evalc at this fixed call site,
% never inside a handle (R2025b crash, memory matlab-evalc-handle-crash).
evalc('sol = lpisolve(prog,struct(''solver'',''mosek''));');
I = sol.solinfo.info;
[v.rel_b,v.psd_min,v.psd_relmin,v.trivial] = cx_resid(sol);
v.t_mosek = I.cpusec;   v.iter = I.iter;
v.deficient = I.iter==0 && I.feasratio==-1;
v.prosta = sprintf('pinf=%d dinf=%d numerr=%d',I.pinf,I.dinf,I.numerr);
if I.numerr==0 && I.pinf==1
    v.st = -1;  v.why = 'sossolve pinf=1';
elseif I.numerr==0 && I.pinf==0 && v.rel_b<=opts.btol && v.psd_relmin>=-opts.psdtol && ~v.trivial
    v.st = +1;  v.why = 'sossolve feasible, residuals certified';
else
    v.why = 'uncertain';
end
end


function mb = peakMB()
mb = NaN;
try
    p = System.Diagnostics.Process.GetCurrentProcess();     p.Refresh();
    mb = double(p.PeakPagedMemorySize64)/2^20;
catch
end
end
