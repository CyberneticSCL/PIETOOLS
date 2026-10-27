function R = cx_compare(cases,opts)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% R = CX_COMPARE(CASES,OPTS) solves each executive LPI twice - through the
% stock opvar/opvar2d executive and through its container transcription
% cx_<exec> - and compares them.
%
% CASES: struct array (or cell of structs) with fields
%   id      label
%   exec    executive name without 'PIETOOLS_', e.g. 'Hinf_gain'; the
%           container transcription is cx_<exec>
%   setname lpisettings name, e.g. 'light'
%   solver  'mosek' (1-D default) or 'sedumi' (2-D)
%   kind    'obj'  a min-gamma executive: container posed at fixed gamma
%           'feas' a feasibility executive (stability, well-posedness)
%   plant   cell of cx_plant arguments, e.g. {'io1'} or {'rd',0.5}
%   thresh  (feas, optional) true to also bisect the plant parameter
%           plant{2} (frac) for the stability threshold on both paths
%
% Per 'obj' case:
%   stock   the unmodified executive, real solve: gam*, shape, solver time;
%   pinned  the SAME stock program captured before the solve and posed at
%           fixed gamma (cx_stock_pinned), bisected on certified verdicts:
%           the stock feasibility threshold, independent of how well the
%           objective solve converged (several stock cases end numerr = 2);
%   cx      the container program at fixed gamma, bisected the same way.
% Per 'feas' case: stock and container verdicts at the plant, and with
% thresh, both thresholds in frac.
% Always: the container program's shape (decision variables, free block, PSD
% block sizes, equality rows) against the stock one, and assembly / solve
% times. The container program has no gamma variable, so its free block is
% one smaller than the stock's when the stock gamma is free.
%
% OPTS: rtol (bisection, default 1e-3), verbose (default true).
%
% Initial coding MMP, 09/25/2026
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if nargin<2,    opts = struct();    end
if ~isfield(opts,'rtol'),       opts.rtol = 1e-3;       end
if ~isfield(opts,'verbose'),    opts.verbose = true;    end
if iscell(cases),   cases = [cases{:}];     end
warning('off','sopvar:noncanonicalMultiplier');
warning('off','sdopvar:noncanonicalMultiplier');

R = cell(1,numel(cases));
for ic = 1:numel(cases)
    c = cases(ic);
    if ~isfield(c,'solver') || isempty(c.solver),   c.solver = 'mosek';   end
    st = cx_settings(c.setname,c.solver);
    % Optional per case: SOSTOOLS' psimplify (off by default, as in the
    % baseline). The 2-D equality systems are ~22% rank-deficient on both
    % paths and SeDuMi does not certify them without it.
    if isfield(c,'simplify') && ~isempty(c.simplify),  st.sos_opts.simplify = c.simplify;  end
    % Optional per case: a function handle applied to the settings (e.g. a
    % different eppos), identically on both paths.
    if isfield(c,'stmod') && ~isempty(c.stmod),  st = c.stmod(st);  end
    sopts = st.sos_opts;
    PIE = cx_plant(c.plant{:});
    cxfn = ['cx_' c.exec];
    r = struct('id',c.id,'exec',c.exec,'kind',c.kind,'solver',c.solver,'err','');
    try
        switch c.kind
            case 'obj'
                % stock, real objective solve
                [sol,gam,tw] = stock_run(c.exec,PIE,st);
                r.stock = struct('gam',gam,'numerr',sol.solinfo.info.numerr, ...
                    'cpusec',sol.solinfo.info.cpusec,'wall',tw,'shape',cx_shape(sol));
                % stock program at fixed gamma, bisected
                prog0 = cx_stock_capture(c.exec,PIE,st);
                g0 = gam;   if ~(isfinite(g0) && g0>0),  g0 = 1;   end
                r.pinned = cx_bisect(@(g) cx_stock_pinned(prog0,g,sopts),g0,opts.rtol);
                % container at fixed gamma, bisected; assembly timed apart
                f = @(g) cx_at(cxfn,PIE,st,g,sopts);
                r.cx = cx_bisect(f,g0,opts.rtol);
                gs = r.pinned.hi;   if ~isfinite(gs),   gs = g0;    end
                [rc,ta] = cx_at(cxfn,PIE,st,gs,sopts);
                r.cxshape = rc.shape;   r.cxassembly = ta;  r.cxcpusec = rc.cpusec;
                r.agree = agree_brackets(r.pinned,r.cx);
            case 'feas'
                [sol,~,tw] = stock_run(c.exec,PIE,st);
                rs = verdict_of(sol);
                r.stock = struct('st',rs,'numerr',sol.solinfo.info.numerr, ...
                    'cpusec',sol.solinfo.info.cpusec,'wall',tw,'shape',cx_shape(sol));
                [r.stock.rel_b,~,r.stock.psd_relmin,r.stock.trivial] = cx_resid(sol);
                t = tic;    prog = feval(cxfn,PIE,st);  ta = toc(t);
                rc = cx_solve(prog,sopts);
                r.cxst = rc.st;     r.cxshape = rc.shape;   r.cxassembly = ta;
                r.cxcpusec = rc.cpusec;     r.cxnumerr = rc.numerr;
                r.cxrel_b = rc.rel_b;   r.cxpsd_relmin = rc.psd_relmin;     r.cxtrivial = rc.trivial;
                % Contradiction = opposite CERTIFIED verdicts; an uncertified
                % side contradicts nothing (as for the gamma brackets).
                r.agree = rs*rc.st ~= -1;
                if isfield(c,'thresh') && ~isempty(c.thresh) && c.thresh
                    % frac destabilizes (feasible BELOW the threshold) and
                    % cx_bisect wants feasible above, so bisect g = 1/frac
                    % and report the bracket back in frac.
                    pl = c.plant;
                    fs = @(g) stock_feas(c.exec,pl,1/g,st);
                    fc = @(g) cx_solve(feval(cxfn,cx_plant(pl{1},1/g,pl{3:end}),st),sopts);
                    r.sthresh = inv_bracket(cx_bisect(fs,1/pl{2},1e-2));
                    r.cthresh = inv_bracket(cx_bisect(fc,1/pl{2},1e-2));
                    r.agree = r.agree && agree_brackets(r.sthresh,r.cthresh);
                end
        end
    catch ME
        r.err = sprintf('[%s] %s (%s:%d)',ME.identifier,ME.message, ...
                        ME.stack(1).name,ME.stack(1).line);
    end
    if opts.verbose,    print_row(r);   end
    R{ic} = r;
end
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [r,ta] = cx_at(cxfn,PIE,st,g,sopts)
t = tic;    prog = feval(cxfn,PIE,st,g);    ta = toc(t);
r = cx_solve(prog,sopts);
end

function [sol,gam,tw] = stock_run(exec,PIE,st)
% The unmodified executive; gam is its first numeric scalar output (the
% executives do not agree on the output position, see bl_b_syn1).
fn = ['PIETOOLS_' exec];
nout = abs(nargout(fn));    out = cell(1,nout);
t = tic;    evalc('[out{1:nout}] = feval(fn,PIE,st);');  tw = toc(t);
sol = out{1};   gam = NaN;
for k = 2:nout
    if isnumeric(out{k}) && isscalar(out{k}),   gam = double(out{k});   break,  end
end
end

function r = stock_feas(exec,pl,fr,st)
[sol,~,~] = stock_run(exec,cx_plant(pl{1},fr,pl{3:end}),st);
r = struct('st',verdict_of(sol));
end

function s = verdict_of(sol)
q = sol.solinfo.info;   s = 0;
if q.numerr==0 && q.pinf==0,    s = +1;     end
if q.numerr==0 && q.pinf==1,    s = -1;     end
end

function B = inv_bracket(B)
% A bracket in g = 1/frac restated in frac: [1/hi, 1/lo], i.e. the
% certified-FEASIBLE frac first and the certified-infeasible one second.
lo = B.lo;  hi = B.hi;
B.lo = 1/hi;    B.hi = 1/lo;    B.est = sqrt(B.lo*B.hi);
B.trace(:,1) = 1./B.trace(:,1);
end

function tf = agree_brackets(A,B)
% Two certified brackets CONTRADICT iff one side's certified-infeasible end
% lies at or above the other's certified-feasible end. A missing end (NaN:
% no certified point on that side, e.g. an uncertified band or a plant whose
% threshold is 0) cannot contradict anything; whether both brackets are
% two-sided is reported separately (field 'resolved' of each bracket).
tf = (isnan(A.lo) || isnan(B.hi) || A.lo < B.hi) && ...
     (isnan(B.lo) || isnan(A.hi) || B.lo < A.hi);
end

function print_row(r)
if ~isempty(r.err)
    fprintf('  %-16s ERROR %s\n',r.id,r.err);     return
end
S = r.stock.shape;  C = r.cxshape;
sh = sprintf('ndv %d/%d  Kf %d/%d  m %d/%d  Ks %s / %s',S.ndv,C.ndv,S.Kf,C.Kf,S.m,C.m, ...
             mat2str(S.Ks),mat2str(C.Ks));
switch r.kind
    case 'obj'
        fprintf('  %-16s gam* %.6g (numerr %d) | pinned [%.6g %.6g] | cx [%.6g %.6g] | agree %d | %s\n', ...
            r.id,r.stock.gam,r.stock.numerr,r.pinned.lo,r.pinned.hi,r.cx.lo,r.cx.hi,r.agree,sh);
    case 'feas'
        s = sprintf(['  %-16s stock %+d (numerr %d, rel_b %.1e, psd %+.1e%s) | cx %+d (numerr %d, '...
            'rel_b %.1e, psd %+.1e%s) | agree %d'],r.id,r.stock.st,r.stock.numerr, ...
            r.stock.rel_b,r.stock.psd_relmin,triv(r.stock.trivial),r.cxst,r.cxnumerr, ...
            r.cxrel_b,r.cxpsd_relmin,triv(r.cxtrivial),r.agree);
        if isfield(r,'sthresh')
            s = [s sprintf(' | frac* stock [%.4g %.4g] cx [%.4g %.4g]',r.sthresh.lo, ...
                 r.sthresh.hi,r.cthresh.lo,r.cthresh.hi)];
        end
        fprintf('%s | %s\n',s,sh);
end
end

function s = triv(t)
if t,   s = ', X=0';    else,   s = '';     end
end
