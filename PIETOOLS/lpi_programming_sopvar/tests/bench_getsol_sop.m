function T = bench_getsol_sop(opts)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% T = BENCH_GETSOL_SOP(OPTS) measures the time and memory of
% getsol_lpivar_sop over a sweep of the number of decision variables q and
% of spatial variables nv (CLAUDE.md sec. 2), and compares it with the
% legacy paths at small q. Returns a struct of tables, prints them.
%
% Operators are real LPI decision containers, from 'lpivar_cdopvar' (one
% nonzero of B per variable) and 'poscopvar' (a Gram form: several
% nonzeros per variable, names listed twice in decvartable). The values
% are random: getsol reads only decvartable and solinfo.RRx, so no solve
% is needed. Each size is timed on
%   'own'  the program's decvartable as declared (q entries for
%          lpivar_cdopvar, n^2 for a Gram of size n);
%   'pad'  that table plus 2q absent names, shuffled (one operator among
%          several); column N is this table's length,
% and the two results must be bit-identical. Peak memory is the MATLAB
% profiler's PeakMem of getsol_lpivar_sop on the 'pad' path (a separate
% call, since profiling perturbs the timing).
%
% Legacy comparison, 1-D only, same q: getsol_lpivar on the 'dopvar' from
% 'lpivar' (for which lpivar_cdopvar reproduces the family exactly), and
% the route of test_lpi_eq_sdopvar_endtoend's recover_dvals, lpigetsol on
% a 'dpvar' of the q names.
%
% OPTS fields (defaults): q [1e4 1e5 1e6 3e6], nv [1 2 3], deg 1, reps 3,
% gram_q [1e4 1e5 1e6] (targets; q ~ m^2 is rounded), gram_nv 1,
% legacy_q [1e3 1e4 1e5 1e6 3e6], legacy_cap 120 (s), multi true. An empty
% size list skips its section. Measured 2026-09-29 (i9-14900KF, R2025b):
% 567 s for the defaults, 3.7 GB MATLAB memory at the end, with another
% MATLAB session running (getsol times then ~1.6x those of a quiet run).
%
% Initial coding MMP, 09/29/2026
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

dflt = struct('q',[1e4 1e5 1e6 3e6],'nv',[1 2 3],'deg',1,'reps',3, ...
              'gram_q',[1e4 1e5 1e6],'gram_nv',1,'legacy_q',[1e3 1e4 1e5 1e6 3e6],'legacy_cap',120,'multi',true);
if nargin<1,    opts = struct();    end
fn = fieldnames(dflt);
for k = 1:numel(fn),    if ~isfield(opts,fn{k}),    opts.(fn{k}) = dflt.(fn{k});    end,    end
V = {'s1','s2','s3'};
T = struct();

% % % 1. lpivar_cdopvar sweep over (nv, q), one square block.
rows = {};
for nv = opts.nv
    vars = V(1:nv);     dom = repmat([0 1],nv,1);
    q1 = count_q(vars,dom,1,opts.deg);
    for qt = opts.q
        m = max(1,round(sqrt(qt/q1)));
        prog = mkprog(vars,dom);
        t0 = tic;   [prog,P] = lpivar_cdopvar(prog,m,{vars},dom,opts.deg);  tb = toc(t0);
        rows(end+1,:) = measure(sprintf('lpivar nv=%d',nv),nv,prog,P,tb,opts.reps); %#ok<AGROW>
        print_row(rows(end,:));
        clear prog P
    end
end
if ~isempty(rows),  T.lpivar = cell2table(rows,'VariableNames',hdr());  end

% % % 2. Multi-block: R^m x L2^m[vars], 2x2 blocks on one list.
if opts.multi
    rows = {};
    for nv = opts.nv
        vars = V(1:nv);     dom = repmat([0 1],nv,1);
        q1 = count_q_multi(vars,dom,1,opts.deg);
        qt = opts.q(end);
        m = max(1,round(sqrt(qt/q1)));
        prog = mkprog(vars,dom);
        sp = {{},vars};
        t0 = tic;
        [prog,P] = lpivar_cdopvar(prog,struct('out',[m;m],'in',[m;m]), ...
                                  struct('out',{sp},'in',{sp}),dom,opts.deg);
        tb = toc(t0);
        rows(end+1,:) = measure(sprintf('2x2 multi nv=%d',nv),nv,prog,P,tb,opts.reps); %#ok<AGROW>
        print_row(rows(end,:));
        clear prog P
    end
    T.multi = cell2table(rows,'VariableNames',hdr());
end

% % % 3. Gram form (poscopvar): several nonzeros per variable, names
% % % listed twice in decvartable. q ~ m^2 for a Gram of m*(basis) rows.
rows = {};
for nv = opts.gram_nv
    vars = V(1:nv);     dom = repmat([0 1],nv,1);
    if isempty(opts.gram_q),    break,  end
    prog = mkprog(vars,dom);
    [~,P1] = poscopvar(prog,1,{vars},dom,opts.deg);    q1 = numel(P1.Zd);
    for qt = opts.gram_q
        m = max(1,round(sqrt(qt/q1)));
        prog = mkprog(vars,dom);
        t0 = tic;   [prog,P] = poscopvar(prog,m,{vars},dom,opts.deg);   tb = toc(t0);
        rows(end+1,:) = measure(sprintf('poscopvar nv=%d',nv),nv,prog,P,tb,opts.reps); %#ok<AGROW>
        print_row(rows(end,:));
        clear prog P
    end
end
if ~isempty(rows),  T.gram = cell2table(rows,'VariableNames',hdr());    end

% % % 4. Legacy paths, 1-D. A route stops where its projected cost at the
% % % next size exceeds legacy_cap: lpivar build and getsol_lpivar ~q,
% % % the dpvar-of-names route ~q^2.25 (measured 1e3..3e4, 09/29).
rows = {};  skipA = false;  skipB = false;
for iq = 1:numel(opts.legacy_q)
    qt = opts.legacy_q(iq);
    rq = inf;   if iq<numel(opts.legacy_q),  rq = opts.legacy_q(iq+1)/qt;   end
    q1 = count_q({'s1'},[0 1],1,opts.deg);
    m = max(1,round(sqrt(qt/q1)));
    pvar s1 s1_dum
    progL = lpiprogram(s1,s1_dum,[0 1]);
    d3 = opts.deg*[1 1 1];
    tA = NaN;   tB = NaN;   tbL = NaN;  qd = NaN;
    if ~skipA
        t0 = tic;   [progL,Pd] = lpivar(progL,[0 0;m m],d3);    tbL = toc(t0);  % legacy dopvar
        qd = numel(progL.decvartable);
    end
    [progL,Pc] = lpivar_cdopvar(progL,m,{{'s1'}},[0 1],d3);   % same family, own variables
    NL = numel(progL.decvartable);
    progL.solinfo.RRx = randn(NL,1);    progL.solinfo.info = struct('pinf',0);
    qc = numel(Pc.Zd);
    t0 = tic;   Sc = getsol_lpivar_sop(progL,Pc);   tn = toc(t0); %#ok<NASGU>
    if ~skipA
        t0 = tic;   Sd = getsol_lpivar(progL,Pd);   tA = toc(t0); %#ok<NASGU>
        skipA = rq*max(tA,tbL) > opts.legacy_cap;
    end
    if ~skipB
        t0 = tic;   dv = double(lpigetsol(progL,dpvar(cellstr(Pc.Zd(:)))));   tB = toc(t0);
        skipB = tB*rq^2.25 > opts.legacy_cap;
        assert(isequal(dv(:),progL.solinfo.RRx(end-qc+1:end)),'recover_dvals route read other values')
    end
    rows(end+1,:) = {qd,qc,tn,tA,tB,tbL}; %#ok<AGROW>
    fprintf('legacy q=%8d (cdopvar q=%8d): getsol_lpivar_sop %8.4f s | getsol_lpivar(dopvar) %8.3f s | lpigetsol(dpvar of q names) %8.3f s | lpivar build %7.2f s\n',...
        qd,qc,tn,tA,tB,tbL);
    clear progL Pd Pc Sd Sc
end
if ~isempty(rows),  T.legacy = cell2table(rows,'VariableNames',{'q_dopvar','q_cdopvar','t_sop','t_getsol_lpivar','t_dpvar_names','t_lpivar_build'});  end
end


% ------------------------------------------------------------------------
function h = hdr()
h = {'family','nv','q','N','nnzB','ncell','build_s','t_own','t_pad','peak_MB','alloc_MB','in_MB'};
end

function print_row(r)
fprintf('%-16s nv=%d q=%9d N=%9d nnzB=%10d cells=%4d | build %7.1f s | getsol own %8.4f s, pad %8.4f s | peak %8.1f MB, alloc %8.1f MB (B+Zd %8.1f MB)\n',r{:});
end

function row = measure(lbl,nv,prog,P,tb,reps)
% 'own': the program as declared, N = its decvartable. 'pad': the same
% table plus 2q names the operator lacks, shuffled, N = numel + 2q, as for
% one operator among several. Values are drawn per NAME, so that a name
% the table lists twice (Gram (i,j),(j,i)) has one value and both tables
% must give bit-identical results.
dt = prog.decvartable(:);   q = numel(P.Zd);
[u,~,j] = unique(dt);   val = randn(numel(u),1);
prog.solinfo.RRx = val(j);      prog.solinfo.info = struct('pinf',0);
pad = strcat('zz_pad_',cellstr(string((1:2*q)')));
dt2 = [dt; pad];    rr2 = [val(j); randn(2*q,1)];   pm = randperm(numel(dt2));
prog2 = prog;   prog2.decvartable = dt2(pm);    prog2.solinfo.RRx = rr2(pm);
N = numel(dt2);
[nnzB,ncell,inMB] = op_size(P);
tc = inf;   ts = inf;
for r = 1:reps
    t0 = tic;   S1 = getsol_lpivar_sop(prog,P);    tc = min(tc,toc(t0));
    t0 = tic;   S2 = getsol_lpivar_sop(prog2,P);   ts = min(ts,toc(t0));
end
assert(isequal(S1,S2),'%s: own and padded/shuffled tables disagree',lbl)
clear S1 S2
profile('-memory','on');
S2 = getsol_lpivar_sop(prog2,P); %#ok<NASGU>
profile('off');
info = profile('info');     ft = info.FunctionTable;
k = find(strcmp({ft.FunctionName},'getsol_lpivar_sop'),1);
pk = ft(k).PeakMem/2^20;    al = ft(k).TotalMemAllocated/2^20;
profile('clear');
row = {lbl,nv,q,N,nnzB,ncell,tb,tc,ts,pk,al,inMB};
end

function [nnzB,ncell,MB] = op_size(P)
nnzB = 0;   ncell = 0;  by = 0;
for k = 1:numel(P.C)
    b = P.C{k};
    if ~isa(b,'sdopvar'),   continue,   end
    prm = b.params;
    for g = 1:numel(prm.B)
        nnzB = nnzB + nnz(prm.B{g});    ncell = ncell + 1;
        Bg = prm.B{g};  w = whos('Bg');     by = by + w.bytes;
    end
end
Zd = P.Zd;  w = whos('Zd');     MB = (by + w.bytes)/2^20;
end

function q = count_q(vars,dom,m,deg)
prog = mkprog(vars,dom);
[~,P] = lpivar_cdopvar(prog,m,{vars},dom,deg);
q = numel(P.Zd);
end

function q = count_q_multi(vars,dom,m,deg)
prog = mkprog(vars,dom);    sp = {{},vars};
[~,P] = lpivar_cdopvar(prog,struct('out',[m;m],'in',[m;m]),struct('out',{sp},'in',{sp}),dom,deg);
q = numel(P.Zd);
end

function prog = mkprog(vars,dom)
% lpiprogram refuses more than 2 variables (lpiprogram.m); for 3, build
% what it would return, as heatNd_lpi does.
nv = numel(vars);
if nv<=2
    prog = lpiprogram(polynomial(vars(:)),[],dom);
else
    prog = sosprogram(polynomial([]),dpvar(zeros(0,1)));
    prog.vartable = [prog.vartable; polynomial(vars(:)); polynomial(strcat(vars(:),'_dum'))];
    prog.dom = dom;
end
end
