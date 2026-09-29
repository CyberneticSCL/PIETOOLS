function L = heatNd_ladder(rungs,lopts)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% L = HEATND_LADDER(RUNGS,LOPTS) build-only sweep of the heatNd benchmark
% (HEATND_PIE -> HEATND_LPI -> HEATND_SDP) over rungs (N, d, preset),
% recording per rung: q = numel(decvartable), nP, m, K.f, K.s, nnz(At), nx,
% build time by stage, and memory; optionally saving the pre-solve SDP at
% kappa = LOPTS.dump_k * lambda1 (absolute kappa = r + k; the SDP depends on
% kappa only, so rungs carry no r). These dumps are UNVERIFIED (no
% reference verdict, not chosen to straddle a measured kappa-hat): the
% verified straddling pairs come from HEATND_VALIDATE (VOPTS.dumpdir).
%
% PRICED BEFORE BUILDING. Rungs are run in increasing (N, d, preset), and a
% rung with N >= 2 is priced from the measured rungs below it with the same
% (preset, d), by the Kronecker model the bases obey (Z_d and the R, Q
% bases are Kronecker products over directions, paper eq. 14): each count
% X in {m, nx, nnz} grows by a per-direction factor, so
%   X(N) ~ X(N-1)^2 / X(N-2)          (N >= 3; measured 3-D m 31770 against
%                                       30.2k predicted from 1-D and 2-D),
%   X(2) ~ X(1)^2                      (conservative: a DN factor < DD);
% build time t(N) ~ 10 us x nnz(N) (2x the measured 3-D build + two
% finalizations, 81 s for 17.8M nonzeros); memory mem(N) ~ mem0 + 300 B x
% nnz(N) (2x the 135 B per nonzero of the measured 3-D build: 3790 MB in
% use, 3970 MB peak private, 17.8M nonzeros), plus 32 B per nonzero per
% dump.
% A rung projected above LOPTS.memlim (MB) or LOPTS.tlim (s) is SKIPPED
% and logged with its projection. A rung whose smaller rungs are missing is
% skipped too (it cannot be priced). Memory figures are process-lifetime
% peaks; in increasing size order each is the rung's own peak.
% Also recorded (not run): a MOSEK solve price per route (SPRICE below),
% calibrated on MEASURED 3-D d = 0 process peaks; its m^2 scaling to other
% m is INFERRED.
%
% INPUT
% - rungs: cell of {N, d, preset} (preset optional: 'bench'); default
%          N = 1:3 x d = 0:2 'bench', and N = 1:3, d = 0 'heavy'. A 4th
%          element r is accepted only if 0 (r only shifts k);
% - lopts: memlim (24000), tlim (1200), ep (0.1), dump_k ([]: none;
%          fractions of lambda1), dumpdir (fullfile(tempdir,'heatNd',
%          'dumps') - NOT durable), tag (''), savefile ('': L is saved there
%          after every rung).
% OUTPUT
% - L: struct array, one per rung: N, d, preset, r (0), bc, kstar (=
%      lambda1), status ('built' | 'skipped: ...'), price (projection), q,
%      nP, m, Kf, Ks, nnz, nx, t (stages + total), mem (MB), dumps (files),
%      solve_price (per route, SPRICE).
%
% Initial coding MMP, 09/27/2026
% MMP, 09/27/2026 (review fixes): rungs without r (r = 12 rungs were
%   duplicates at kappa = r + k); dumps named by absolute kappa, off by
%   default; tempdir default (it was this session's scratch path); solve
%   price per route (all rows / independent rows) on measured peaks.
% MMP, 09/27/2026 (final reviews): the independent-row price recalibrated
%   to the measured maxima, 697 s / 11450 MB (it used 586 s / 11100 MB,
%   15-20% below the later keep-route solves).
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if nargin<1 || isempty(rungs)
    rungs = {};
    for d = 0:2,    for N = 1:3,    rungs{end+1} = {N,d,'bench',0};    end,    end %#ok<AGROW>
    for N = 1:3,    rungs{end+1} = {N,0,'heavy',0};    end                          %#ok<AGROW>
end
if nargin<2,    lopts = struct();   end
df = struct('memlim',24000,'tlim',1200,'ep',0.1,'dump_k',[],'dumpdir',fullfile(tempdir,'heatNd','dumps'), ...
            'tag','','savefile','');
fn = fieldnames(df);
for i = 1:numel(fn),    if ~isfield(lopts,fn{i}),   lopts.(fn{i}) = df.(fn{i});     end,    end
if ~isempty(lopts.dump_k) && ~exist(lopts.dumpdir,'dir'),   mkdir(lopts.dumpdir);   end
info = heatNd_path();
fprintf('heatNd_ladder: maxNumCompThreads = %d, memlim %d MB, tlim %d s\n',info.threads,lopts.memlim,lopts.tlim);
% Normalize and order: preset, d, N.
R = cellfun(@(c) [c, repmat({[]},1,4-numel(c))],rungs,'UniformOutput',false);
for i = 1:numel(R)
    if isempty(R{i}{3}),    R{i}{3} = 'bench';  end
    if isempty(R{i}{4}),    R{i}{4} = 0;        end
    if R{i}{4}~=0           % the SDP depends on kappa = r + k only
        error('heatNd_ladder:r','Rungs carry no r: it only shifts k (kappa = r + k).')
    end
end
% N first, so every smaller rung is measured before a larger one is priced
% and the lifetime peak grows with the rung.
key = cellfun(@(c) sprintf('%d|%d|%s|%g',c{1},c{2},c{3},c{4}),R,'UniformOutput',false);
[~,ord] = sort(key);    R = R(ord);
m0 = memnow();  mem0 = m0.matlab;
L = struct('N',{},'d',{},'preset',{},'r',{},'bc',{},'kstar',{},'status',{},'price',{}, ...
           'q',{},'nP',{},'m',{},'Kf',{},'Ks',{},'nnz',{},'nx',{},'t',{},'mem',{}, ...
           'dumps',{},'solve_price',{});
for i = 1:numel(R)
    [N,d,preset,r] = R{i}{:};
    pie = heatNd_pie(N,r);
    e = struct('N',N,'d',d,'preset',preset,'r',r,'bc',{pie.bc},'kstar',pie.exact.kstar, ...
               'status','','price',[],'q',NaN,'nP',NaN,'m',NaN,'Kf',NaN,'Ks',[],'nnz',NaN, ...
               'nx',NaN,'t',[],'mem',[],'dumps',{{}},'solve_price',[]);
    % % % Price.
    if N>=2
        f1 = find_rung(L,N-1,d,preset,r);   f2 = find_rung(L,N-2,d,preset,r);
        if isempty(f1) || (N>=3 && isempty(f2))
            e.status = 'skipped: smaller rungs missing, cannot price';
            L(end+1) = e;   report(e);                                  %#ok<AGROW>
            if ~isempty(lopts.savefile),    save(lopts.savefile,'L');   end
            continue
        end
        a = L(f1);
        if N>=3
            b = L(f2);  pm = a.m^2/b.m;     pnx = a.nx^2/b.nx;  pnz = a.nnz^2/b.nnz;
        else
            pm = a.m^2;     pnx = a.nx^2;   pnz = a.nnz^2;
        end
        pt = 1e-5*pnz;              % 2x the measured 3-D 4.6 us per nonzero (81 s / 17.8M)
        % Memory per nonzero of At, calibrated on the one 3-D build measured
        % (3790 MB in use, ~2.3 GB above idle, for 17.8M nonzeros: 135 B),
        % doubled; a small rung's own increment is noise, so not scaled.
        pmem = mem0 + (300 + 32*numel(lopts.dump_k))*pnz/2^20;
        e.price = struct('m',pm,'nx',pnx,'nnz',pnz,'t',pt,'mem',pmem);
        if pmem>lopts.memlim || pt>lopts.tlim
            e.status = sprintf('skipped: projected %.0f MB / %.0f s (limits %d MB / %d s)', ...
                               pmem,pt,lopts.memlim,lopts.tlim);
            e.solve_price = sprice(pm);
            L(end+1) = e;   report(e);                                  %#ok<AGROW>
            if ~isempty(lopts.savefile),    save(lopts.savefile,'L');   end
            continue
        end
    end
    % % % Build: base, then (15b) at each dump rate, dump.
    o = struct('preset',preset);
    t0 = tic;   [~,meta] = heatNd_lpi(pie,d,[],lopts.ep,o);     tb = toc(t0);
    o.base = meta.base;
    kd = lopts.dump_k;  if isempty(kd),     kd = 1;     end         % shape needs one kappa
    lam1 = pie.exact.lambda1;                   % r = 0: k = kappa
    tf = zeros(1,numel(kd));    td = tf;
    for j = 1:numel(kd)
        t0 = tic;   [prog,mk] = heatNd_lpi(pie,d,kd(j)*lam1,lopts.ep,o);   tf(j) = toc(t0);
        if ~isempty(lopts.dump_k)
            f = fullfile(lopts.dumpdir,sprintf('heatNd%s_N%d_d%d_%s_kap%.6f.mat', ...
                         lopts.tag,N,d,preset,mk.kappa));
            t0 = tic;   D = heatNd_sdp(prog,f,mk);     td(j) = toc(t0);
            e.dumps{end+1} = f;
            clear D
        end
    end
    s = mk.sdp;
    e.q = mk.ndec;  e.nP = mk.nP;   e.m = s.m;  e.Kf = s.Kf;    e.Ks = s.Ks;
    e.nnz = s.nnzAt;    e.nx = s.nx;
    e.t = meta.t;   e.t.base = tb;  e.t.finalize = tf;  e.t.dump = td;  e.t.total = tb+sum(tf);
    e.mem = memnow();   e.status = 'built';     e.solve_price = sprice(e.m);
    clear prog meta mk o
    L(end+1) = e;   report(e);                                          %#ok<AGROW>
    if ~isempty(lopts.savefile),    save(lopts.savefile,'L');   end
end
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function f = find_rung(L,N,d,preset,r)
f = [];
if N<1,     f = 0;  return,     end         % no rung needed (N = 2 has no N-2)
for i = 1:numel(L)
    if L(i).N==N && L(i).d==d && strcmp(L(i).preset,preset) && L(i).r==r && strcmp(L(i).status,'built')
        f = i;
    end
end
end

function p = sprice(m)
% MOSEK solve price per route, calibrated on MEASURED 3-D bench d = 0
% solves (24 MOSEK threads, upper ends of the measured ranges):
%   all rows, m 31770:        MOSEK 499-921 s, process peak 20.1-20.9 GB;
%   independent rows, 22226:  MOSEK 371-697 s, process peak 10.7-11.45 GB
%                             (a process holding the solve only; maxima:
%                             697 s tight at 14.60, 666.6 s defaults in the
%                             dump re-solve, both 11.45 GB);
%   plus, once per family, the HEATND_LINDEP QR: 363 s (d = 0), 476 s (d = 1).
% INFERRED: the m^2 scaling (IPM Schur complement) to other m, and the
% independent fraction 0.70 of the rows (measured 0.70 at d = 0, 0.67 at
% d = 1). Not calibrated: d = 1's keep-route solve peak (its 22.3 GB
% lifetime peak includes the QR in the same process).
sa = (m/31770)^2;   mk = 0.70*m;    sk = (mk/22226)^2;
p = struct('t',921*sa,'mem',21000*sa, ...                       % all rows (legacy fields)
           'all',struct('t',921*sa,'mem',21000*sa), ...
           'keep',struct('m',mk,'t',697*sk,'mem',11450*sk,'t_qr',363*sa));
end

function report(e)
if strcmp(e.status,'built')
    fprintf(['N=%d d=%d %-6s r=%-4g %-12s | q %9d nP %5d m %7d Kf %5d nnz %10d nx %9d | Ks %s | ' ...
             'build %.1f s (R %.1f Q %.1f 15a %.1f) | MATLAB %.0f MB, peak private %.0f MB | ' ...
             'solve price (inferred) all rows %.0f s %.0f MB, keep rows %.0f s %.0f MB\n'], ...
        e.N,e.d,e.preset,e.r,strjoin(e.bc,'x'),e.q,e.nP,e.m,e.Kf,e.nnz,e.nx,ks_str(e.Ks), ...
        e.t.total,e.t.R,e.t.Q,e.t.eq15a,e.mem.matlab,e.mem.peak_private,e.solve_price.all.t, ...
        e.solve_price.all.mem,e.solve_price.keep.t,e.solve_price.keep.mem);
else
    fprintf('N=%d d=%d %-6s r=%-4g %-12s | %s\n',e.N,e.d,e.preset,e.r,strjoin(e.bc,'x'),e.status);
end
end

function s = ks_str(Ks)
% Run-length Gram list, e.g. 288x7 552x7.
[u,~,j] = unique(Ks,'stable');  c = accumarray(j(:),1);
s = strjoin(arrayfun(@(a,b) sprintf('%dx%d',a,b),u(:),c(:),'UniformOutput',false),' ');
end

function m = memnow()
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
