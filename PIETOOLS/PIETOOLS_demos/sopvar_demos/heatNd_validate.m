function T = heatNd_validate(stage,vopts)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% T = HEATND_VALIDATE(STAGE,VOPTS) the validation runs of the heatNd
% benchmark: certified kappa-hat by HEATND_BISECT against the exact
% threshold kappa* = lambda1, and against the stock / paper numbers below.
%
% KAPPA. The Cor. 35 SDP depends on kappa = r + k only (HEATND_PIE), so the
% runs are at r = VOPTS.r (default 0, where k = kappa); another r only
% shifts k. The r-invariant reach is the gap k* - k-hat = lambda1 -
% kappa-hat, reported as the bracket [lambda1 - hi, lambda1 - lo].
%
% STAGE
%   '1d'  N = 1 (DD), d in VOPTS.ds (0:1), presets bench, heavy, listing
%         (Sec. 7.1 as printed), 'product' (listing with the product Psatz)
%         and 'productT' (the same with the stock 'lpivar' P basis: the
%         container analogue of stock lpi_ineq psatz = 1);
%   '2d'  N = 2 (DD x DN, paper Ex. 2 / Table 1), d in VOPTS.ds, presets
%         bench, listingL (the stock listing with linear generators, the
%         only stock configuration found to reproduce a Table 1 row), heavy,
%         listing, listingLT (listingL with the tensor 'lpivar_cdopvar' P);
%         VOPTS.presets (cellstr) restricts either list;
%   '3d'  N = 3 (DD x DN x DN), bench, d = VOPTS.d (0): ONE bisection,
%         resumable across runs (each <= 20 min): the k-affine SDP is saved
%         to VOPTS.kaff_file (default in VOPTS.dir) on the first run, and the
%         trace of final verdicts [kappa st] to VOPTS.dir/trace_<tag>.mat
%         after each run; a later call with the same VOPTS continues.
% After each bisection, straddling dumps are written to VOPTS.dumpdir (''
% for none): the SDP at lo (certified +1) and at hi (certified -1), each
% with its independent row set (info.keep) and the route that produced its
% verdict (info.ref), so HEATND_SOLVE(FILE) repeats it. Only verdicts SOLVED
% in this call are dumped (a reused trace row carries no route).
%
% VOPTS (defaults): dir (fullfile(tempdir,'heatNd','val') - NOT durable),
%   dumpdir (fullfile(dir,'dumps')), r (0), ds (0:1), presets ({} = all),
%   rtol (1e-5 for '1d'/'2d'), tight (true: every verdict solve uses MOSEK
%   tolerances 1e-10; MEASURED 2-D rel_b 1e-10 to 1e-8 instead of up to
%   1.5e-6), retry ({'rows','loose'}: an uncertain verdict is re-solved on
%   the other row set, then also with MOSEK defaults, HEATND_BISECT), tmax
%   (Inf); '3d': d (0), rtol3 (1e-3), tmax3
%   (1080 s), kappa0 ([0.8 1.0] lambda1), lindep3 (true: solve on the
%   independent rows; verdicts on all rows), retry3 ({}: a 3-D retry costs
%   8-15 min), tight3 (true), klist (absolute kappas to solve in order
%   instead of searching), kaff_file ('': dir/kaff_<tag>.mat).
%
% Reference numbers (scout ref, stock PIETOOLS on this tree, MOSEK 11,
% certified brackets [lo, hi] in k at the stated r; kappa = r + k):
%   1-D DD, r = 0 (k* 9.869604): listing d=0 none, d=1 [0, 0.0193];
%       listing+psatz d=0 [7.99995, 8.00010], d=1 [9.869526, 9.879604].
%   1-D DD, r = 9 (k* 0.869604): listing infeasible at k = 0 (d = 0, 1);
%       listing+psatz d=0 infeasible at 0, d=1 [0.869526, 0.879604].
%   2-D, r = 0 (k* 12.337006): paper d=0 12.336, d=1 12.336; stock listing
%       d=0 infeasible at 0; listing + linear psatz d=0 [10.124, 11.062],
%       d=1 [12.084, 12.297].
%   2-D, r = 12 (k* 0.337006): paper d=0 0.33578, d=1 0.33666; stock
%       listing d=0 infeasible at 0; listing + linear d=0 infeasible at 0,
%       d=1 [0.13345, 0.29701].
%
% OUTPUT T: struct array, one row per bisection (N, r, d, preset, lambda1,
% kstar, lo, hi (kappa), khat, gap, unc, resolved, stop, nsolve, m, nx,
% Ks, t, trace, dumps), also saved to VOPTS.dir.
%
% Initial coding MMP, 09/27/2026
% MMP, 09/27/2026 (review fixes): kappa = r + k throughout (r = 0 default;
%   the r = 9 / r = 12 runs were duplicates); gap instead of relative reach;
%   directories are arguments with tempdir defaults (they were this
%   session's scratch path); straddling dumps with row set and route;
%   tight tolerances and a retry by default.
% MMP, 09/27/2026 (final reviews): a dump's info.ref records cert_viol
%   (HEATND_BISECT trace column 9, now the gated absolute Farkas violation).
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if nargin<2,    vopts = struct();   end
df = struct('dir',fullfile(tempdir,'heatNd','val'),'dumpdir',[],'r',0,'ds',0:1,'presets',{{}}, ...
            'rtol',1e-5,'tight',true,'retry',{{'rows','loose'}},'tmax',Inf, ...
            'd',0,'rtol3',1e-3,'tmax3',1080,'kappa0',[],'lindep3',true,'retry3',{{}}, ...
            'tight3',true,'klist',[],'kaff_file','');
fn = fieldnames(df);
for i = 1:numel(fn),    if ~isfield(vopts,fn{i}),   vopts.(fn{i}) = df.(fn{i});     end,    end
if isempty(vopts.dumpdir) && ~ischar(vopts.dumpdir),    vopts.dumpdir = fullfile(vopts.dir,'dumps');    end
if ~exist(vopts.dir,'dir'),     mkdir(vopts.dir);   end
if ~isempty(vopts.dumpdir) && ~exist(vopts.dumpdir,'dir'),  mkdir(vopts.dumpdir);   end
info = heatNd_path();
fprintf('heatNd_validate(%s): maxNumCompThreads = %d (MOSEK uses its own count, per trace row), %s\n', ...
        stage,info.threads,char(datetime('now')));
T = struct('N',{},'r',{},'d',{},'preset',{},'lambda1',{},'kstar',{},'lo',{},'hi',{},'khat',{}, ...
           'gap',{},'unc',{},'resolved',{},'stop',{},'nsolve',{},'m',{},'nx',{},'Ks',{},'t',{}, ...
           'trace',{},'dumps',{});
ep = 0.1;   r = vopts.r;
so = struct('tight',vopts.tight);
switch stage
    case {'1d','2d'}
        if strcmp(stage,'1d')
            N = 1;
            pre = {{'bench',struct()},{'heavy',struct()},{'listing',struct()}, ...
                   {'product',struct('preset','listing','psatz','product')}, ...
                   {'productT',struct('preset','listing','psatz','product','Pbasis','tensor')}};
        else
            N = 2;
            pre = {{'bench',struct()},{'listingL',struct()},{'heavy',struct()},{'listing',struct()}, ...
                   {'listingLT',struct('preset','listingL','Pbasis','tensor')}};
        end
        if ~isempty(vopts.presets),     pre = pre(cellfun(@(c) any(strcmp(c{1},vopts.presets)),pre));  end
        pie = heatNd_pie(N,r);
        for d = vopts.ds
            for p = 1:numel(pre)
                o = pre{p}{2};  if ~isfield(o,'preset'),  o.preset = pre{p}{1};  end
                fprintf('\n-- N=%d d=%d %s (r = %g, lambda1 = kappa* = %.6f)\n',N,d,pre{p}{1},r,pie.exact.lambda1);
                B = heatNd_bisect(pie,d,ep,o,struct('rtol',vopts.rtol,'tmax',vopts.tmax, ...
                                  'solve',so,'retry',{vopts.retry}));
                dumps = write_dumps(B,pre{p}{1},vopts.dumpdir);
                T(end+1) = row(N,r,d,pre{p}{1},B,dumps);                             %#ok<AGROW>
                report(T(end));
            end
        end
        save(fullfile(vopts.dir,sprintf('val_%s.mat',stage)),'T');
    case '3d'
        N = 3;  d = vopts.d;
        pie = heatNd_pie(N,r);  lam1 = pie.exact.lambda1;
        tag = sprintf('N3_d%d_bench',d);
        kf = vopts.kaff_file;   if isempty(kf),     kf = fullfile(vopts.dir,['kaff_' tag '.mat']);  end
        tf = fullfile(vopts.dir,['trace_' tag '.mat']);
        tr0 = zeros(0,2);   hist = {};
        if exist(tf,'file'),    S = load(tf);   tr0 = S.tr0;    hist = S.hist;  end
        k0 = vopts.kappa0;  if isempty(k0),     k0 = [0.8 1.0]*lam1;    end
        bo = struct('rtol',vopts.rtol3,'tmax',vopts.tmax3,'trace0',tr0,'kaff_file',kf,'kappa0',k0, ...
                    'lindep',vopts.lindep3,'klist',vopts.klist,'solve',struct('tight',vopts.tight3), ...
                    'retry',{vopts.retry3});
        B = heatNd_bisect(pie,d,ep,struct('preset','bench'),bo);
        n0 = size(tr0,1);   new = B.trace(n0+1:end,:);  whynew = B.why(n0+1:end);
        fin = final_rows(new);                          % one final verdict per solved kappa
        tr0 = [tr0; new(fin,1:2)];  hist(end+1,:) = {new, whynew};         %#ok<AGROW>
        save(tf,'tr0','hist');                          % verdicts so far, for resuming
        dumps = write_dumps(B,'bench',vopts.dumpdir);
        B = rmfield(B,'kaff');                          % the SDP itself is in kf
        T(end+1) = row(N,r,d,'bench',B,dumps);
        report(T(end));
        save(fullfile(vopts.dir,['val_3d_' tag '.mat']),'T','B');
    otherwise
        error('heatNd_validate:stage','STAGE is ''1d'', ''2d'' or ''3d''.')
end
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function s = row(N,r,d,preset,B,dumps)
s = struct('N',N,'r',r,'d',d,'preset',preset,'lambda1',B.lambda1,'kstar',B.kstar, ...
           'lo',B.lo,'hi',B.hi,'khat',B.khat,'gap',B.gap,'unc',B.unc,'resolved',B.resolved, ...
           'stop',B.stop,'nsolve',B.nsolve,'m',B.sdp.m,'nx',B.sdp.nx,'Ks',B.sdp.Ks,'t',B.t_total, ...
           'trace',B.trace,'dumps',{dumps});
end


function report(s)
fprintf(['   => kappa in [%.6f, %.6f)  (k-hat %.6f at r = %g)  gap k*-khat in [%.2e, %.2e]  ' ...
         'resolved %d (%s)  uncertain in bracket: %s  (%d attempts, %.1f s)\n'], ...
        s.lo,s.hi,s.khat,s.r,s.gap(1),s.gap(2),s.resolved,s.stop,mat2str(s.unc,8),s.nsolve,s.t);
end


function i = final_rows(tr)
% Index of the last attempt at each solved kappa (the one the bracket uses).
i = [];
for k = unique(tr(:,1))'
    j = find(tr(:,1)==k & tr(:,12)>0,1,'last');     i(end+1) = j;                   %#ok<AGROW>
end
i = sort(i);
end


function files = write_dumps(B,preset,ddir)
% SDP at lo (+1) and at hi (-1), each with the route of its verdict, from
% the k-affine family the verdict was solved on. Only attempts of this call.
files = {};
if isempty(ddir) || isempty(B.kaff) || ~isfield(B.kaff,'D1'),   return,     end
K = B.kaff;     tr = B.trace;
for side = [1 -1]
    if side==1,     kap = B.lo;     else,   kap = B.hi;     end
    if isnan(kap),  continue,   end
    j = find(tr(:,1)==kap & tr(:,12)>0,1,'last');
    if isempty(j) || tr(j,2)~=side,     continue,   end     % reused or not final
    D = K.D1;   D.At = K.D1.At + (kap-K.kap1)*K.dA;     D.kappa = kap;
    if isfield(K,'keep'),   D.keep = K.keep;
    else,                   D.keep = heatNd_lindep(D.At);   % stored for 'keep' re-solves
    end
    rows = {'all','keep'};
    D.ref = struct('st',side,'rows',rows{tr(j,10)+1},'tight',logical(tr(j,11)),'rel_b',tr(j,3), ...
                   'psd_relmin',tr(j,4),'cert_viol',tr(j,9),'t_mosek',tr(j,5),'iter',tr(j,6), ...
                   'mosek_threads',tr(j,8),'attempt',tr(j,12),'why',B.why{j}, ...
                   'when',char(datetime('now')),'source','heatNd_validate / heatNd_bisect');
    mt = B.meta;    mt.kappa = kap;     mt.k = kap - B.r;   mt.r = B.r;
    f = fullfile(ddir,sprintf('heatNd_N%d_d%d_%s_kap%.6f.mat',mt.N,mt.d,preset,kap));
    heatNd_sdp(D,f,mt);     files{end+1} = f;                                       %#ok<AGROW>
end
end
