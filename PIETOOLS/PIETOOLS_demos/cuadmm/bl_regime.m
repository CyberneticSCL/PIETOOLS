function bl_regime(block)                                                   % CC, 09/26/2026
% BL_REGIME  One block of the 2026-09-26 cuADMM testing regime, run as its own
% MATLAB process by gpu/regime.sh (so a crash or out-of-memory loses one block).
%
%   bl_regime('B1a')
%
% Every item is one (case, solver, purpose) run that completes on its own.
% Items already 'ok' in <cuadmm_outdir>/regime/regime.tsv are skipped, so a
% block can be re-run after a crash.  An item is not started if its estimated
% end is past $REGIME_DEADLINE (posix s); every bl_bisect call also gets that
% deadline, the STOP file (<cuadmm_outdir>/STOP) and a per-call timeout, so no
% single solver call can overrun it.  Each item's result struct is saved as
% regime/<key>.mat; later blocks read the Mosek brackets from there.
%
% BLOCKS (question -> readout; see the proposal in the session log / memory
% note cuadmm-10h-regime-proposal):
%   B0   fresh 1-D (+nl, stab2) dumps on the current tree
%   B1a  Mosek (rtol 1e-6) and SeDuMi bisections on the 7 1-D objective cases;
%        est_rd1 also with the pin row rescaled.  Gold standard for B1c/B1d.
%   B2S  lambda_LPI of the heavy 1-D stability family by rebuild; sentinels just
%        past it and at 1.01 lambda*: an F there voids c = 0 F verdicts.
%   B1c  cuADMM bisection on hinf_rd1_hv/hinfdu_rd1/hinfco_rd1, psd_tol held at
%        1e-7, plus a held-out certification run at Mosek's tight gamma_I.
%   B1d  cuADMM bisection on ctrl/hinfduco/est only if B1a certified them.
%   B2   certificates on the 9 1-D feasibility cases, Mosek and cuADMM.
%   B3b  Mosek bisection on nl_fisher_opt (m ~ 5k).
%   B3a  can a cuADMM point be certified at m ~ 5k: probes at 1.2 gamma_F.
%   B3c  cuADMM bisection on nl_fisher_opt, only if B3a certified.
%   B5   2-D stability with the four linear psatz generators: Mosek and cuADMM.
%   B4   scale_hinf n02/n04/n08: Mosek bisection and cuADMM ms/iteration only.
% CC, 09/27/2026: Sol-readiness blocks (night 1 of the 09-27 plan; Sol has no
%   Mosek, so verdicts at size must be checked where the answer is known):
%   N1   the driver itself: resume skips a done item, the deadline and STOP
%        refuse to start one, bl_bisect's STOP check fires before any solver
%        call.  Results in regime/N1_checks.txt; any FAIL -> DRIVER_TRIP.
%   N2b  infeasible seeds: hi at 0.99/0.999 x Mosek's certified gamma_I, no
%        doubling, SeDuMi and cuADMM.  The first certification has no
%        norm_ratio reference, so this is the unguarded case.  F -> SOUNDNESS_TRIP.
%   S1   scale ladder: scale_stab n08/16/24/32 at 0.5 lambda* (feasible by
%        construction) and 1.01 lambda* sentinels at n16/24/32 (infeasible by
%        theory), cuADMM under the default rule.  Readouts: verdicts, iterations
%        to certificate, ms/it, t_init, t_cert, memory vs m, to price Sol jobs.
%        Decoupled copies: iteration counts are expected n-invariant (519 / 4023
%        to 1e-4 / 1e-6 at every n on 09-24), so this measures cost per
%        iteration and soundness at size, not convergence of a coupled problem.
%   N4   the harness at m = 80550 (Th_n3, 2-D, 26.1M nonzeros): re-certify the
%        kept 1e-6 cuADMM point (certification cost, and whether it passes),
%        then one live run capped at the 1e-4 level (dump, setup, certification
%        end to end).

cuadmm_path;
OUT = cuadmm_outdir();
RG  = fullfile(OUT,'regime');  if ~exist(RG,'dir'), mkdir(RG); end
MAN = fullfile(RG,'regime.tsv');
DL  = str2double(getenv('REGIME_DEADLINE'));  if isnan(DL), DL = Inf; end
G.OUT = OUT;  G.RG = RG;  G.MAN = MAN;  G.DL = DL;
G.DMP = fullfile(OUT,'baseline','dumps');
G.BO  = fullfile(RG,'bisect');  if ~exist(G.BO,'dir'), mkdir(G.BO); end
if ~exist(MAN,'file')
    fid = fopen(MAN,'w');
    fprintf(fid,'key\tstatus\tstart\twall_s\tgamma_I\tgamma_F\tgamma_i\tgamma_f\tF\tf\ti\tI\titers\tt_solver\tflags\tnote\n');
    fclose(fid);
end
fprintf('REGIME %s start %s  deadline %s\n',block,datestr(now,31),stamp(DL));

ONE = {'hinf_rd1','hinf_rd1_hv','hinfdu_rd1','hinfco_rd1','hinfduco_rd1','ctrl_rd1','est_rd1'};
FEAS = {'stab_rd1','stabdual_rd1','stabpde_rd1','stabpded_rd1','stab_rd1_hv', ...
        'stab_rd1_tight','stab_tr1','stab_wave1','wellposed_rd1'};
%cu = @(extra) merge(struct('solver','cuadmm','deadline',DL,'outdir',G.BO),extra); % CC, 09/27/2026 (was)
% CC, 09/27/2026: CUADMM_REQUIRE_FULL_GPU=1 (set by sol/sol_harness.slurm) makes
% every cuADMM run refuse a MIG slice
cu = @(extra) merge(struct('solver','cuadmm','deadline',DL,'outdir',G.BO, ...
     'require_full_gpu',strcmp(getenv('CUADMM_REQUIRE_FULL_GPU'),'1')),extra); % CC, 09/27/2026
ms = @(extra) merge(struct('solver','mosek','deadline',DL,'outdir',G.BO),extra);
sd = @(extra) merge(struct('solver','sedumi','deadline',DL,'outdir',G.BO),extra);

switch block
case 'B0'
    ids = [FEAS, ONE, {'nl_fisher_opt','nl_fisher_nobnd','stab2_rd'}];
    item(G,'B0|dumps',10*60,@() dumps(ids));

case 'B1a'
    for k = 1:numel(ONE)
        id = ONE{k};  hi = hi_of(G,id);
        item(G,['B1a|mosek|' id],120,@() bl_bisect(dmp(G,id),ms(struct('hi',hi,'rtol',1e-6, ...
             'max_probes',80,'hi_doublings',12,'tag','B1a'))));
        item(G,['B1a|sedumi|' id],300,@() bl_bisect(dmp(G,id),sd(struct('hi',hi,'rtol',1e-4, ...
             'max_probes',60,'hi_doublings',12,'tag','B1a'))));
    end
    hi = hi_of(G,'est_rd1');
    item(G,'B1a|mosek|est_rd1_pin',120,@() bl_bisect(dmp(G,'est_rd1'),ms(struct('hi',hi, ...
         'rtol',1e-6,'max_probes',80,'hi_doublings',12,'pin_scale','auto','tag','B1a_pin'))));
    item(G,'B1a|sedumi|est_rd1_pin',300,@() bl_bisect(dmp(G,'est_rd1'),sd(struct('hi',hi, ...
         'rtol',1e-4,'max_probes',60,'hi_doublings',12,'pin_scale','auto','tag','B1a_pin'))));

case 'B2S'
    item(G,'B2S|lambda',10*60,@() bl_lambda('rd','heavy','PIE2PDEstability',0.95,1.0,1e-4,'lamhv',DL));
    Lr = getres(G,'B2S|lambda');
    if ~isempty(Lr)
        fr = [Lr.lam_I, 1.01];  nm = {'sentLPI','sent101'};
        for q = 1:2
            id = strrep(sprintf('%s_f%.6f',nm{q},fr(q)),'.','p');
            c  = struct('id',id,'cls','sentinel','dim','1D','kind','feas', ...
                        'builder','bl_b_stab1','args',{{'rd',fr(q),'heavy','PIE2PDEstability'}});
            item(G,['B2S|build|' id],120,@() dumps(c));
            item(G,['B2S|mosek|' id],60,@() bl_bisect(dmp(G,id),ms(struct('mode','feas','tag','B2S'))));
            item(G,['B2S|cuadmm|' id],7*60,@() bl_bisect(dmp(G,id),cu(struct('mode','feas', ...
                 'probe_tol',1e-7,'probe_cap',50000,'run_timeout',400,'tag','B2S'))));
            R = getres(G,['B2S|cuadmm|' id]);
            if ~isempty(R) && R.counts.F > 0
                fid = fopen(fullfile(RG,'SOUNDNESS_TRIP'),'w');
                fprintf(fid,'cuADMM F on sentinel %s (frac %.6f)\n',id,fr(q));  fclose(fid);
            end
        end
    end

case {'B1c','B1d'}
    if strcmp(block,'B1c'), L1 = {'hinf_rd1_hv','hinfdu_rd1','hinfco_rd1'};
    else,                   L1 = {'ctrl_rd1','hinfduco_rd1','est_rd1_pin'}; end
    for k = 1:numel(L1)
        nm = L1{k};  id = strrep(nm,'_pin','');  pin = 1;  if ~strcmp(nm,id), pin = 'auto'; end
        Rm = getres(G,['B1a|mosek|' nm]);
        if isempty(Rm) || ~isfinite(Rm.gamma_F)
            skiprow(G,[block '|cuadmm|' nm],'no Mosek-certified upper end in B1a');  continue
        end
        if strcmp(block,'B1d') && ~(Rm.gamma_I > 0)
            skiprow(G,[block '|cuadmm|' nm],'B1a Mosek bracket has no certified lower end');  continue
        end
        item(G,[block '|cuadmm|' nm],21*60,@() bl_bisect(dmp(G,id),cu(struct('hi',1.2*Rm.gamma_F, ...
             'cert_iter',30000,'max_wall',20*60,'run_timeout',300,'pin_scale',pin,'tag',block))));
        if Rm.gamma_I > 0
            item(G,[block '|heldout|' nm],3*60,@() bl_bisect(dmp(G,id),cu(struct('mode','probe', ...
                 'gammas',Rm.gamma_I,'probe_tol',1e-7,'probe_cap',30000,'run_timeout',300, ...
                 'pin_scale',pin,'tag',[block '_heldout']))));
        end
    end

case 'B1e'
    % CC, 09/26/2026: follow-ups to three B1c/B1d outcomes, measured this run:
    %   ctrl_rd1 was run before the lean_kmin fix and with the unscaled pin (its
    %     cuADMM bound is only the seed, 4700 vs Mosek 3917);
    %   hinfco_rd1 certified nothing at 30k iterations (repaired lambda_min
    %     -3e-7..-7e-7 relative) -- does a long, tight run certify at all?
    %   hinfdu_rd1 came out 1.7% loose because near-boundary feasible probes did
    %     not reach tau inside the 20k cap -- do larger caps close it?
    Rm = getres(G,'B1a|mosek|ctrl_rd1');
    if ~isempty(Rm) && isfinite(Rm.gamma_F)
        item(G,'B1e|cuadmm|ctrl_rd1',26*60,@() bl_bisect(dmp(G,'ctrl_rd1'),cu(struct('hi',1.2*Rm.gamma_F, ...
             'cert_iter',30000,'max_wall',25*60,'run_timeout',300,'pin_scale','auto','tag','B1e'))));
    end
    Rm = getres(G,'B1a|mosek|hinfco_rd1');
    if ~isempty(Rm) && isfinite(Rm.gamma_F)
        item(G,'B1e|probe|hinfco_rd1',6*60,@() bl_bisect(dmp(G,'hinfco_rd1'),cu(struct('mode','probe', ...
             'gammas',1.2*Rm.gamma_F,'probe_tol',1e-8,'probe_cap',100000,'run_timeout',330,'tag','B1e'))));
    end
    Rm = getres(G,'B1a|mosek|hinfdu_rd1');
    if ~isempty(Rm) && isfinite(Rm.gamma_F)
        item(G,'B1e|cuadmm|hinfdu_rd1',21*60,@() bl_bisect(dmp(G,'hinfdu_rd1'),cu(struct('hi',1.2*Rm.gamma_F, ...
             'cert_iter',60000,'kmax',60000,'max_wall',20*60,'run_timeout',300,'tag','B1e'))));
    end

case 'B2'
    for k = 1:numel(FEAS)
        id = FEAS{k};
        item(G,['B2|mosek|' id],60,@() bl_bisect(dmp(G,id),ms(struct('mode','feas','tag','B2'))));
        item(G,['B2|cuadmm|' id],7*60,@() bl_bisect(dmp(G,id),cu(struct('mode','feas', ...
             'probe_tol',1e-7,'probe_cap',50000,'run_timeout',400,'tag','B2'))));
    end

case 'B3b'
    hi = hi_of(G,'nl_fisher_opt');
    item(G,'B3b|mosek|nl_fisher_opt',40*60,@() bl_bisect(dmp(G,'nl_fisher_opt'),ms(struct('hi',hi, ...
         'rtol',1e-4,'max_probes',40,'hi_doublings',4,'max_wall',40*60,'tag','B3b'))));

case 'B3a'
    Rm = getres(G,'B3b|mosek|nl_fisher_opt');
    if isempty(Rm) || ~isfinite(Rm.gamma_F), skiprow(G,'B3a','no Mosek bracket from B3b'); return; end
    tl = [1e-4 1e-5 1e-6];  cp = [8000 15000 25000];
    for q = 1:3
        item(G,sprintf('B3a|cuadmm|tol%g',tl(q)),19*60,@() bl_bisect(dmp(G,'nl_fisher_opt'), ...
             cu(struct('mode','probe','gammas',1.2*Rm.gamma_F,'probe_tol',tl(q),'probe_cap',cp(q), ...
             'run_timeout',18*60,'tag',sprintf('B3a_%g',tl(q))))));
    end
    if Rm.gamma_I > 0
        item(G,'B3a|cuadmm|control',19*60,@() bl_bisect(dmp(G,'nl_fisher_opt'),cu(struct( ...
             'mode','probe','gammas',Rm.gamma_I,'probe_tol',1e-5,'probe_cap',15000, ...
             'run_timeout',18*60,'tag','B3a_control'))));
    end

case 'B3c'
    Rm = getres(G,'B3b|mosek|nl_fisher_opt');  tstar = NaN;  spi = NaN;
    for t = [1e-4 1e-5 1e-6]
        Ra = getres(G,sprintf('B3a|cuadmm|tol%g',t));
        if ~isempty(Ra) && Ra.counts.F > 0 && isnan(tstar), tstar = t;  spi = Ra.spi; end
    end
    Rc = getres(G,'B3a|cuadmm|control');
    if isnan(tstar) || (~isempty(Rc) && Rc.counts.F > 0)
        skiprow(G,'B3c','gate failed: no certified cuADMM point in B3a, or the control certified');  return
    end
    ci = floor(17*60/max(spi,1e-3));
    item(G,'B3c|cuadmm|nl_fisher_opt',76*60,@() bl_bisect(dmp(G,'nl_fisher_opt'),cu(struct( ...
         'hi',1.2*Rm.gamma_F,'hi_doublings',1,'cert_tol',tstar,'cert_iter',ci,'rtol',1e-3, ...
         'max_wall',75*60,'run_timeout',18*60,'pin_scale','auto','tag','B3c'))));    % CC, 09/26/2026: pin 'auto' (ctrl_rd1: unscaled pin saturates b over a wide bracket)

case 'B5'
    item(G,'B5|build|stab2_rd_psz',8*60,@() dumps({'stab2_rd_psz'}));
    item(G,'B5|mosek|stab2_rd_psz',10*60,@() bl_bisect(dmp(G,'stab2_rd_psz'),ms(struct('mode','feas','tag','B5'))));
    item(G,'B5|cuadmm|stab2_rd_psz',16*60,@() bl_bisect(dmp(G,'stab2_rd_psz'),cu(struct('mode','feas', ...
         'probe_tol',1e-7,'probe_cap',10000,'run_timeout',15*60,'tag','B5'))));

case 'B4'
    for n = [2 4 8]
        id = sprintf('scale_hinf_n%02d',n);
        item(G,['B4|build|' id],10*60,@() dumps({id}));
        hi = hi_of(G,id);
        item(G,['B4|mosek|' id],20*60,@() bl_bisect(dmp(G,id),ms(struct('hi',hi,'rtol',1e-3, ...
             'max_probes',30,'hi_doublings',4,'max_wall',20*60,'tag','B4'))));
        Rm = getres(G,['B4|mosek|' id]);
        if ~isempty(Rm) && isfinite(Rm.gamma_F)
            item(G,['B4|cuadmm|' id],5*60,@() bl_bisect(dmp(G,id),cu(struct('mode','probe', ...
                 'gammas',1.2*Rm.gamma_F,'probe_tol',1e-12,'probe_cap',2000,'run_timeout',240,'tag','B4'))));
        end
    end

% CC, 09/27/2026 (start): Sol-readiness blocks (header).
case 'N1'
    chk = {};  ck = @(name,ok) sprintf('%s %s',passfail(ok),name);
    item(G,'N1|ok',5,@() struct('probe',1));
    n0 = nrows(G.MAN,'N1|ok');  item(G,'N1|ok',5,@() struct('probe',2));
    chk{end+1} = ck('resume: a done item is not re-run',nrows(G.MAN,'N1|ok') == n0);
    G2 = G;  G2.DL = posixtime(datetime('now')) + 30;
    item(G2,'N1|deadline',600,@() struct('probe',3));
    chk{end+1} = ck('deadline: an item that would overrun is refused',lastnote(G.MAN,'N1|deadline',"SKIP_DEADLINE"));
    sf = fullfile(G.OUT,'STOP');  fid = fopen(sf,'w');  fclose(fid);
    cl = onCleanup(@() delete_if(sf));         % never leave STOP behind
    item(G,'N1|stop',5,@() struct('probe',4));
    chk{end+1} = ck('STOP: the driver refuses a new item',lastnote(G.MAN,'N1|stop',"STOP file"));
    try
        Rs = bl_bisect(dmp(G,'stab_tr1'),cu(struct('mode','feas','probe_tol',1e-7, ...
             'probe_cap',2000,'run_timeout',60,'stopfile',sf,'tag','N1stop')));
        nolog = isempty(dir(fullfile(G.BO,'stab_tr1_N1stop_cuadmm','run_*.log')));
        chk{end+1} = ck('STOP: bl_bisect makes no solver call',~isempty(Rs.stopped) && nolog);
    catch ME
        chk{end+1} = ck(['STOP: bl_bisect (' ME.message ')'],false);
    end
    delete_if(sf);
    fid = fopen(fullfile(G.RG,'N1_checks.txt'),'w');  fprintf(fid,'%s\n',chk{:});  fclose(fid);
    fprintf('%s\n',chk{:});
    if any(startsWith(chk,'FAIL')), fid = fopen(fullfile(G.RG,'DRIVER_TRIP'),'w'); fprintf(fid,'%s\n',chk{:}); fclose(fid); end

case 'N2b'
    for id = {'hinf_rd1','hinfdu_rd1'}
        Rm = getres(G,['B1a|mosek|' id{1}]);
        if isempty(Rm) || ~(Rm.gamma_I > 0), skiprow(G,['N2b|' id{1}],'no Mosek-certified gamma_I'); continue; end
        for fr = [0.99 0.999]
            hi = fr*Rm.gamma_I;  t = strrep(sprintf('N2b_%g',fr),'.','p');
            item(G,sprintf('N2b|sedumi|%s|%g',id{1},fr),60,@() bl_bisect(dmp(G,id{1}), ...
                 sd(struct('hi',hi,'hi_doublings',0,'rtol',1e-3,'max_probes',3,'tag',t))));
            item(G,sprintf('N2b|cuadmm|%s|%g',id{1},fr),6*60,@() bl_bisect(dmp(G,id{1}), ...
                 cu(struct('hi',hi,'hi_doublings',0,'cert_iter',30000,'run_timeout',300, ...
                 'max_probes',3,'tag',t))));
            for sv = {'sedumi','cuadmm'}
                R = getres(G,sprintf('N2b|%s|%s|%g',sv{1},id{1},fr));
                if ~isempty(R) && R.counts.F > 0, trip(G,sprintf('%s F at %g x Mosek gamma_I of %s',sv{1},fr,id{1})); end
            end
        end
    end

case 'S1'
    % worst-case caps priced from the 09-24 ms/it (contended, so conservative):
    % n08 15, n16 21, n24 38, n32 55 ms/it
    spi = containers.Map({8,16,24,32},{0.015,0.021,0.038,0.055});
    cap = 30000;
    for n = [8 16 24 32]
        id = sprintf('scale_stab_n%02d',n);  rt = ceil(1.3*cap*spi(n)) + 120;
        item(G,['S1|cuadmm|' id],rt + 300,@() bl_bisect(dmp(G,id),cu(struct('mode','feas', ...
             'probe_tol',1e-7,'probe_cap',cap,'run_timeout',rt,'tag','S1'))));
    end
    item(G,'S1|build|sentinels',15*60,@() bigbuild([16 24 32],1.01));
    for n = [16 24 32]
        id = sprintf('scale_sent_f1p01_n%02d',n);  rt = ceil(1.3*cap*spi(n)) + 120;
        item(G,['S1|cuadmm|' id],rt + 300,@() bl_bisect(dmp(G,id),cu(struct('mode','feas', ...
             'probe_tol',1e-7,'probe_cap',cap,'run_timeout',rt,'tag','S1'))));
        R = getres(G,['S1|cuadmm|' id]);
        if ~isempty(R) && R.counts.F > 0, trip(G,sprintf('cuADMM F on sentinel %s (1.01 lambda*)',id)); end
    end

case 'N4'
    xk = fullfile(G.OUT,'keep','Th_n3_X_1e-6.txt');
    if exist(xk,'file')
        item(G,'N4|file|Th_n3',45*60,@() bl_bisect(dmp(G,'Th_n3'),struct('mode','feas','solver','file', ...
             'xfile',xk,'outdir',G.BO,'deadline',DL,'tag','N4f')));
    else
        skiprow(G,'N4|file|Th_n3',['no kept point ' xk]);
    end
    item(G,'N4|cuadmm|Th_n3',60*60,@() bl_bisect(dmp(G,'Th_n3'),cu(struct('mode','feas', ...
         'probe_tol',1e-4,'probe_cap',4000,'run_timeout',35*60,'tag','N4'))));

case 'X0'
    % Sol-0 smoke, run identically on the workstation and on Sol: a discarded
    % warm-up (GPU clock ramp: a first timed run can be ~2.6x slow), then three
    % certified-standard feasibility runs -- the cheapest case, the ladder
    % anchor, and the coupled ladder at n = 4 (bl_big(4,0.5,true)).
    ff = @(t) struct('mode','feas','probe_tol',1e-7,'probe_cap',30000,'run_timeout',240,'tag',t);
    item(G,'X0|warm|stab_tr1',120,@() bl_bisect(dmp(G,'stab_tr1'),cu(ff('X0w'))));
    for id = {'stab_tr1','stab_rd1_hv','scale_rot_f0p50_n04'}
        item(G,['X0|cuadmm|' id{1}],300,@() bl_bisect(dmp(G,id{1}),cu(ff('X0'))));
    end

case 'X1'
    % Sol vs desktop at large m (CC, 09/28/2026): a warm-up, then scale_stab_n32
    % (m 138,240) exactly as S1 ran it on the workstation 09-27: 14,092 it,
    % 354.8 s solve, 25.2 ms/it, init 0.4 s, certify 4.2 s, F.
    item(G,'X1|warm|stab_tr1',120,@() bl_bisect(dmp(G,'stab_tr1'),cu(struct('mode','feas', ...
         'probe_tol',1e-7,'probe_cap',30000,'run_timeout',240,'tag','X1w'))));
    item(G,'X1|cuadmm|scale_stab_n32',10*60,@() bl_bisect(dmp(G,'scale_stab_n32'),cu(struct('mode','feas', ...
         'probe_tol',1e-7,'probe_cap',30000,'run_timeout',15*60,'tag','X1'))));

case 'S2'
    % the COUPLED ladder (bl_big rot = true, header there): a large problem that
    % is not decoupled copies.  Measured 09-27 at 1e-6: 3937 / 3261 / 6519
    % iterations at n = 2 / 4 / 8 against 4023 at every n for scale_stab, and
    % Mosek-F at all three (the invariance argument, checked at small n).
    % Order: truth checks and the n16 sentinel first, so the deadline cuts the
    % largest rungs, not the soundness items.  Caps priced at ~25/45/70 ms/it
    % for n16/24/32 (INFERRED: n08 measured 9 ms/it, decoupled 21/38/55).
    item(G,'S2|build|rot05',20*60,@() bigbuild_rot([8 16 24 32],0.5));
    item(G,'S2|build|rot101',20*60,@() bigbuild_rot([16 24 32],1.01));
    for n = [8 16]
        id = sprintf('scale_rot_f0p50_n%02d',n);
        item(G,['S2|mosek|' id],15*60,@() bl_bisect(dmp(G,id),ms(struct('mode','feas','tag','S2'))));
    end
    cap = 40000;  spi = containers.Map({8,16,24,32},{0.012,0.025,0.045,0.070});
    L2 = {'scale_rot_f0p50_n08','scale_rot_f0p50_n16','scale_rot_f1p01_n16', ...
          'scale_rot_f0p50_n24','scale_rot_f1p01_n24','scale_rot_f0p50_n32','scale_rot_f1p01_n32'};
    for k = 1:numel(L2)
        id = L2{k};  n = str2double(id(end-1:end));  rt = ceil(1.3*cap*spi(n)) + 120;
        item(G,['S2|cuadmm|' id],rt + 300,@() bl_bisect(dmp(G,id),cu(struct('mode','feas', ...
             'probe_tol',1e-7,'probe_cap',cap,'run_timeout',rt,'tag','S2'))));
        R = getres(G,['S2|cuadmm|' id]);
        if contains(id,'f1p01') && ~isempty(R) && R.counts.F > 0
            trip(G,sprintf('cuADMM F on coupled sentinel %s (1.01 lambda*)',id));
        end
    end
% CC, 09/27/2026 (end)

otherwise
    error('bl_regime:block','unknown block %s',block);
end
fprintf('REGIMEDONE %s %s\n',block,datestr(now,31));
end


% CC, 09/27/2026 (start): helpers for the Sol-readiness blocks.
function n = nrows(MAN,key)
% manifest rows with this key (resume test: a skipped item adds none)
n = 0;  fid = fopen(MAN,'r');  if fid < 0, return; end
fgetl(fid);
while true
    l = fgetl(fid);  if ~ischar(l), break; end
    if startsWith(l,[key sprintf('\t')]), n = n + 1; end
end
fclose(fid);
end

function ok = lastnote(MAN,key,pat)
% the last row for key is a SKIP whose note contains pat
ok = false;  fid = fopen(MAN,'r');  if fid < 0, return; end
fgetl(fid);  last = '';
while true
    l = fgetl(fid);  if ~ischar(l), break; end
    if startsWith(l,[key sprintf('\t')]), last = l; end
end
fclose(fid);
ok = contains(last,sprintf('\tSKIP\t')) && contains(last,pat);
end

function delete_if(f), if exist(f,'file'), delete(f); end, end

function s = passfail(ok), if ok, s = 'PASS'; else, s = 'FAIL'; end, end

function trip(G,msg)
fid = fopen(fullfile(G.RG,'SOUNDNESS_TRIP'),'a');  fprintf(fid,'%s %s\n',datestr(now,31),msg);  fclose(fid);
fprintf('REGIME SOUNDNESS_TRIP %s\n',msg);
end

function R = bigbuild(ns,frac)
bl_big(ns,frac);
R = struct('built',{arrayfun(@(n) sprintf('n%02d',n),ns,'UniformOutput',false)});
end

function R = bigbuild_rot(ns,frac)
bl_big(ns,frac,true);
R = struct('built',{arrayfun(@(n) sprintf('rot n%02d',n),ns,'UniformOutput',false)});
end
% CC, 09/27/2026 (end)


% =========================================================================
function item(G,key,est_s,fn)
% run one item unless done, or unless it would end past the deadline
if isdone(G.MAN,key), fprintf('REGIME skip %s (done)\n',key); return; end
if exist(fullfile(G.OUT,'STOP'),'file'), skiprow(G,key,'STOP file'); return; end
if posixtime(datetime('now')) + est_s > G.DL, skiprow(G,key,sprintf('SKIP_DEADLINE (est %.0f s)',est_s)); return; end
fprintf('REGIME start %s (est %.0f s) %s\n',key,est_s,datestr(now,31));
st = datestr(now,31);  t = tic;  status = 'ok';  note = '';  R = [];
try
    R = fn();
catch ME
    status = 'ERR';  note = regexprep(ME.message,'[\t\r\n]+',' ');
    fprintf('REGIME ERR %s: %s\n',key,note);
end
w = toc(t);
if ~isempty(R), save(fullfile(G.RG,[fsafe(key) '.mat']),'R'); end
row = {key,status,st,sprintf('%.1f',w)};
if isstruct(R) && isfield(R,'gamma_F')
    row = [row, cellfun(@(v) sprintf('%.10g',v),{R.gamma_I,R.gamma_F,R.gamma_i,R.gamma_f},'UniformOutput',false), ...
           cellfun(@num2str,{R.counts.F,R.counts.f,R.counts.i,R.counts.I},'UniformOutput',false), ...
           {num2str(R.iters),sprintf('%.1f',R.t_solver),strjoin(R.flags,' | '),note}];
else
    row = [row, repmat({''},1,11), {note}];
end
fid = fopen(G.MAN,'a');  fprintf(fid,'%s\n',strjoin(row,sprintf('\t')));  fclose(fid);
fprintf('REGIME done %s %s %.1f s\n',key,status,w);
end


function skiprow(G,key,why)
fid = fopen(G.MAN,'a');
fprintf(fid,'%s\tSKIP\t%s\t0%s%s\n',key,datestr(now,31),repmat(sprintf('\t'),1,12),why);  % why -> 'note' (col 16)
fclose(fid);
fprintf('REGIME skip %s: %s\n',key,why);
end


function d = isdone(MAN,key)
d = false;  fid = fopen(MAN,'r');  if fid < 0, return; end
fgetl(fid);
while true
    l = fgetl(fid);  if ~ischar(l), break; end
    p = strsplit(l,sprintf('\t'),'CollapseDelimiters',false);
    if numel(p) >= 2 && strcmp(p{1},key) && strcmp(p{2},'ok'), d = true; end
end
fclose(fid);
end


function R = getres(G,key)
f = fullfile(G.RG,[fsafe(key) '.mat']);  R = [];
if exist(f,'file'), S = load(f);  R = S.R; end
end


function R = dumps(ids)
% (re)build and dump through bl_run; returns a small struct for the manifest
bl_run(ids);
R = struct('built',{ids});
end


function f = dmp(G,id), f = fullfile(G.DMP,[id '.mat']); end


function hi = hi_of(G,id)
% 1.2 x the objective-form gamma from this run's bl_mosek.tsv; for the cases
% whose objective form returns nonsense (hinfduco 8251.9, ctrl 15433, measured
% 09-24) start low and let the certified seed double upward instead.
switch id
    case 'hinfduco_rd1', hi = 0.3;  return
    case 'ctrl_rd1',     hi = 1;    return
end
T = fullfile(G.OUT,'baseline','bl_mosek.tsv');
fid = fopen(T,'r');  h = strsplit(fgetl(fid),sprintf('\t'));  jg = find(strcmp(h,'gam'));  hi = NaN;
while true
    l = fgetl(fid);  if ~ischar(l), break; end
    p = strsplit(l,sprintf('\t'),'CollapseDelimiters',false);
    if strcmp(p{1},id), hi = 1.2*str2double(p{jg}); end
end
fclose(fid);
% CC, 09/26/2026: no gam recorded (the PIESOS builders return none) -> solve the
% dump's objective form once with Mosek; gamma = c'x * bscl, since bl_fix pins
% x_j = gamma in the UNnormalised units (b_un = b*bscl).  Measured need:
% nl_fisher_opt's bl_mosek row has an empty gam, which ERRed B3b.
if ~(hi > 0)
    try
        S  = load(dmp(G,id));  Mt = load(strrep(dmp(G,id),'.mat','_meta.mat'));
        [~,res] = mosekopt('minimize info echo(0)',Sedumi2Mosek(S.At',full(S.b),S.c,S.K));
        x  = MosekSol2SedumiSol(S.K,res);
        hi = 1.2*full(S.c(:)'*x(:))*Mt.bscl;
        fprintf('REGIME hi_of %s: objective form solved, gamma_obj*1.2 = %.8g (%s)\n',id,hi,res.sol.itr.prosta);
    catch ME
        fprintf('REGIME hi_of %s: objective solve failed: %s\n',id,ME.message);
    end
end
% NaN, not an error: this runs outside item(), so an error here would end the
% whole block; bl_bisect rejects hi = NaN inside the item and logs it ERR
if ~(hi > 0), warning('bl_regime:hi','no objective-form gamma for %s in %s',id,T); hi = NaN; end
end


function s = merge(a,b)
s = a;  f = fieldnames(b);
for i = 1:numel(f), s.(f{i}) = b.(f{i}); end
end


function s = fsafe(k), s = regexprep(k,'[^A-Za-z0-9_]','_'); end


function s = stamp(p)
if isfinite(p), s = datestr(datetime(p,'ConvertFrom','posixtime','TimeZone','local'),31); else, s = 'none'; end
end
