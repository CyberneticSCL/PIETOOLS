function bl_run(chunk)
% bl_run(chunk) -- Mosek arm of the overnight baseline, one chunk of classes.
%
% Writes ONE TSV line per case, flushed immediately, so a hard MATLAB crash
% loses only the case in flight. Re-running skips ids already recorded 'ok', so
% the suite is resumable with no bookkeeping. Chunking by class keeps the ~45 s
% of startup+path cost to a few payments instead of one per case, at the price
% of losing the rest of one chunk to a crash.
%
% THREE TIMES ARE RECORDED, because they measure different things and only one
% of them is comparable to the GPU arm:
%   t_wall    wall clock around the whole builder: PDE->PIE conversion, LPI
%             construction AND solve. What a user actually waits for.
%   t_mosek   Mosek's OWN optimizer clock, MSK_DINF_OPTIMIZER_TIME, which
%             sossolve.m:367 stores as solinfo.info.cpusec. Excludes assembly.
%   t_dump    Mosek re-solving the dumped bytes directly via Sedumi2Mosek.
%             THIS is the like-for-like number: cuADMM is handed exactly these
%             At/b/c/K, so this and CUADMM_TIMING solve measure the same work on
%             the same data. t_mosek and t_dump differing would itself be a
%             finding -- it would mean the dump is not the program that was
%             solved.
%
% CAPTURE IS POST-SOLVE, deliberately. Every executive calls lpisolve
% internally, and for programs carrying an 'ineq' expression sossolve inserts
% blocks via addextrasosvar DURING the solve (sossolve.m:177-183). A pre-solve
% capture of an Hinf program therefore misses the gamma>=0 block -- measured
% earlier as gamma = -0.0 with numerr 1. Reading sol.expr AFTER the solve gets
% the complete system.
%
% NO ACCEPT/REJECT IS STORED, only quantities. Three separate gates have already
% manufactured a wrong answer in this campaign, so the file records enough to
% diagnose a bad verdict afterwards WITHOUT re-solving: normx together with
% rel_b exposes the X=0 trap (rel_b is exactly 1 there), feasratio is recorded
% but never acted on (meaningless at c=0: measured -0.0612 at a feasible point,
% +0.956 with numerr=2 at an infeasible one), and the PSD margin is normalised
% by the largest eigenvalue ACROSS blocks, since a per-block ratio reads -1 for
% a numerically zero block whatever the solution.

cuadmm_path;
HERE = fileparts(mfilename('fullpath'));
OUT  = fullfile(cuadmm_outdir(),'baseline');  if ~exist(OUT,'dir'), mkdir(OUT); end
DMP  = fullfile(OUT,'dumps');      if ~exist(DMP,'dir'), mkdir(DMP); end
TSV  = fullfile(OUT,'bl_mosek.tsv');

COLS = {'id','cls','dim','kind','status','m','Kf','nblk','Ks','nnzAt','nvar', ...
        'vec_len','normb','c_nnz','t_wall','t_mosek','t_dump','rel_b','rel_dump', ...
        'gam','gam_exact','feasratio','numerr','pinf','dinf','psd_min','psd_relmin', ...
        'normx','trivial','psd_rel_bad','note'};
if ~exist(TSV,'file')
    fid=fopen(TSV,'w'); fprintf(fid,'%s\n',strjoin(COLS,sprintf('\t'))); fclose(fid);
end
done = readdone(TSV);
% Resume state lives under cuadmm_outdir(), which PERSISTS across sessions and
% checkouts. Skipping silently would mix rows measured on old code with new ones
% after a change, so say how many will be reused and how to start fresh.
if done.Count > 0
    warning('bl_run:resume', ['RESUMING: %d case(s) already recorded in\n  %s\n' ...
        'will be SKIPPED, not re-measured. If the code under test changed since, ' ...
        'delete that file or setenv(''CUADMM_OUT'',<new folder>) first.'], done.Count, TSV);
end

C = bl_cases();
for i = 1:numel(C)
    c = C(i);
    if nargin>=1 && ~isempty(chunk) && ~any(strcmpi(c.cls,cellstr(chunk))), continue; end
    if isKey(done,c.id), fprintf('BL skip %s\n',c.id); continue; end
    row = emptyrow(COLS);
    row.id=c.id; row.cls=c.cls; row.dim=c.dim; row.kind=c.kind;
    % announce on START, not only on completion: a runner that prints only when a
    % case finishes tells you nothing about which one is in flight, and reading the
    % skip lines instead cost a wrong attribution of a 19.7 GB run to the Hinf class
    % when it was scale_stab_n24's Mosek solve.
    fprintf('BL start %s\n',c.id);
    try
        tw = tic;
        evalc('[sol,M] = feval(c.builder,c.args{:});');
        row.t_wall = toc(tw);

        Atf=[]; bf=[];
        for q=1:sol.expr.num, Atf=[Atf,sol.expr.At{q}]; bf=[bf;sol.expr.b{q}]; end %#ok<AGROW>
        x  = sol.solinfo.RRx(:);
        nb = norm(full(bf));
        row.rel_b = norm(full(Atf'*x-bf))/max(nb,eps);
        row.normx = norm(x);
        row.trivial = double((abs(row.rel_b-1)<=1e-6) || row.normx<=1e-12);

        S = sdpshape(sol);
        row.m=S.m; row.Kf=S.Kf; row.nblk=S.nblk; row.Ks=mat2str(S.Ks);
        row.nnzAt=S.nnzAt; row.normb=nb; row.c_nnz=S.c_nnz;
        row.nvar = S.Kf + sum(double(S.Ks).^2);
        row.vec_len = S.Kf + sum(double(S.Ks).*(double(S.Ks)+1)/2);
        % PSD MUST be read from solinfo.x (SeDuMi CONE order), not RRx (decvartable
        % order). MEASURED: slicing RRx by cone blocks gives relmin -8.96e-02 on
        % this plant, IDENTICAL at eppos2 = 1e-6,1e-4,1e-2,1 -- scale-invariant, so
        % it is a slicing artefact and not a noise floor. The cone vector gives
        % +1.7e-11. rel_b above is unaffected: RRx and Atf are both decvartable.
        xc = sol.solinfo.x(:);
        [row.psd_min,row.psd_relmin] = psdmargin(xc,S);
        [~,row.psd_rel_bad]         = psdmargin(x,S);   % the wrong-vector value, kept visible

        I = sol.solinfo.info;
        row.t_mosek = getf(I,'cpusec');  row.feasratio = getf(I,'feasratio');
        row.numerr  = getf(I,'numerr');  row.pinf = getf(I,'pinf'); row.dinf = getf(I,'dinf');
        if isfield(M,'gam') && isnumeric(M.gam) && isscalar(M.gam), row.gam = double(M.gam); end
        if isfield(M,'gam_exact'), row.gam_exact = M.gam_exact; end

        % ---- dump exactly what sossolve solved, then re-solve those bytes
        cvec = sparse(size(Atf,1),1);
        if S.c_nnz>0 && isfield(sol,'objective') && ~isempty(sol.objective)
            cvec(1:numel(sol.objective)) = sol.objective;
        end
        RR = mkRR(sol);  bscl = nb;  if bscl==0||~isfinite(bscl), bscl=1; end
        D.At = RR'*Atf;  D.b = bf/bscl;  D.c = RR'*cvec;
        D.K  = struct('f',S.Kf,'l',0,'q',[],'s',S.Ks);
        D.Ns = S.Ks;  D.Kf = S.Kf;      % dump2cuadmm reads these off the .mat
        save(fullfile(DMP,[c.id '.mat']),'-struct','D','-v7.3');
        save(fullfile(DMP,[c.id '_meta.mat']),'RR','bscl','S');
        dump2cuadmm(fullfile(DMP,[c.id '.mat']),fullfile(DMP,c.id));

        try
            prob = Sedumi2Mosek(D.At',full(D.b),D.c,D.K);
            [~,res] = mosekopt('minimize info echo(0)',prob);
            row.t_dump = res.info.MSK_DINF_OPTIMIZER_TIME;
            xd = MosekSol2SedumiSol(D.K,res);
            row.rel_dump = norm(full(D.At'*xd(:)-D.b))/max(norm(full(D.b)),eps);
        catch ME2
            row.note = ['dumpsolve: ' regexprep(ME2.message,'[\t\r\n]+',' ')];
        end
        row.status = 'ok';
    catch ME
        row.status = 'ERR';
        row.note = regexprep(ME.message,'[\t\r\n]+',' ');
    end
    appendrow(TSV,COLS,row);
    fprintf('BL %-18s %-4s m=%-7s tM=%-8s tD=%-8s rel=%-10s %s\n', row.id,row.status, ...
        num2str(row.m),num2str(row.t_mosek),num2str(row.t_dump),num2str(row.rel_b),row.note);
end
fprintf('BLDONE %s\n', strjoin(cellstr(chunk),','));
end


function [pmin,prel] = psdmargin(x,S)
% Normalise by the largest eigenvalue ACROSS ALL BLOCKS. A per-block ratio was
% measured to report -1.000 at every setting, because a numerically zero block's
% own maximum is rounding noise -- a property of the normalisation, not the fit.
off = S.Kf;  pmin = inf;  pmax = -inf;
for k = 1:numel(S.Ks)
    N = double(S.Ks(k));
    Xk = reshape(x(off+(1:N^2)),N,N);  ev = eig((Xk+Xk')/2);
    pmin = min(pmin,min(ev));  pmax = max(pmax,max(ev));  off = off + N^2;
end
prel = pmin/max(pmax,eps);
end

function RR = mkRR(prog)
% Replicates processvars (sossolve.m:1156-1190). The 'poly' arm is kept and the
% otherwise-error arm too: 'poly' is unreachable for the LPI executives but LIVE
% for polyopvar LocalStability, and it PREPENDS columns via spantiblkdiag, which
% silently invalidates any hand-computed index into x.
RR = speye(prog.var.idx{1}-1);
for i=1:prog.var.num
    sz = prog.var.idx{i+1}-prog.var.idx{i};
    switch prog.var.type{i}
        case 'poly', RR = spantiblkdiag(RR,speye(sz));
        case 'sos',  RR = spblkdiag(RR,speye(sz));
        otherwise,   error('bl_run:vartype','unexpected var type %s',prog.var.type{i});
    end
end
for i=1:prog.extravar.num
    RR = spblkdiag(RR,speye(prog.extravar.idx{i+1}-prog.extravar.idx{i}));
end
end

function v = getf(I,f), if isstruct(I) && isfield(I,f), v = I.(f); else, v = NaN; end, end

function r = emptyrow(COLS)
r = struct();
for i=1:numel(COLS), r.(COLS{i}) = NaN; end
r.id=''; r.cls=''; r.dim=''; r.kind=''; r.status=''; r.Ks=''; r.note='';
end

function appendrow(TSV,COLS,row)
fid = fopen(TSV,'a');  s = cell(1,numel(COLS));
for i=1:numel(COLS)
    v = row.(COLS{i});
    if ischar(v), s{i}=v; elseif isempty(v), s{i}=''; else, s{i}=num2str(v,'%.6g'); end
end
fprintf(fid,'%s\n',strjoin(s,sprintf('\t')));  fclose(fid);
end

function d = readdone(TSV)
d = containers.Map('KeyType','char','ValueType','logical');
fid = fopen(TSV,'r');  fgetl(fid);
while true
    l = fgetl(fid);  if ~ischar(l), break; end
    p = strsplit(l,sprintf('\t'),'CollapseDelimiters',false);
    if numel(p)>=5 && ~isempty(p{1}) && strcmp(p{5},'ok'), d(p{1}) = true; end
end
fclose(fid);
end
