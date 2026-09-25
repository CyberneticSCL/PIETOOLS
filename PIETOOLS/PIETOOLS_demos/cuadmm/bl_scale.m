function bl_scale()
% bl_scale -- the size ladder, with the solve GATED ON MEASURED SHAPE.
%
% This is the arm that tests the only claim that really matters for cuADMM:
% that it keeps going where an interior-point method runs out of memory. n
% decoupled copies of the same state leave the physics and the analysis
% untouched and grow only the SDP, so every rung is the same problem at a
% different size rather than a different problem.
%
% WHY A SEPARATE RUNNER. bl_run's builders call the executive, which solves
% internally, so m is unknown until the solve is already paid for. Mosek's Schur
% complement is dense m x m: at m = 138,000 that is 152 GB against this machine's
% 64 GB, so an ungated rung would not fail cleanly -- it would thrash the page
% file for hours in an unattended window and take the rest of the suite with it.
% stab_mirror builds the identical program WITHOUT solving (faithfulness checked
% against the stock executive: matching m, Kf, Ns, nnz), so the shape can be
% measured and the rung gated before anything expensive happens.
%
% Three outcomes are recorded and they are different findings:
%   ok        solved
%   skip_big  built, shape recorded, solve declined by the gate. The SHAPE is
%             still data -- it is how the ladder's m(n) curve is extended past
%             the point where this machine can solve.
%   ERR       failed to build

cuadmm_path;
HERE = fileparts(mfilename('fullpath'));
OUT  = fullfile(cuadmm_outdir(),'baseline');  if ~exist(OUT,'dir'), mkdir(OUT); end
DMP  = fullfile(OUT,'dumps');      if ~exist(DMP,'dir'), mkdir(DMP); end
TSV  = fullfile(OUT,'bl_scale.tsv');

MAXM_SOLVE = 0;        % SHAPE-ONLY pass: n<=16 was already solved cleanly by the
                       % ungated run (m=135n^2 exactly, 0.037/0.200/3.82/117.8 s).
                       % Mosek hit 22.6 GB RSS on n=24 (m=77760) with 18.6 GB free and
                       % was stopped before it thrashed; its dense Schur needs 48 GB on a
                       % 64 GB machine. So the large rungs contribute SHAPE, not timing.
MAXM_DUMP  = 0;
NS = [1 2 3 4 6 8 12 16 24 32];

COLS = {'id','n','status','m','Kf','nblk','Ks','nnzAt','nvar','vec_len','normb', ...
        't_build','t_mosek','rel_b','psd_min','psd_relmin','numerr','normx','note'};
if ~exist(TSV,'file')
    fid=fopen(TSV,'w'); fprintf(fid,'%s\n',strjoin(COLS,sprintf('\t'))); fclose(fid);
end
done = readdone(TSV);

for n = NS
    id = sprintf('scale_stab_n%02d',n);
    if isKey(done,id), fprintf('SC skip %s\n',id); continue; end
    row = emptyrow(COLS); row.id=id; row.n=n;
    fprintf('SC start %s\n',id);
    try
        clear stateNameGenerator
        pvar s t
        x = pde_var(n,s,[0,1]);
        PIE = initialize(convert([diff(x,t,1)==diff(x,s,2)+0.5*pi^2*x; ...
                                  subs(x,s,0)==0; subs(x,s,1)==0]));
        st = lpisettings('heavy');
        st.sos_opts.solver='mosek'; st.sos_opts.simplify=false;
        tb = tic;  evalc('out = stab_mirror(PIE,st);');  row.t_build = toc(tb);

        S = sdpshape(out.prog);
        row.m=S.m; row.Kf=S.Kf; row.nblk=S.nblk; row.Ks=mat2str(S.Ks);
        row.nnzAt=S.nnzAt;
        row.nvar = S.Kf + sum(double(S.Ks).^2);
        row.vec_len = S.Kf + sum(double(S.Ks).*(double(S.Ks)+1)/2);
        bf=[]; for q=1:out.prog.expr.num, bf=[bf;out.prog.expr.b{q}]; end %#ok<AGROW>
        row.normb = norm(full(bf));

        if S.m > MAXM_SOLVE
            row.status = 'skip_big';
            row.note = sprintf('m=%d over solve gate %d (dense Schur %.0f GB)', ...
                               S.m,MAXM_SOLVE,8*S.m^2/2^30);
        else
            ts = tic;  evalc('sol = lpisolve(out.prog,st.sos_opts);');  row.t_mosek = toc(ts);
            Atf=[]; bf2=[];
            for q=1:sol.expr.num, Atf=[Atf,sol.expr.At{q}]; bf2=[bf2;sol.expr.b{q}]; end %#ok<AGROW>
            xr = sol.solinfo.RRx(:);  xc = sol.solinfo.x(:);
            row.rel_b = norm(full(Atf'*xr-bf2))/max(norm(full(bf2)),eps);
            row.normx = norm(xr);
            [row.psd_min,row.psd_relmin] = psdmargin(xc,S);   % CONE vector, not RRx
            I = sol.solinfo.info;  row.numerr = getf(I,'numerr');
            if isfield(I,'cpusec'), row.t_mosek = I.cpusec; end
            if S.m <= MAXM_DUMP
                cvec = sparse(size(Atf,1),1);
                RR = mkRR(sol);  bscl = norm(full(bf2)); if bscl==0, bscl=1; end
                D.At = RR'*Atf;  D.b = bf2/bscl;  D.c = RR'*cvec;
                D.K = struct('f',S.Kf,'l',0,'q',[],'s',S.Ks);
                D.Ns = S.Ks;  D.Kf = S.Kf;
                save(fullfile(DMP,[id '.mat']),'-struct','D','-v7.3');
                save(fullfile(DMP,[id '_meta.mat']),'RR','bscl','S');
                dump2cuadmm(fullfile(DMP,[id '.mat']),fullfile(DMP,id));
            end
            row.status = 'ok';
        end
    catch ME
        row.status='ERR'; row.note=regexprep(ME.message,'[\t\r\n]+',' ');
    end
    appendrow(TSV,COLS,row);
    fprintf('SC %-18s %-9s n=%-3d m=%-8s nvar=%-9s t=%-9s rel=%-11s %s\n', ...
        row.id,row.status,n,num2str(row.m),num2str(row.nvar), ...
        num2str(row.t_mosek),num2str(row.rel_b),row.note);
end
fprintf('SCDONE\n');
end

function [pmin,prel] = psdmargin(x,S)
off = S.Kf;  pmin = inf;  pmax = -inf;
for k = 1:numel(S.Ks)
    N = double(S.Ks(k));
    if off+N^2 > numel(x), pmin=NaN; prel=NaN; return; end
    Xk = reshape(x(off+(1:N^2)),N,N);  ev = eig((Xk+Xk')/2);
    pmin = min(pmin,min(ev));  pmax = max(pmax,max(ev));  off = off + N^2;
end
prel = pmin/max(pmax,eps);
end

function RR = mkRR(prog)
RR = speye(prog.var.idx{1}-1);
for i=1:prog.var.num
    sz = prog.var.idx{i+1}-prog.var.idx{i};
    switch prog.var.type{i}
        case 'poly', RR = spantiblkdiag(RR,speye(sz));
        case 'sos',  RR = spblkdiag(RR,speye(sz));
        otherwise,   error('bl_scale:vartype','unexpected var type %s',prog.var.type{i});
    end
end
for i=1:prog.extravar.num
    RR = spblkdiag(RR,speye(prog.extravar.idx{i+1}-prog.extravar.idx{i}));
end
end

function v = getf(I,f), if isstruct(I)&&isfield(I,f), v=I.(f); else, v=NaN; end, end

function r = emptyrow(COLS)
r = struct(); for i=1:numel(COLS), r.(COLS{i}) = NaN; end
r.id=''; r.status=''; r.Ks=''; r.note='';
end

function appendrow(TSV,COLS,row)
fid=fopen(TSV,'a'); s=cell(1,numel(COLS));
for i=1:numel(COLS)
    v=row.(COLS{i});
    if ischar(v), s{i}=v; elseif isempty(v), s{i}=''; else, s{i}=num2str(v,'%.6g'); end
end
fprintf(fid,'%s\n',strjoin(s,sprintf('\t'))); fclose(fid);
end

function d = readdone(TSV)
d = containers.Map('KeyType','char','ValueType','logical');
fid=fopen(TSV,'r'); fgetl(fid);
while true
    l=fgetl(fid); if ~ischar(l), break; end
    p=strsplit(l,sprintf('\t'),'CollapseDelimiters',false);
    if numel(p)>=3 && ~isempty(p{1}) && any(strcmp(p{3},{'ok','skip_big'})), d(p{1})=true; end
end
fclose(fid);
end
