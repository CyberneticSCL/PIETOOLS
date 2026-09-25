% bl_big.m -- the rungs Mosek cannot reach.
%
% Unlocked by the lpi_eq fix (nnz instead of ~all(all(C.C==0))): before it, the
% n=32 program could not even be BUILT -- construction requested a 9613344x1536
% dense array, 123.8 GB. The wall was in PIETOOLS, not in any solver.
%
% No Mosek arm here, deliberately. Its Schur complement is dense m x m: 45 GB at
% n=24 and 153 GB at n=32 against this machine's 64 GB. A run was started at
% n=24 and reached 22.6 GB resident with 18.6 GB free before being stopped, so
% the ladder's Mosek side ends at n=16 (m=34560, 117.8 s) by memory, not by
% patience. These rungs exist to ask whether cuADMM keeps going past that.
cuadmm_path;
HERE = fileparts(mfilename('fullpath'));
DMP  = fullfile(cuadmm_outdir(),'baseline','dumps');
TSV  = fullfile(cuadmm_outdir(),'baseline','bl_big.tsv');
if ~exist(TSV,'file')
    fid=fopen(TSV,'w');
    fprintf(fid,'id\tn\tstatus\tm\tKf\tKs\tnnzAt\tnvar\tvec_len\tt_build\tt_dump\tschur_GB\tnote\n');
    fclose(fid);
end
for n = [24 32]
    id = sprintf('scale_stab_n%02d',n);
    st_s=''; note=''; m=NaN; Kf=NaN; Ks=''; nz=NaN; nv=NaN; vl=NaN; tb=NaN; td=NaN; sg=NaN;
    try
        clear stateNameGenerator
        pvar s t
        x = pde_var(n,s,[0,1]);
        PIE = initialize(convert([diff(x,t,1)==diff(x,s,2)+0.5*pi^2*x; ...
                                  subs(x,s,0)==0; subs(x,s,1)==0]));
        stg = lpisettings('heavy');
        stg.sos_opts.solver='mosek'; stg.sos_opts.simplify=false;
        t0=tic; evalc('out = stab_mirror(PIE,stg);'); tb=toc(t0);
        S = sdpshape(out.prog);
        m=S.m; Kf=S.Kf; Ks=mat2str(S.Ks); nz=S.nnzAt;
        nv = S.Kf+sum(double(S.Ks).^2);
        vl = S.Kf+sum(double(S.Ks).*(double(S.Ks)+1)/2);
        sg = 8*double(m)^2/2^30;
        Atf=[]; bf=[];
        for q=1:out.prog.expr.num, Atf=[Atf,out.prog.expr.At{q}]; bf=[bf;out.prog.expr.b{q}]; end
        RR = mkRR(out.prog);  bscl = norm(full(bf)); if bscl==0, bscl=1; end
        D.At = RR'*Atf;  D.b = bf/bscl;  D.c = sparse(size(Atf,1),1);
        D.K = struct('f',S.Kf,'l',0,'q',[],'s',S.Ks);  D.Ns=S.Ks; D.Kf=S.Kf;
        t1=tic;
        save(fullfile(DMP,[id '.mat']),'-struct','D','-v7.3');
        save(fullfile(DMP,[id '_meta.mat']),'RR','bscl','S');
        dump2cuadmm(fullfile(DMP,[id '.mat']),fullfile(DMP,id));
        td=toc(t1);
        st_s='ok';
    catch ME
        st_s='ERR'; note=regexprep(ME.message,'[\t\r\n]+',' ');
    end
    fid=fopen(TSV,'a');
    fprintf(fid,'%s\t%d\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%.1f\t%s\n', ...
        id,n,st_s,num2str(m),num2str(Kf),Ks,num2str(nz),num2str(nv),num2str(vl), ...
        num2str(tb),num2str(td),sg,note);
    fclose(fid);
    fprintf('BIG %s %s n=%d m=%s nvar=%s t_build=%s t_dump=%s schur=%.0fGB %s\n', ...
        id,st_s,n,num2str(m),num2str(nv),num2str(tb),num2str(td),sg,note);
end
fprintf('BIGDONE\n');

function RR = mkRR(prog)
RR = speye(prog.var.idx{1}-1);
for i=1:prog.var.num
    sz = prog.var.idx{i+1}-prog.var.idx{i};
    switch prog.var.type{i}
        case 'poly', RR = spantiblkdiag(RR,speye(sz));
        case 'sos',  RR = spblkdiag(RR,speye(sz));
        otherwise,   error('bl_big:vartype','unexpected var type %s',prog.var.type{i});
    end
end
for i=1:prog.extravar.num
    RR = spblkdiag(RR,speye(prog.extravar.idx{i+1}-prog.extravar.idx{i}));
end
end
