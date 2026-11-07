% dump1d.m -- dump the 1-D stability LPI at Dup=1 and Dup=2 so the predicted
% cuADMM cost ratio (1.15x) can be MEASURED rather than modelled.
cuadmm_path;
HERE = fileparts(mfilename('fullpath'));
cuadmm_shadow('lpisolve');
global CENSUS_PROG
ROOT = fullfile(cuadmm_outdir(),'cu_1d'); if ~exist(ROOT,'dir'), mkdir(ROOT); end
fprintf('D1 label|m|nvar|Kf|nblk|Ns|eigcost|nnzAt|vec_len\n');
for Dup = [1 2]
    lab = sprintf('1d_Dup%d',Dup);
    CENSUS_PROG = [];
    pvar s t
    x = pde_var('state',1,s,[0,1]);
    PIE = convert([diff(x,t,1)==diff(x,s,2)+2*x; subs(x,s,0)==0; subs(x,s,1)==0]);
    st = lpisettings('light'); st.eppos2=1e-2; st.eppos=1e-2;
    n1=1;n2=1;n3=1;n4=2;
    st.dd2 = {n1+Dup,   [n2+Dup-1, n3+Dup,   n4+Dup  ], [n2+Dup-1, n3+Dup,   n4+Dup  ]};
    st.dd3 = {n1+Dup-1, [n2+Dup-2, n3+Dup-1, n4+Dup-1], [n2+Dup-2, n3+Dup-1, n4+Dup-1]};
    evalc('PIETOOLS_PIE2PDEstability(PIE,st);');
    prog = CENSUS_PROG;  S = cuadmm_private('sdpshape',prog);
    Atf=[];bf=[];
    for i=1:prog.expr.num, Atf=[Atf,prog.expr.At{i}]; bf=[bf;prog.expr.b{i}]; end
    RR = speye(prog.var.idx{1}-1);
    for i=1:prog.var.num
        sz = prog.var.idx{i+1}-prog.var.idx{i};
        switch prog.var.type{i}
            case 'poly', RR = spantiblkdiag(RR,speye(sz));
            case 'sos',  RR = spblkdiag(RR,speye(sz));
        end
    end
    for i=1:prog.extravar.num
        RR = spblkdiag(RR,speye(prog.extravar.idx{i+1}-prog.extravar.idx{i}));
    end
    D.At = RR'*Atf;  D.b = bf/max(norm(full(bf)),eps);
    D.c = sparse(size(Atf,1),1);
    D.K = struct('f',S.Kf,'l',0,'q',[],'s',S.Ks);  D.Ns=S.Ks; D.Kf=S.Kf;
    dfile = fullfile(ROOT,[lab '.mat']); save(dfile,'-struct','D');
    info = dump2cuadmm(dfile,fullfile(ROOT,lab));
    fprintf('D1 %s|%d|%d|%d|%d|[%s]|%g|%d|%d\n',lab,S.m,S.nvar,S.Kf,S.nblk, ...
            strtrim(num2str(S.Ks)),S.eigcost,info.nnz_At,info.vec_len);
end
fprintf('D1DONE\n');
cuadmm_shadow('off');   % remove the shadow: it stays active for the whole session otherwise
