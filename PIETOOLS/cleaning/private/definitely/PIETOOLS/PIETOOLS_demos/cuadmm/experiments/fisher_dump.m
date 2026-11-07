% fisher_dump.m -- dump the NONLINEAR PIESOS Fisher SDP for cuADMM.
%
% Safe to capture PRE-SOLVE here: every expression in this program is type
% 'eq' (verified), so the addextrasosvar hazard that invalidates a pre-solve
% capture of the Hinf executives does not apply.
% Three modes so the objective's effect is separable from the nonlinearity:
%   nobnd  c = 0, V upper bound dropped        (weakest, pure feasibility)
%   gamfix c = 0, bound kept with gamma numeric (faithful to Thm 1)
%   opt    c != 0, gamma minimised              (the shipped default)
cuadmm_path;
HERE = fileparts(mfilename('fullpath'));
ROOT = fullfile(cuadmm_outdir(),'cu_fish'); if ~exist(ROOT,'dir'), mkdir(ROOT); end
R = 4.0;
MODES = {'nobnd',false; 'gamfix',2.0; 'opt',true};
fprintf('FD mode|m|Kf|nblk|Ns|nnz|vec|c_nnz|sedumi_ok|rel_b|t_sedumi\n');
for k = 1:size(MODES,1)
    nm = MODES{k,1};
    [prog,I0] = cuadmm_private('fisher_prog',R,MODES{k,2},false);       % build only
    Atf=[]; bf=[];
    for i=1:prog.expr.num, Atf=[Atf,prog.expr.At{i}]; bf=[bf;prog.expr.b{i}]; end
    RR = mkRR(prog);
    S = cuadmm_private('sdpshape',prog);
    nvar = size(Atf,1);
    c = sparse(nvar,1);
    if I0.c_nnz>0, c(1:prog.var.idx{end}-1) = prog.objective; end
    bscl = norm(full(bf)); if bscl==0||~isfinite(bscl), bscl=1; end
    D.At = RR'*Atf;  D.b = bf/bscl;  D.c = RR'*c;
    D.K = struct('f',S.Kf,'l',0,'q',[],'s',S.Ks);  D.Ns=S.Ks; D.Kf=S.Kf;
    lab = ['fish_' nm];
    save(fullfile(ROOT,[lab '.mat']),'-struct','D','-v7.3');
    info = dump2cuadmm(fullfile(ROOT,[lab '.mat']),fullfile(ROOT,lab));
    save(fullfile(ROOT,[lab '_meta.mat']),'RR','bscl','S','R','nm');
    t0=tic; pars.fid=0; [xs,~,i2]=sedumi(D.At,full(D.b),D.c,D.K,pars); ts=toc(t0);
    rb = norm(full(D.At'*xs-D.b))/norm(full(D.b));
    ok = (i2.numerr==0)&&(rb<1e-5)&&(abs(rb-1)>1e-6);
    fprintf('FD %s|%d|%d|%d|[%s]|%d|%d|%d|%d|%.3e|%.1f\n', nm,S.m,S.Kf,S.nblk, ...
            strtrim(num2str(S.Ks)),info.nnz_At,info.vec_len,nnz(D.c),ok,rb,ts);
end
fprintf('FDDONE\n');

function RR = mkRR(prog)
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
end
