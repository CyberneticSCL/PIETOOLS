% bis_dump.m -- dump FIXED-gamma H-infinity feasibility programs for cuADMM.
% With gamma numeric there is no lpidecvar, no lpi_ineq and no objective, so
% c = 0 and a pre-solve capture is valid (the addextrasosvar hazard only
% applies to 'ineq' expressions, and there are none).  Self-checked with
% SeDuMi on the dumped data.
cuadmm_path;
HERE = fileparts(mfilename('fullpath'));
ROOT = fullfile(cuadmm_outdir(),'cu_bis'); if ~exist(ROOT,'dir'), mkdir(ROOT); end
CASES = { 'light',0.0600 ; 'light',0.0523 ; 'heavy',0.0600 ; 'heavy',0.0523 };
fprintf('BD lab|m|Kf|Ns|nnz|vec|c_nnz|sedumi_feas|rel_b\n');
for k = 1:size(CASES,1)
    pre = CASES{k,1}; g = CASES{k,2};
    [prog,S] = mkfixed(pre,g);
    Atf=[];bf=[];
    for i=1:prog.expr.num, Atf=[Atf,prog.expr.At{i}]; bf=[bf;prog.expr.b{i}]; end
    RR = mkRR(prog);
    bscl = norm(full(bf)); if bscl==0||~isfinite(bscl), bscl=1; end
    D.At = RR'*Atf;  D.b = bf/bscl;  D.c = sparse(size(Atf,1),1);
    D.K = struct('f',S.Kf,'l',0,'q',[],'s',S.Ks);  D.Ns=S.Ks; D.Kf=S.Kf;
    lab = sprintf('bis_%s_g%04.0f',pre,g*1e5);
    save(fullfile(ROOT,[lab '.mat']),'-struct','D');
    info = dump2cuadmm(fullfile(ROOT,[lab '.mat']),fullfile(ROOT,lab));
    save(fullfile(ROOT,[lab '_meta.mat']),'RR','bscl','S','g','pre');
    pars.fid=0; [xs,~,i2]=sedumi(D.At,full(D.b),D.c,D.K,pars);
    rb = norm(full(D.At'*xs-D.b))/norm(full(D.b));
    ok = (i2.numerr==0)&&(rb<1e-6)&&(abs(rb-1)>1e-6);
    fprintf('BD %s|%d|%d|[%s]|%d|%d|%d|%d|%.3e\n',lab,S.m,S.Kf, ...
            strtrim(num2str(S.Ks)),info.nnz_At,info.vec_len,nnz(D.c),ok,rb);
end
fprintf('BDDONE\n');

function [prog,S] = mkfixed(pre,gam)
st = lpisettings(pre);
pvar s t
x = pde_var('state',1,s,[0,1]);
w = pde_var('input',1); z = pde_var('output',1);
PIE = initialize(convert([diff(x,t,1)==diff(x,s,2)+2*x+s*w;
      z == int(x,s,[0,1]); subs(x,s,0)==0; subs(x,s,1)==0]));
Top=PIE.T; Aop=PIE.A; Bwop=PIE.Bw; Czop=PIE.Cz; Dzwop=PIE.Dzw;
prog = lpiprogram(PIE.vars(:,1),PIE.vars(:,2),PIE.dom);
[prog,R1op] = poslpivar(prog,Top.dim,st.dd1,st.options1);
if st.override1~=1
    [prog,P2op] = poslpivar(prog,Top.dim,st.dd12,st.options12);
    Rop = R1op+P2op;
else, Rop = R1op; end
Qdeg = get_lpivar_degs(Rop,Top);
[prog,Qop] = lpivar(prog,Top.dim,Qdeg);
prog = lpi_eq(prog,Top'*Qop-Rop);
Iw = mat2opvar(eye(size(Bwop,2)),Bwop.dim(:,2),PIE.vars,PIE.dom);
Iz = mat2opvar(eye(size(Czop,1)),Czop.dim(:,1),PIE.vars,PIE.dom);
Dop = [-gam*Iw,    Dzwop',   Bwop'*Qop;
        Dzwop,     -gam*Iz,  Czop;
        Qop'*Bwop, Czop',    Aop'*Qop+Qop'*Aop];
[prog,De1op] = poslpivar(prog,Dop.dim,st.dd2,st.options2);
if st.override2~=1
    [prog,De2op] = poslpivar(prog,Dop.dim,st.dd3,st.options3);
    Deop = De1op+De2op;
else, Deop = De1op; end
prog = lpi_eq(prog,Deop+Dop,'symmetric');
S = cuadmm_private('sdpshape',prog);
end

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
