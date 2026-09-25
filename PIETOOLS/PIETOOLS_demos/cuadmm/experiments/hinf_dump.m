% hinf_dump.m -- dump the H-infinity gain SDP for cuADMM, WITH A NONZERO
% OBJECTIVE.
%
% CAPTURED PRE-SOLVE IS WRONG HERE.  sossolve.m:177-183 calls addextrasosvar
% INSIDE the solve to create the slack cone variable for every 'ineq'
% expression, and only then assembles Atf/bf.  PIETOOLS_Hinf_gain imposes
% gam >= 0 with lpi_ineq, so a pre-solve capture omits that block entirely:
% measured Ns = [10 17 8] against the solved program's [10 17 8 1], and
% SeDuMi on that dump returns gamma = -0.0 with numerr=1, because minimising
% gamma without gam >= 0 is unbounded.  So extract from the SOLVED program,
% which retains the modification.
%
% Two further things this program needs that the feasibility ones did not:
%   c   = RR'*c with sos.objective on the leading block
%   bscl: b -> b/bscl scales the argmin by 1/bscl, so gamma scales BACK by
%         bscl on import (sossolve.m:678 does solinfo.x = x*bscl).
% SELF-CHECK: re-solve the dumped data with SeDuMi directly and compare gamma
% against the executive's.
cuadmm_path;
HERE = fileparts(mfilename('fullpath'));
ROOT = fullfile(cuadmm_outdir(),'cu_hinf'); if ~exist(ROOT,'dir'), mkdir(ROOT); end

pvar s t
x = pde_var('state',1,s,[0,1]);
w = pde_var('input',1);
z = pde_var('output',1);
PIE = convert([diff(x,t,1)==diff(x,s,2)+2*x+s*w;
               z == int(x,s,[0,1]);
               subs(x,s,0)==0; subs(x,s,1)==0]);

for pre = {'light','heavy'}
    st = lpisettings(pre{1});  st.sos_opts.solver = 'mosek';
    evalc('[prog,~,gam] = PIETOOLS_Hinf_gain(PIE,st);');
    g_ref = double(gam);
    S = cuadmm_private('sdpshape',prog);

    Atf=[]; bf=[];
    for i=1:prog.expr.num, Atf=[Atf,prog.expr.At{i}]; bf=[bf;prog.expr.b{i}]; end
    RR = mkRR(prog);
    nvar = size(Atf,1);
    c = sparse(nvar,1);
    c(1:prog.var.idx{end}-1) = prog.objective;
    obj_idx = find(prog.objective ~= 0);  obj_idx = obj_idx(1);
    bscl = norm(full(bf)); if bscl==0 || ~isfinite(bscl), bscl = 1; end

    D.At = RR'*Atf;  D.b = bf/bscl;  D.c = RR'*c;
    D.K  = struct('f',S.Kf,'l',0,'q',[],'s',S.Ks);
    D.Ns = S.Ks;  D.Kf = S.Kf;
    lab = ['hinf_' pre{1}];
    save(fullfile(ROOT,[lab '.mat']),'-struct','D');
    info = dump2cuadmm(fullfile(ROOT,[lab '.mat']), fullfile(ROOT,lab));
    save(fullfile(ROOT,[lab '_meta.mat']),'RR','bscl','obj_idx','S','g_ref');

    pars.fid = 0;
    [xs,~,i2] = sedumi(D.At,full(D.b),D.c,D.K,pars);
    RRx = (RR*xs)*bscl;
    g_dump = RRx(obj_idx);
    fprintf(['HD %-5s m=%d Kf=%d Ns=[%s] nnz=%d vec=%d bscl=%.4f | gam_ref=%.8f ' ...
             'gam_dump=%.8f  reldiff=%.2e  numerr=%d\n'], ...
        pre{1},S.m,S.Kf,strtrim(num2str(S.Ks)),info.nnz_At,info.vec_len,bscl, ...
        g_ref,g_dump,abs(g_ref-g_dump)/max(abs(g_ref),eps),i2.numerr);
end
fprintf('HDDONE\n');

function RR = mkRR(prog)
RR = speye(prog.var.idx{1}-1);
for i = 1:prog.var.num
    sz = prog.var.idx{i+1}-prog.var.idx{i};
    switch prog.var.type{i}
        case 'poly', RR = spantiblkdiag(RR,speye(sz));
        case 'sos',  RR = spblkdiag(RR,speye(sz));
    end
end
for i = 1:prog.extravar.num
    RR = spblkdiag(RR,speye(prog.extravar.idx{i+1}-prog.extravar.idx{i}));
end
end
