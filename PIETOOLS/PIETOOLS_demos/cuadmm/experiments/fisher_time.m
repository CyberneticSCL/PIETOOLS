% fisher_time.m -- SOLVE-ONLY times, so the cuADMM comparison is like-for-like.
% The ladder's t column wrapped PDE conversion + LPI build + solve; cuADMM's
% CUADMM_TIMING solve is the solve alone.  Build once here, time lpisolve only.
% Also settles whether sol_opts.simplify (hardcoded true in PIESOS_Fisher.m:218)
% changes the system Mosek actually sees -- if it does, the shipped path does
% not solve the SDP that was dumped.
cuadmm_path;
R = 4.0;
MODES = {'nobnd',false; 'gamfix',2.0; 'opt',true};
fprintf('FT mode|simplify|m_pre|m_post|Ks_post|t_solve1|t_solve2|rel_b|feasratio|numerr\n');
for k = 1:size(MODES,1)
    for sp = [false true]
        prog = cuadmm_private('fisher_prog',R,MODES{k,2},false);
        S0 = cuadmm_private('sdpshape',prog);
        so = struct('solver','mosek','simplify',sp);
        t1=tic; evalc('sol = lpisolve(prog,so);'); tA=toc(t1);
        prog2 = cuadmm_private('fisher_prog',R,MODES{k,2},false);
        t2=tic; evalc('sol2 = lpisolve(prog2,so);'); tB=toc(t2);
        Atf=[];bf=[];
        for i=1:sol.expr.num, Atf=[Atf,sol.expr.At{i}]; bf=[bf;sol.expr.b{i}]; end
        xv = sol.solinfo.RRx(:);
        rb = norm(full(Atf'*xv-bf))/norm(full(bf));
        S1 = cuadmm_private('sdpshape',sol);  I = sol.solinfo.info;
        fprintf('FT %s|%d|%d|%d|[%s]|%.1f|%.1f|%.3e|%+.4f|%d\n', MODES{k,1},sp, ...
            S0.m,S1.m,strtrim(num2str(S1.Ks)),tA,tB,rb,I.feasratio,I.numerr);
    end
end
fprintf('FTDONE\n');
