% run_fisher.m -- the shipped PIESOS Fisher example, inspected after the fact.
% NOTE the script begins with `clear`, so nothing declared before run() survives;
% everything reported here is read from the workspace the script leaves behind.
cuadmm_path;
run(which('PIESOS_Fisher'));   % this checkout's copy, not a hardcoded workstation path
if exist('prog_sol','var')
    S = cuadmm_private('sdpshape',prog_sol);  I = prog_sol.solinfo.info;
    fprintf('RF R=%g  m=%d nvar=%d Kf=%d nblk=%d Nmax=%d Ns=[%s]\n', ...
            R,S.m,S.nvar,S.Kf,S.nblk,S.Nmax,strtrim(num2str(S.Ks)));
    fprintf('RF nnzAt=%d  c_nnz=%d  pinf=%d dinf=%d numerr=%d feasratio=%+.4f\n', ...
            nnz([prog_sol.expr.At{:}]),nnz(prog_sol.objective), ...
            I.pinf,I.dinf,I.numerr,I.feasratio);
    Atf=[];bf=[];
    for i=1:prog_sol.expr.num, Atf=[Atf,prog_sol.expr.At{i}]; bf=[bf;prog_sol.expr.b{i}]; end
    xv = prog_sol.solinfo.RRx(:);
    fprintf('RF |b|=%.4e  rel_b=%.3e\n', norm(full(bf)), ...
            norm(full(Atf'*xv-bf))/norm(full(bf)));
    fprintf('RF expr types: %s\n', strjoin(unique(prog_sol.expr.type),','));
    if exist('gam_sol','var'), fprintf('RF gam_sol=%.6f\n', double(gam_sol)); end
else
    fprintf('RF prog_sol MISSING\n');
end
fprintf('RFDONE\n');
