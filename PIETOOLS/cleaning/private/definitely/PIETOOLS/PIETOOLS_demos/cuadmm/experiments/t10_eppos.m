% t10_eppos.m -- MAINTAINER'S CORRECTION: ||b|| ~ 1e-6 is not a property of the
% stability class, it is eppos2, set too low.  Small eppos2 makes every
% solution small and confounds the ABSOLUTE metrics.  Sweep it and report both
% absolute and relative measures so the invariance is visible.
cuadmm_path;
pvar s t
x   = pde_var('state',1,s,[0,1]);
PIE = convert([diff(x,t,1)==diff(x,s,2)+2*x; subs(x,s,0)==0; subs(x,s,1)==0]);

fprintf('CEP eppos2|normb|normx|normDop|rel_b|rel_cu|op_ind|psd_min|psd_relmin|Pmin|blocks_mineig\n');
for e2 = [1e-6 1e-4 1e-2 1]
    st = lpisettings('light');  st.eppos2 = e2;
    try
        M = cuadmm_private('stab_mirror',PIE,st);
        opts.solver='sedumi'; opts.params.fid=0;
        evalc('sol = lpisolve(M.prog,opts);');
        Atf=[];bf=[];
        for i=1:sol.expr.num, Atf=[Atf,sol.expr.At{i}]; bf=[bf;sol.expr.b{i}]; end
        r = cuadmm_private('lpi_resid',sol,M,Atf,bf,16);
        bm = sprintf('%+.1e ', r.psd_blocks(:,1));
        fprintf('CEP %.0e|%.3e|%.3e|%.3e|%.3e|%.3e|%.3e|%+.3e|%+.3e|%+.3e|%s\n', ...
            e2,r.normb,r.normx,r.normDop,r.rel_b,r.rel_cu,r.op_ind, ...
            r.psd_min,r.psd_relmin,r.pop_mineig,strtrim(bm));
    catch ME
        fprintf('CEP %.0e|FAILED|%s\n',e2,strrep(ME.message,newline,' '));
    end
end
fprintf('CEPDONE\n');
