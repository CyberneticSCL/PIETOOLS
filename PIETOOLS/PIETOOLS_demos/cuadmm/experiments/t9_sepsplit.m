% t9_sepsplit.m -- sep on the LYAPUNOV operator and sep on the NEGATIVITY
% operator are different restrictions.  t8 set both at once and got a clean
% infeasibility certificate; this separates them, so any part of the saving
% that IS available can be identified.
cuadmm_path;
pvar s t
x   = pde_var('state',1,s,[0,1]);
PIE = convert([diff(x,t,1)==diff(x,s,2)+2*x; subs(x,s,0)==0; subs(x,s,1)==0]);

combos = [0 0; 1 0; 0 1; 1 1];
for k = 1:size(combos,1)
    sLF = combos(k,1);  sD = combos(k,2);
    st = lpisettings('light');
    st.options1.sep=sLF; st.options12.sep=sLF;
    st.options2.sep=sD;  st.options3.sep=sD;
    try
        M = cuadmm_private('stab_mirror',PIE,st);  S = cuadmm_private('sdpshape',M.prog);
        opts.solver='sedumi'; opts.params.fid=0;
        evalc('sol = lpisolve(M.prog,opts);');
        I = sol.solinfo.info;
        ok = isfield(sol.solinfo,'RRx') && ~isempty(sol.solinfo.RRx);
        extra = '';
        if ok
            Atf=[];bf=[];
            for i=1:sol.expr.num, Atf=[Atf,sol.expr.At{i}]; bf=[bf;sol.expr.b{i}]; end
            r = cuadmm_private('lpi_resid',sol,M,Atf,bf,16);
            extra = sprintf(' | rel_b=%.2e op_ind=%.2e op_indP=%.2e op_eqQ=%.2e psd_min=%+.1e Pmin=%+.1e normDop=%.2e TRIV=%d', ...
                            r.rel_b,r.op_ind,r.op_indP,r.op_eqQ,r.psd_min,r.pop_mineig,r.normDop,r.trivial);
        end
        fprintf('T9 sepLF=%d sepD=%d m=%-4d eig=%-6g Ns=[%s] | pinf=%d feasratio=%+.3f sol=%d%s\n',...
                sLF,sD,S.m,S.eigcost,strtrim(num2str(S.Ks)),I.pinf,I.feasratio,ok,extra);
    catch ME
        fprintf('T9 sepLF=%d sepD=%d ERROR %s\n',sLF,sD,strrep(ME.message,newline,' '));
    end
end
fprintf('T9DONE\n');
