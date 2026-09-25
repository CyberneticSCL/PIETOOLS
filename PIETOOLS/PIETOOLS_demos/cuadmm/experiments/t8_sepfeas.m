% t8_sepfeas.m -- is sep=1 INFEASIBLE, or does the solver merely fail?
% SeDuMi's info distinguishes them: pinf=1 is a proof of primal infeasibility,
% numerr>0 is a numerical failure, and feasratio near -1 also indicates
% infeasibility.  Swept over lam so a feasibility MARGIN effect is separable
% from a structural one (lam* = pi^2 = 9.8696 for this plant).
cuadmm_path;
pvar s t
fprintf('T8 lam*=pi^2=%.4f\n',pi^2);
for lam = [0 2 5]
  x   = pde_var('state',1,s,[0,1]);
  PIE = convert([diff(x,t,1)==diff(x,s,2)+lam*x; subs(x,s,0)==0; subs(x,s,1)==0]);
  for sepv = [0 1]
    st = lpisettings('light');
    st.options1.sep=sepv; st.options12.sep=sepv;
    st.options2.sep=sepv; st.options3.sep=sepv;
    try
      M = cuadmm_private('stab_mirror',PIE,st);
      S = cuadmm_private('sdpshape',M.prog);
      opts.solver='sedumi'; opts.params.fid=0;
      evalc('sol = lpisolve(M.prog,opts);');
      I = sol.solinfo.info;
      has = isfield(sol.solinfo,'RRx') && ~isempty(sol.solinfo.RRx);
      fprintf(['T8 lam=%-4g sep=%d m=%-4d eig=%-6g | pinf=%d dinf=%d numerr=%d ' ...
               'feasratio=%+.4f | solution=%d\n'], ...
              lam,sepv,S.m,S.eigcost,I.pinf,I.dinf,I.numerr,I.feasratio,has);
    catch ME
      fprintf('T8 lam=%-4g sep=%d ERROR %s\n',lam,sepv,strrep(ME.message,newline,' '));
    end
  end
end
fprintf('T8DONE\n');
