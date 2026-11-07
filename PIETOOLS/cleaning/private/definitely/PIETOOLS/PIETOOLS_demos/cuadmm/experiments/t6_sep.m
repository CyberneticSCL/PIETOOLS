% t6_sep.m -- THE sep=1 QUESTION.  The shape map says options.sep=1 is cheaper
% on every axis for both solver families and is off by default.  The thing a
% shape map cannot tell us is whether it still CERTIFIES: sep=1 restricts the
% operator class to R1=R2, so it may simply lose feasibility.
%
% Solve with SeDuMi (reference solver) at sep=0 and sep=1, on the same plant,
% and report the full residual panel for each.
cuadmm_path;
pvar s t
x   = pde_var('state',1,s,[0,1]);
lam = 2;                                    % well inside the pi^2 limit
PIE = convert([diff(x,t,1)==diff(x,s,2)+lam*x; subs(x,s,0)==0; subs(x,s,1)==0]);

for pre = {'light','heavy'}
for sepv = [0 1]
    st = lpisettings(pre{1});
    st.options1.sep = sepv;  st.options12.sep = sepv;
    if isfield(st,'options2'), st.options2.sep = sepv; end
    if isfield(st,'options3'), st.options3.sep = sepv; end
    try
        M  = cuadmm_private('stab_mirror',PIE,st);
        S  = cuadmm_private('sdpshape',M.prog);
        t0 = tic;
        opts.solver = 'sedumi';  opts.params.fid = 0;
        evalc('sol = lpisolve(M.prog,opts);');
        tw = toc(t0);
        Atf=[]; bf=[];
        for i=1:sol.expr.num, Atf=[Atf, sol.expr.At{i}]; bf=[bf; sol.expr.b{i}]; end
        r = cuadmm_private('lpi_resid',sol,M,Atf,bf,16);
        fprintf(['T6 %-6s sep=%d m=%-4d eig=%-7g | rel_b=%.3e rel_cu=%.3e | ' ...
                 'psd_min=%+.3e | op_ind=%.3e op_hs=%.3e op_coef=%.3e | ' ...
                 'Pmin=%+.3e | t=%.1fs\n'], ...
                pre{1},sepv,S.m,S.eigcost,r.rel_b,r.rel_cu,r.psd_min, ...
                r.op_ind,r.op_hs,r.op_coef,r.pop_mineig,tw);
    catch ME
        fprintf('T6 %-6s sep=%d FAILED %s\n',pre{1},sepv,strrep(ME.message,newline,' '));
    end
end
end
fprintf('T6DONE\n');
