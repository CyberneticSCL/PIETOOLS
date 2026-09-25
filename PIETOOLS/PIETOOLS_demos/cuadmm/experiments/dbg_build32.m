% dbg_build32.m -- where does the LPI BUILD die at n=32?
% The scaling ladder showed the wall is not in the SDP solver: n=32 never
% reaches one. stab_mirror requests a 9613344x1536 DENSE array (123.8 GB) during
% construction. CLAUDE.md item 2 forbids exactly this -- a dense allocation on a
% dimension involving the decision-variable count. Locate the call site.
cuadmm_path;
for n = [24 32]
    fprintf('B32 n=%d start\n',n);
    try
        clear stateNameGenerator
        pvar s t
        x = pde_var(n,s,[0,1]);
        PIE = initialize(convert([diff(x,t,1)==diff(x,s,2)+0.5*pi^2*x; ...
                                  subs(x,s,0)==0; subs(x,s,1)==0]));
        st = lpisettings('heavy');
        st.sos_opts.solver='mosek'; st.sos_opts.simplify=false;
        t0=tic; out = cuadmm_private('stab_mirror',PIE,st); tb=toc(t0);
        S = cuadmm_private('sdpshape',out.prog);
        fprintf('B32 n=%d OK m=%d nvar=%d t_build=%.1f\n',n,S.m,S.Kf+sum(double(S.Ks).^2),tb);
    catch ME
        fprintf('B32 n=%d ERR %s\n',n,regexprep(ME.message,'[\r\n\t]+',' '));
        for k=1:numel(ME.stack)
            fprintf('B32   at %s line %d\n',ME.stack(k).name,ME.stack(k).line);
        end
    end
end
fprintf('B32DONE\n');
