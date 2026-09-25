% t13_mosek.m -- is Mosek usable through PIETOOLS, and how much faster?
% Worth a minute: the settings search needs dozens of 2-D solves.
cuadmm_path;
pvar s1 s2 t
x = pde_var('state',1,[s1;s2],[0,1;0,1]);
PIE = convert([diff(x,t,1)==diff(x,s1,2)+diff(x,s2,2)+2*x;
               subs(x,s1,0)==0; subs(x,s1,1)==0;
               subs(x,s2,0)==0; subs(x,s2,1)==0]);
st = lpisettings('heavy');
st.settings_2d.eppos = 1e-2*[1;1;1;1]; st.eppos=1e-2; st.eppos2=1e-2;
for slv = {'sedumi','mosek'}
    try
        st2 = st; st2.sos_opts.solver = slv{1};
        t0=tic; evalc('prog = PIETOOLS_stability_2D(PIE,st2);'); tw=toc(t0);
        xv=prog.solinfo.RRx(:); Atf=[];bf=[];
        for i=1:prog.expr.num, Atf=[Atf,prog.expr.At{i}]; bf=[bf;prog.expr.b{i}]; end
        r=norm(full(Atf'*xv-bf))/norm(full(bf));
        fprintf('T13 %-7s t=%.1fs rel_b=%.3e\n', slv{1}, tw, r);
    catch ME
        fprintf('T13 %-7s FAILED %s\n', slv{1}, strrep(ME.message,newline,' '));
    end
end
fprintf('T13DONE\n');
