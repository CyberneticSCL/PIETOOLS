% t12_2dfeas.m -- is 2-D stability certifiable AT ALL here?  Discriminating
% test, predicting which cases pass AND which fail:
%   lam=0 (pure 2-D heat, V=int x^2 works trivially) must pass if the
%   machinery is sound; if even that fails the problem is the settings or the
%   solver, not the stability margin.
% Crossed with the settings preset, because 'light' may simply be too weak in
% 2-D (the 2-D degree structure is far richer than 1-D).
cuadmm_path;
pvar s1 s2 t
fprintf('F2 lam|preset|m|Nmax|pinf|numerr|feasratio|rel_b|trivial|t\n');
for lam = [0 2]
for pre = {'light','heavy'}
    x = pde_var('state',1,[s1;s2],[0,1;0,1]);
    PIE = convert([diff(x,t,1)==diff(x,s1,2)+diff(x,s2,2)+lam*x;
                   subs(x,s1,0)==0; subs(x,s1,1)==0;
                   subs(x,s2,0)==0; subs(x,s2,1)==0]);
    st = lpisettings(pre{1});
    st.settings_2d.eppos = 1e-2*[1;1;1;1]; st.eppos=1e-2; st.eppos2=1e-2;
    try
        t0=tic; evalc('prog = PIETOOLS_stability_2D(PIE,st);'); tw=toc(t0);
        S=cuadmm_private('sdpshape',prog); I=prog.solinfo.info; xv=prog.solinfo.RRx(:);
        Atf=[];bf=[];
        for i=1:prog.expr.num, Atf=[Atf,prog.expr.At{i}]; bf=[bf;prog.expr.b{i}]; end
        res=norm(full(Atf'*xv-bf)); nb=norm(full(bf));
        triv = (abs(res/nb-1)<=1e-6) || norm(xv)<=1e-12;
        fprintf('F2 %g|%s|%d|%d|%d|%d|%+.4f|%.3e|%d|%.1f\n', ...
                lam,pre{1},S.m,S.Nmax,I.pinf,I.numerr,I.feasratio,res/nb,triv,tw);
    catch ME
        fprintf('F2 %g|%s|ERROR|%s\n',lam,pre{1},strrep(ME.message,newline,' '));
    end
end
end
fprintf('F2DONE\n');
