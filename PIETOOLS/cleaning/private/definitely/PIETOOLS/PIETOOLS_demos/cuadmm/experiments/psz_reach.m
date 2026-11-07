% psz_reach.m -- is "reach 0.00 for the shipped settings, searched P" a real
% property or a Mosek artefact?
%
% poslpivar_2d.m line 158 (my patch header) asserts the psatz result "carries
% over to a searched Lyapunov operator (reach 0.99 at m = 3912, against 0.00 for
% the shipped settings)". Two things now contradict the 0.00:
%   - my own note records "shipped heavy with a FULL P search reaches only
%     0.25*lam* at m=8950" -- 0.25, not 0.00;
%   - another session's SeDuMi reference certifies this exact program at
%     0.10*lam*, psatz off, searched P, rel 9.46e-07;
% and today's Mosek run failed at 0.10, 0.50 AND 0.90 with psatz off, which is
% precisely the "Mosek is the weak solver here" pattern.
%
% If 0.00 was measured under Mosek it is an artefact, and the header line would
% tell a reader the psatz is MANDATORY for a searched P when it is not. Measured
% here with BOTH solvers on the same programs so the claim can be rewritten to
% say what was actually observed, and with which solver.
%
% Acceptance is deliberately not a single number: rel_b, the trivial-point test
% (rel_b == 1 with ||x|| ~ 0) and the PSD margin from the CONE vector are all
% printed, because numerr has already passed a program with rel_b = 1.000 in
% this campaign.
cuadmm_path;
LAMSTAR = 2*pi^2;
FR = [0.10 0.25 0.50];
fprintf('PR solver|frac|m|rel_b|psd_relmin|normx|trivial|numerr|t_s\n');
for slv = {'sedumi','mosek'}
    for f = FR
        try
            clear stateNameGenerator
            pvar s1 s2 t
            x = pde_var(1,[s1;s2],[0,1;0,1]);
            PIE = initialize(convert([diff(x,t,1)==diff(x,s1,2)+diff(x,s2,2)+f*LAMSTAR*x;
                   subs(x,s1,0)==0; subs(x,s1,1)==0; subs(x,s2,0)==0; subs(x,s2,1)==0]));
            st = lpisettings('light');
            st.sos_opts.solver = slv{1};  st.sos_opts.simplify = false;
            st.settings_2d.eq_use_psatz = [0;0];     % psatz OFF = the contested config
            t0 = tic;
            evalc('[sol,Pop,Qop] = PIETOOLS_stability_2D(PIE,st);');
            tw = toc(t0);
            S = cuadmm_private('sdpshape',sol);
            Atf=[];bf=[];
            for q=1:sol.expr.num, Atf=[Atf,sol.expr.At{q}]; bf=[bf;sol.expr.b{q}]; end
            xr = sol.solinfo.RRx(:);  xc = sol.solinfo.x(:);
            rb = norm(full(Atf'*xr-bf))/max(norm(full(bf)),eps);
            off=S.Kf; pmin=inf; pmax=-inf;
            for k=1:numel(S.Ks)
                N=double(S.Ks(k)); Xk=reshape(xc(off+(1:N^2)),N,N);
                ev=eig((Xk+Xk')/2); pmin=min(pmin,min(ev)); pmax=max(pmax,max(ev)); off=off+N^2;
            end
            triv = double((abs(rb-1)<=1e-6) || norm(xr)<=1e-12);
            I = sol.solinfo.info;
            fprintf('PR %s|%.2f|%d|%.4e|%+.3e|%.4e|%d|%d|%.1f\n', ...
                slv{1},f,S.m,rb,pmin/max(pmax,eps),norm(xr),triv,gf(I,'numerr'),tw);
        catch ME
            fprintf('PR %s|%.2f|ERR|%s\n',slv{1},f,regexprep(ME.message,'[\t\r\n]+',' '));
        end
    end
end
fprintf('PRDONE\n');

function v = gf(I,f), if isstruct(I)&&isfield(I,f), v=I.(f); else, v=-1; end, end
