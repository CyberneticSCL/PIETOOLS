% hinf_ref.m -- H-infinity gain reference, and VALIDATE that I can recover
% gamma from the raw solution vector before trusting any cuADMM answer.
%
% gamma is read by lpigetsol from solinfo.RRx (decvartable order), NOT from
% solinfo.x (the solver's cone vector).  Independently: sossolve puts
% sos.objective on the leading block of c, so the objective's index is
% find(prog.objective ~= 0).  If RRx at that index equals lpigetsol's gamma,
% the extraction is trustworthy.
%
% Also records the SDP shape and, critically for cuADMM, whether c is nonzero
% -- this is the FIRST program in this campaign that is an optimisation rather
% than a feasibility problem.
cuadmm_path;
pvar s t
x = pde_var('state',1,s,[0,1]);
w = pde_var('input',1);
z = pde_var('output',1);
PIE = convert([diff(x,t,1)==diff(x,s,2)+2*x+s*w;
               z == int(x,s,[0,1]);
               subs(x,s,0)==0; subs(x,s,1)==0]);

for pre = {'light','heavy'}
  for slv = {'sedumi','mosek'}
    st = lpisettings(pre{1});  st.sos_opts.solver = slv{1};
    try
        t0=tic; [prog,~,gam] = PIETOOLS_Hinf_gain(PIE,st); tw=toc(t0);
        g_exec = double(gam);

        % independent recovery
        oi = find(prog.objective ~= 0);
        RRx = prog.solinfo.RRx(:);
        g_mine = RRx(oi);

        Atf=[];bf=[];
        for i=1:prog.expr.num, Atf=[Atf,prog.expr.At{i}]; bf=[bf;prog.expr.b{i}]; end
        S = cuadmm_private('sdpshape',prog);
        rb = norm(full(Atf'*RRx-bf))/norm(full(bf));
        I = prog.solinfo.info;
        fprintf(['HR %-5s %-6s gam_exec=%.8f gam_RRx=%.8f match=%.2e | m=%d Kf=%d ' ...
                 'Ns=[%s] |b|=%.3e obj_nnz=%d obj_idx=%d | rel_b=%.3e numerr=%d t=%.2f\n'], ...
            pre{1},slv{1},g_exec,g_mine,abs(g_exec-g_mine)/max(abs(g_exec),eps), ...
            S.m,S.Kf,strtrim(num2str(S.Ks)),norm(full(bf)),nnz(prog.objective),oi(1), ...
            rb,I.numerr,tw);
    catch ME
        fprintf('HR %-5s %-6s FAILED %s\n',pre{1},slv{1},strrep(ME.message,newline,' '));
    end
  end
end
fprintf('HRDONE\n');
