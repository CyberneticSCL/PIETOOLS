% gate_popcheck.m -- does an operator gate normalised by ||Pop|| accept a point
% the SDP solver has certified INFEASIBLE?
%
% An audit agent reported it does, on the stock PIE2PDEstability build at
% lam = 50. That is a claim about MY panel: lpi_resid gates on op_indP =
% ||Dop+Deop|| / ||Pop||, adopted precisely because ||Pop|| carries the eppos
% constant and so cannot collapse to zero the way ||Dop|| does. If the report is
% right, the denominator recommendation I sent to the low-rank session is not
% merely incomplete but actively wrong.
%
% Campaign rule: reproduce a reported defect before acting on it. The sweep runs
% ACROSS the Poincare threshold lam* = pi^2, where feasibility is settled by
% theorem rather than by solver opinion.
%
% stab_mirror is used rather than the stock executive because the quantity under
% test is ||Dop + Deop||, and Deop -- the positive slack the executive adds -- is
% internal to PIETOOLS_PIE2PDEstability, which returns only (prog,P). The mirror
% reproduces that build line for line and was checked to pose the same SDP
% (matching m, Kf, Ns, nnz). Rebuilding the residual from the returned P alone
% would silently make op_ind identically 1.
cuadmm_path;
LAMSTAR = pi^2;
LAMS = [0.5 0.9 1.0 1.05 1.5 5.0]*LAMSTAR;
THR  = 1e-4;
fprintf('GP lam_over_lamstar|numerr|pinf|dinf|rel_b|op_ind|op_indP|normRes|normDop|normPop|normx|acc_ind|acc_indP|acc_relb\n');
for lam = LAMS
    clear stateNameGenerator
    pvar s t
    x = pde_var(1,s,[0,1]);
    PIE = initialize(convert([diff(x,t,1)==diff(x,s,2)+lam*x; ...
                              subs(x,s,0)==0; subs(x,s,1)==0]));
    st = lpisettings('light');
    st.sos_opts.solver='mosek'; st.sos_opts.simplify=false;
    try
        out = cuadmm_private('stab_mirror',PIE,st);
        evalc('sol = lpisolve(out.prog,st.sos_opts);');
        Atf=[];bf=[];
        for q=1:sol.expr.num, Atf=[Atf,sol.expr.At{q}]; bf=[bf;sol.expr.b{q}]; end
        xr = sol.solinfo.RRx(:);
        rb = norm(full(Atf'*xr-bf))/max(norm(full(bf)),eps);
        I  = sol.solinfo.info;

        Ds = getsol_lpivar(sol,out.Dop);
        Es = getsol_lpivar(sol,out.Deop);
        Ps = getsol_lpivar(sol,out.Pop);
        GR = cuadmm_private('opnorm_pi',Ds+Es,16);        % the operator that must be zero
        GD = cuadmm_private('opnorm_pi',Ds,16);
        GP = cuadmm_private('opnorm_pi',Ps,16);
        oi  = norm(GR)/max(norm(GD),eps);
        oip = norm(GR)/max(norm(GP),eps);
        fprintf('GP %.3f|%d|%d|%d|%.4e|%.4e|%.4e|%.4e|%.4e|%.4e|%.4e|%d|%d|%d\n', ...
            lam/LAMSTAR,gf(I,'numerr'),gf(I,'pinf'),gf(I,'dinf'), ...
            rb,oi,oip,norm(GR),norm(GD),norm(GP),norm(xr), ...
            oi<THR,oip<THR,rb<THR);
    catch ME
        fprintf('GP %.3f|ERR|%s\n',lam/LAMSTAR,regexprep(ME.message,'[\r\n\t]+',' '));
    end
end
fprintf('GPDONE\n');

function v = gf(I,f), if isstruct(I)&&isfield(I,f), v=I.(f); else, v=-1; end, end
