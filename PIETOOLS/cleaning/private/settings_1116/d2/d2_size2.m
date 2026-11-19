% D2_SIZE2: D2_SIZE with the triangular cap (lpi.ineq.cap_inc, 11/18/2026) on the slack cap rungs: (cap 1, none) the control, (1, auto), (0, auto), (2, auto), (NaN, auto). MMP, 11/18/2026.
% D2_SIZE  Size (never solve) the executive's direct form 'Q' at the light
% tier on the 2-D delayed heat equation with its 3-D delay state (library
% case 56, literature item 8) at tau = 0.5 tau*, tau* = 1.034726, for the
% slack cap rungs between the light 3-D entry (cap_delta 0: 475,404
% variables, m 33,025, INFEASIBLE on Sol, job 1851319) and no cap (NaN:
% 6,813,852 variables, m 163,847, OUT_OF_MEMORY at 200 GB on Sol, job
% 1851724). Same configuration as sol_d2q.m otherwise: 3-D term entry
% 'faces1', wcross 0 (tier), Q at lpi.var.deg{3} = 1, own-support slacks,
% base raise 0, max_raise_feas(3) = 0 (one program per rung). The size-only
% LPISOLVE shadow (size/shadow) prints ndv, m and errors; this script
% records '@@ D2SIZE' lines. The NaN rung is the control against the Sol
% log. SCRATCH, MMP, 11/17/2026.
tstart = tic;
R = 'C:/Users/mpeet/OneDrive/Documents/GitHub/PIETOOLS/PIETOOLS';
S = 'C:/Users/mpeet/AppData/Local/Temp/claude/363c499c6332a503/af90e0a1-bcff-4e8d-9a24-9aa4c52491cc/scratchpad/sol1117/size';
cd(R);  pietools_path_update;
addpath([S '/shadow']);
for f = {'PIETOOLS_stability_sop','exec_tools_sop','get_lift_degs','lpi_ineq_sop','poscopvar_direct','poscopvar','lpisettings_sop','settings_sopvar_light','lpivar_cdopvar','PIETOOLS_PDE_Ex_2D_Delayed_Heat_Eq_Fridman'}
    wf = strrep(which(f{1}),'\','/');   assert(startsWith(wf,R),f{1}+" does not resolve inside the repo: "+wf);
end
wl = strrep(which('lpisolve'),'\','/');   assert(startsWith(wl,[S '/shadow']),"lpisolve is not the size-only shadow: "+wl);
fprintf('@@ D2SIZE start %s\n',datestr(now,31));
PDE_t = PIETOOLS_PDE_Ex_2D_Delayed_Heat_Eq_Fridman(0,{'tau=0.5*1.034726;'});
PDE = initialize_PIETOOLS_PDE(PDE_t,true);
PIE = initialize(convert(PDE,{'pie','silent'}));
fprintf('@@ D2SIZE pie built %.0f s: dim %d\n',toc(tstart),PIE.dim);
variants = {1,[]; 1,'auto'; 0,'auto'; 2,'auto'; NaN,'auto'};
for v = 1:size(variants,1)
    cd3 = variants{v,1};   ci = variants{v,2};
    st = lpisettings_sop('light');   st.sos_opts.solver = 'mosek';   st.sop.verbose = true;
    st.lpi.ineq.psatz = {'auto','faces1','faces1'};
    st.sop.cap_delta = [NaN NaN cd3];   st.lpi.ineq.cap_inc = ci;   tag = sprintf('cd=%g capinc=%s',cd3,mat2str(ci));
    t0 = tic;
    try
        PIETOOLS_stability_sop(PIE,st,'Q');
        fprintf('@@ D2SIZE cd=%g: the executive returned without the size error (t=%.0f)\n',cd3,toc(t0));
    catch ME
        tok = regexp(ME.message,'LB_SIZE ndv=(\d+) m=(\d+) ineq=(\d+)','tokens','once');
        if isempty(tok)
            fprintf('@@ D2SIZE cd=%g ERROR %s (t=%.0f)\n',cd3,strrep(ME.message,newline,' '),toc(t0));
        else
            fprintf('@@ D2SIZE cd=%g ndv=%s m=%s ineq=%s t=%.0f\n',cd3,tok{:},toc(t0));
        end
    end
end
fprintf('@@ D2SIZE end %s total %.0f s\n',datestr(now,31),toc(tstart));
disp('## D2SIZE DONE');
