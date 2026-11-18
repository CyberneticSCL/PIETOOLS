% PF_RUN_T0  Phase 0 of the term run: the row reader against known results.
% (a) the row storages of S17's PIE: basis lists and Gram sizes for Rdef,
% Rlight0, Rlight, Rlight2, Rheavy, against POSCOPVAR with DEG omitted;
% (b) S18 at 0.5 tau* in the 'P' form with Rlight2 (the storage that
% certified it at 767k variables as 'sep11 raised'); (c) the 0-D-1-D system
% S1 (ODE in a Dirichlet BC of a 1-D heat) with Rlight against the stock
% storage: verdict, variables, m. SCRATCH, MMP, 11/16/2026.
S = pf_setup(false);
tau43 = 1.034726;
Sys = mixed_domain_systems('S17_transport_edge_diffusion',struct('r2',0.5*pi^2));
PIE = exec_pie_sop(initialize(convert(Sys.PDE,{'pie','silent'})));
T = PIE.T;  dom = struct('vars',{T.vars},'dom',T.dom);
st = lpisettings_sop('light');
prog = lpiprogram(T.vars(:),[],T.dom,st);
[sp,dm] = copvar_space_list(T,'out');
[p0,P0] = poscopvar(prog,dm,sp,dom);
fprintf('@@ T0 storage def: ndv=%d\n',size(p0.decvartable,1));
for nm = {'Rdef','Rlight0','Rlight','Rlight2','Rheavy'}
    try
        [p1,P1,info] = pf_row_storage(prog,T,'out',dom,pf_rows(nm{1}));
        fprintf('@@ T0 storage %s: ndv=%d sep=%s basis=%s\n',nm{1},size(p1.decvartable,1),mat2str(info.sep),mat2str(info.basis_list));
    catch ME
        fprintf('@@ T0 storage %s: ERROR %s\n',nm{1},strrep(ME.message,newline,' '));
    end
end
% (b) S18 with Rlight2 in the 'P' form, facesS, dw 0
setappdata(0,'PF_ROW',pf_rows('Rlight2'));
C = pf_cfg({'facesS'},2);
pf_exec_m('S18_lib43_interior_delay',{struct('tau',0.5*tau43,'sgn',-1)},C,0,{'row'},'light','t0_s18_Rlight2');
% (c) S1 (R x L2[s1]) with Rlight against stock, 'P' form, dw 0 and 1, the tier's own terms
setappdata(0,'PF_ROW',pf_rows('Rlight'));
pf_exec_m('S1_ode_heat1d_bc',{struct('a',0.5)},C,[0 1],{'row','stock'},'light','t0_s1');
rmappdata(0,'PF_ROW');
disp('## T0 DONE');
