% PF_RUN_T3  The product term sets (PF_TERMSETS3: the preordering of the box
% quadratics, the faces with the mixed face x box products, the full
% preordering of the faces; each at the per-direction offset 1 and at 0)
% on the mixed suite: T2 first (the direct PIE2PDE test, Q at length 2,
% dw 0, at 0.5 and 0.9 of the margins), then T1 (the lift of the PDE norm,
% bracket [lam, 1.1 lam], 10 steps). Same systems, lengths and tolerances
% as PF_RUN_T1 and PF_RUN_T2. SCRATCH, MMP, 11/16/2026.
S = pf_setup(false);
tau43 = 1.034726;   a26 = pi^2/4 - 0.5;
[optl,names] = pf_termsets3();
qd1 = {struct('mult',1,'int',[1 1],'out',1,'in',1)};
for FRAC = [0.5 0.9]
    sys = {'S3_edge1d_heat2d',struct('r',FRAC*pi^2); 'S8_two_edges',struct('r3',FRAC*pi^2); ...
           'S17_transport_edge_diffusion',struct('r2',FRAC*pi^2); 'S18_lib43_interior_delay',struct('tau',FRAC*tau43,'sgn',-1); ...
           'S26_lib42_ode_heat_delay',struct('a',FRAC*a26)};
    for i = 1:size(sys,1)
        for j = 1:numel(optl)
            try
                pf_qdirect(sys{i,1},{sys{i,2}},qd1,0,sprintf('t3_%s_%s_%g',sys{i,1}(1:3),names{j},FRAC),1e-4,0,optl{j});
            catch ME
                fprintf('@@ T3 ERROR %s %s: %s\n',sys{i,1},names{j},strrep(ME.message,newline,' '));
            end
        end
    end
    fprintf('## T3-T2 DONE frac=%g\n',FRAC);
end
sys = {'S17','S17_transport_edge_diffusion',struct('r2',0.5*pi^2); 'S3','S3_edge1d_heat2d',struct('r',0.9*pi^2); ...
       'S8','S8_two_edges',struct('r3',0.9*pi^2); 'S18','S18_lib43_interior_delay',struct('tau',0.5*tau43,'sgn',-1); ...
       'S26','S26_lib42_ode_heat_delay',struct('a',0.9*a26)};
for i = 1:size(sys,1)
    try
        Sys = mixed_domain_systems(sys{i,2},sys{i,3});
        PIE = exec_pie_sop(initialize(convert(Sys.PDE,{'pie','silent'})));
        [sp,~] = copvar_space_list(PIE.T,'in');
        degs = zeros(1,numel(sp));
        for k = 1:numel(sp),    degs(k) = [0 6 4 3]*(numel(sp{k})==[0 1 2 3])';  end
        fprintf('@@ T1 %s terms=%s\n',sys{i,1},strjoin(names,','));
        pf_poinc(['t3_' sys{i,1}],sys{i,1},PIE,degs,optl,10,1.1);
    catch ME
        fprintf('@@ T1 ERROR %s: %s\n',sys{i,1},strrep(ME.message,newline,' '));
    end
end
disp('## T3 DONE');
