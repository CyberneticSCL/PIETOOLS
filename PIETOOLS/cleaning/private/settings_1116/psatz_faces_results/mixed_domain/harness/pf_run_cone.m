% PF_RUN_CONE  Fixed operators against the storage cones (PF_CONE): the
% PDE-norm forms T'*T of S17, S3, S18 and of S17's x2 alone (local in s1,
% Volterra kernel in s2), the certified S18 storage P* and derivative
% -D*, and the certified S3 storage (uniform, control). Variants: uniform
% ('def') at the light and the raised degrees; every separable basis
% operator (sep11) at both; every basis operator with split kernels
% (sep00) at the raised degrees. SCRATCH, MMP, 11/16/2026.
S = pf_setup(false);
pvar s1 s2 t
tau43 = 1.034726;
light  = struct('mult',[1 0 0],'int',[1 1 1],'joint',[2 2 2],'sep',0);
raised = struct('mult',[1 1 1],'int',[2 2 2],'joint',[4 4 4],'sep',0);
variants = {{'def',[0 0],light,'light'},{'def',[0 0],raised,'raised'},{'sep',[1 1],light,'light'},{'sep',[1 1],raised,'raised'},{'sep',[0 0],raised,'raised'}};
ops = struct('name',{},'X',{});
sys = {'S17_transport_edge_diffusion',struct('r2',0.5*pi^2); 'S3_edge1d_heat2d',struct('r',0.9*pi^2); 'S18_lib43_interior_delay',struct('tau',0.5*tau43,'sgn',-1)};
for i = 1:size(sys,1)
    Sys = mixed_domain_systems(sys{i,1},sys{i,2});
    PIE = exec_pie_sop(initialize(convert(Sys.PDE,{'pie','silent'})));
    ops(end+1) = struct('name',[sys{i,1}(1:3) '_TtT'],'X',PIE.T'*PIE.T);  %#ok<AGROW>
    if i==3,    T18 = PIE.T;    A18 = PIE.A;  end
    if i==2,    PIE3 = PIE;     end
end
clear stateNameGenerator
x2 = pde_var(1,[s1;s2],[0,1;0,1],[0;2]);
PDE = [diff(x2,t)==diff(x2,s2,2); subs(x2,s2,0)==0; subs(x2,s2,1)==0];
PIE = exec_pie_sop(initialize(convert(PDE,{'pie','silent'})));
ops(end+1) = struct('name','S17x2_TtT','X',PIE.T'*PIE.T);
L = load(sprintf('%s/results/s18_Pstar.mat',S));
ops(end+1) = struct('name','S18_Pstar','X',L.Pst);
ops(end+1) = struct('name','S18_negD','X',-(T18'*L.Pst*A18 + A18'*L.Pst*T18));
pf_setup('route');
st3 = lpisettings_sop('light');  st3.sos_opts.solver = 'mosek';  st3.sop.verbose = false;
st3.sop.dw_feas = [1 1 1];  st3.sop.max_raise_feas = [0 0 0];
[~,P3,info3] = PIETOOLS_stability_sop(PIE3,st3,'P');
fprintf('@@ CONE S3 storage solve: %s cert=%d ndv=%d\n',info3.status,info3.certified,info3.shape.ndv);
pf_setup(false);
ops(end+1) = struct('name','S3_Pstar','X',P3);
pf_cone('mixed',ops,variants);
disp('## CONE RUN DONE');
