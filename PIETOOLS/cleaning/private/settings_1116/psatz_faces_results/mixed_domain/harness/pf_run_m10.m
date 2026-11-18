% PF_RUN_M10  The finite-energy form 'FP' (P >= 0, no margin; T'PA + A'PT +
% 2 epneg T'T <= 0: int ||T x(t)||^2 dt <= C ||T x(0)||^2, the Datko
% condition, which the obstruction argument does not exclude) on S17 at
% 0.5 r2* and on D4 at 0.9 pi^2, N-D default storage, facesS and box, dw 0
% and 1, epneg 1e-5 (the executive's finite-energy default) and 1e-2.
% SCRATCH, MMP, 11/16/2026.
C = pf_cfg({'facesS','box'},2);
for ep = [1e-5 1e-2]
    setappdata(0,'PF_EPNEG',ep);
    pf_exec_m('S17_transport_edge_diffusion',{struct('r2',0.5*pi^2)},C,[0 1],{'def'},'light',sprintf('m10_s17_fp_%g',ep),'mosek','FP');
end
setappdata(0,'PF_EPNEG',1e-5);
pf_exec_m('S17_transport_edge_diffusion',{struct('r2',0.5*pi^2)},C(1),[0 1],{'sep11'},'light','m10_s17_fp_sep11','mosek','FP');
rmappdata(0,'PF_EPNEG');
disp('## M10 DONE');
