% PF_RUN_Q  The PIE2PDE form ('Q': P >= eppos T'*T, T'*Q = P, A'*Q + Q'*A <= 0)
% on the mixed-domain systems, which the repo executive poses as 'P' in two
% or more variables; the qform shadow allows it. S3 at 0.9 r* (control),
% S17 at 0.5 r2*, S18 at 0.5 tau*; light, N-D default storage, then the
% separable full-operator storage at the raised degrees; facesS; dw 0, 1.
% SCRATCH, MMP, 11/16/2026.
C = pf_cfg({'facesS'},2);
tau43 = 1.034726;
pf_exec_m('S3_edge1d_heat2d',{struct('r',0.9*pi^2)},C,[0 1],{'def'},'light','q_s3','mosek','Q');
pf_exec_m('S17_transport_edge_diffusion',{struct('r2',0.5*pi^2)},C,[0 1],{'def'},'light','q_s17','mosek','Q');
pf_exec_m('S18_lib43_interior_delay',{struct('tau',0.5*tau43,'sgn',-1)},C,[0 1],{'def'},'light','q_s18','mosek','Q');
setappdata(0,'PF_POS',struct('mult',[1 1 1],'int',[2 2 2],'joint',[4 4 4],'sep',0));
pf_exec_m('S17_transport_edge_diffusion',{struct('r2',0.5*pi^2)},C,[0 1],{'sep11'},'light','q_s17_sep11','mosek','Q');
pf_exec_m('S18_lib43_interior_delay',{struct('tau',0.5*tau43,'sgn',-1)},C,0,{'sep11'},'light','q_s18_sep11','mosek','Q');
rmappdata(0,'PF_POS');
disp('## Q DONE');
