function prog = cx_H2_norm_c(PIE,st,gam)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PROG = CX_H2_NORM_C(PIE,ST,GAM) builds, with copvar/cdopvar only, the LPI
% of PIETOOLS_H2_norm_c (1-D, controllability gramian, NON-coercive Q-form)
% at the FIXED bound GAM, and returns it unsolved:
%
%   R = T*Q >= 0 (R positive, Q free),  W >= 0 on the z space,
%   Dneg = [-gam*Iw, B1'; B1, Q'*A' + A*Q] <= 0,
%   Dpos = [W, C1*Q; Q'*C1', R] >= 0,        gam - trace(W.P) >= 0.
%
% GAM is the executive's decision variable gam, i.e. the H2 bound itself
% (the stock output is gam, not sqrt, L195) - the quantity
% cx_stock_pinned pins. Mirrors the stock options.h2 branch (L112).
%
% Stock lines mirrored (executives/PIETOOLS_H2_norm_c.m):
%   L66-73    operators, Tw == 0 check
%   L76-95    settings read exactly as the stock
%   L100      lpiprogram
%   L123-130  R = poslpivar(dd1,options1) [+ dd12 psatz]   -> cx_h2_lf
%   L133-135  Q = lpivar(get_lpivar_degs(R,T)), T*Q - R = 0 -> cx_h2_qdeg,
%             lpivar_cdopvar, lpi_eq_cdopvar (non-symmetric)
%   L139      W = poslpivar(C1.dim(:,1)), default degrees  -> cx_h2_pos
%   L146-151  Iw, Dneg, Dpos                                -> container
%             concatenation; spaces kept separate
%   L166-179  slacks at dd2 [+ dd3 psatz], lpi_eq 'symmetric' -> cx_h2_slack
%   L183-184  gam - trace(W.P) >= 0 via the stock scalar lpi_ineq
%             -> cx_h2_trace
% NOT mirrored: L106-110 (gam as decision variable, gam >= 0, objective):
% a dpvar cannot scale a container (@cdopvar/mtimes.m:121-124), so gam is
% fixed and bisected outside; L160-163 (sosineq_on): no container lpi_ineq,
% error; L55-63 (2-D dispatch): see cx_H2_norm_2D_c.
%
% Initial coding MMP, 09/25/2026
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if nargin<3,    error('cx_H2_norm_c:gam','The container LPI is posed at fixed gamma.'),  end
if ~isa(PIE,'pie_struct'),  error('cx_H2_norm_c:pie','PIE must be a pie_struct.'),  end
PIE = initialize(PIE);
if PIE.dim==2,  error('cx_H2_norm_c:dim','2-D: use cx_H2_norm_2D_c.'),  end
if ~(PIE.Tw==0)                                                             % L71-73
    error('cx_H2_norm_c:Tw','H2 norm LPI cannot be solved with disturbances at the boundary.')
end
if st.sosineq_on                                                            % L87-88, L160-163
    error('cx_H2_norm_c:sosineq','sosineq_on: no container lpi_ineq (library gap).')
end

Tm = opvar2copvar(PIE.T);   Am = opvar2copvar(PIE.A);
Bw = opvar2copvar(PIE.B1);  Cz = opvar2copvar(PIE.C1);
% L146: Iw = mat2opvar(eye(size(B1op,2)),B1op.dim(:,2),...), then converted
Iw = opvar2copvar(mat2opvar(eye(size(PIE.B1,2)),PIE.B1.dim(:,2),PIE.vars,PIE.dom));

prog = lpiprogram(PIE.vars(:,1),PIE.vars(:,2),PIE.dom);                     % L100
[sp,dm] = cx_space_list(Tm,'out');          % state spaces, zero-dim dropped
[prog,Rm] = cx_h2_lf(prog,dm,sp,PIE.dom,st);                                % L123-130

% L133-135: Q free with lpivar's degree family; T*Q - R = 0, NOT symmetric.
[prog,Qm] = lpivar_cdopvar(prog,dm,sp,PIE.dom,cx_h2_qdeg(PIE,st));
prog = lpi_eq_cdopvar(prog,Tm*Qm - Rm);

% L139: poslpivar(prog,C1op.dim(:,1)) takes the nargin==2 defaults (psatz 0,
% d = {1,[1,1,1],[1,1,1]}); on a finite z it is a PSD nz x nz matrix.
[spz,dmz] = cx_space_list(Cz,'out');
[prog,Wm] = cx_h2_pos(prog,dmz,spz,PIE.dom,[],struct());

% L148-151. gam is numeric, so -(gam*Iw) is the scalar branch of mtimes.
Dneg = [-(gam*Iw),  Bw';
        Bw,         Qm'*Am' + Am*Qm];
Dpos = [Wm,         Cz*Qm;
        Qm'*Cz',    Rm];

% L166-179: Dneg = -Deop, Dpos = Deopp, both 'symmetric'.
[prog,Nneg] = cx_h2_slack(prog,Dneg,PIE.dom,st);
[prog,Npos] = cx_h2_slack(prog,Dpos,PIE.dom,st);
prog = lpi_eq_cdopvar(prog,Nneg + Dneg,'symmetric');
prog = lpi_eq_cdopvar(prog,Npos - Dpos,'symmetric');

% L183-184: gam >= trace(W.P), the stock scalar lpi_ineq on a dpvar.
prog = lpi_ineq(prog,gam - cx_h2_trace(Wm));
end
