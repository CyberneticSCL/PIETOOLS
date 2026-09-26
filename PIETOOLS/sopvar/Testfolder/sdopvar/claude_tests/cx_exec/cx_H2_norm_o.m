function prog = cx_H2_norm_o(PIE,st,gam)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PROG = CX_H2_NORM_O(PIE,ST,GAM) builds, with copvar/cdopvar only, the LPI
% of PIETOOLS_H2_norm_o (1-D, observability gramian, NON-coercive Q-form)
% at the FIXED bound GAM, and returns it unsolved:
%
%   R = T'*Q >= 0 (R positive, Q free),  W >= 0 on the w space,
%   Dneg = [-gam*Iz, C1; C1', Q'*A + A'*Q] <= 0,
%   Dpos = [W, B1'*Q; Q'*B1, R] >= 0,        gam - trace(W.P) >= 0.
%
% GAM is the executive's decision variable gam, the H2 bound itself (the
% stock output is gam, not sqrt, L191) - the quantity cx_stock_pinned pins.
% Mirrors the stock options.h2 branch (L109).
%
% Stock lines mirrored (executives/PIETOOLS_H2_norm_o.m):
%   L65-72    operators, Tw == 0 check;   L74-93 settings;   L99 lpiprogram
%   L120-127  R = poslpivar(dd1,options1) [+ dd12 psatz]   -> cx_h2_lf
%   L130-132  Q = lpivar(get_lpivar_degs(R,T)), T'*Q - R = 0 (non-symmetric)
%   L136      W = poslpivar(B1.dim(:,2)), default degrees  -> cx_h2_pos over
%             B1's INPUT spaces; an L2 operator when w is distributed
%   L143-148  Iz, Dneg, Dpos                                -> concatenation
%   L162-175  slacks at dd2 [+ dd3 psatz], lpi_eq 'symmetric' -> cx_h2_slack
%   L179-180  gam - trace(W.P) >= 0, stock scalar lpi_ineq  -> cx_h2_trace
%             (a constant 0 when w has no R^n part: the stock trace of a
%             0x0 dpvar is likewise a 1x1 dpvar)
% NOT mirrored: L105-107 (gam decision variable and objective); L156-159
% (sosineq_on: no container lpi_ineq, error); L54-62 (2-D dispatch).
%
% Initial coding MMP, 09/25/2026
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if nargin<3,    error('cx_H2_norm_o:gam','The container LPI is posed at fixed gamma.'),  end
if ~isa(PIE,'pie_struct'),  error('cx_H2_norm_o:pie','PIE must be a pie_struct.'),  end
PIE = initialize(PIE);
if PIE.dim==2,  error('cx_H2_norm_o:dim','2-D: use cx_H2_norm_2D_o.'),  end
if ~(PIE.Tw==0)                                                             % L70-72
    error('cx_H2_norm_o:Tw','H2 norm LPI cannot be solved with disturbances at the boundary.')
end
if st.sosineq_on                                                            % L85-86, L156-159
    error('cx_H2_norm_o:sosineq','sosineq_on: no container lpi_ineq (library gap).')
end

Tm = opvar2copvar(PIE.T);   Am = opvar2copvar(PIE.A);
Bw = opvar2copvar(PIE.B1);  Cz = opvar2copvar(PIE.C1);
% L143: Iz = mat2opvar(eye(size(C1op,1)),C1op.dim(:,1),...), then converted
Iz = opvar2copvar(mat2opvar(eye(size(PIE.C1,1)),PIE.C1.dim(:,1),PIE.vars,PIE.dom));

prog = lpiprogram(PIE.vars(:,1),PIE.vars(:,2),PIE.dom);                     % L99
[sp,dm] = cx_space_list(Tm,'out');
[prog,Rm] = cx_h2_lf(prog,dm,sp,PIE.dom,st);                                % L120-127

% L130-132: T'*Q - R = 0, NOT symmetric.
[prog,Qm] = lpivar_cdopvar(prog,dm,sp,PIE.dom,cx_h2_qdeg(PIE,st));
prog = lpi_eq_cdopvar(prog,Tm'*Qm - Rm);

% L136: nargin==2 defaults, on the disturbance (input) spaces of B1.
[spw,dmw] = cx_space_list(Bw,'in');
[prog,Wm] = cx_h2_pos(prog,dmw,spw,PIE.dom,[],struct());

Dneg = [-(gam*Iz),  Cz;                                                     % L145-146
        Cz',        Qm'*Am + Am'*Qm];
Dpos = [Wm,         Bw'*Qm;                                                 % L147-148
        Qm'*Bw,     Rm];

[prog,Nneg] = cx_h2_slack(prog,Dneg,PIE.dom,st);                            % L162-173
[prog,Npos] = cx_h2_slack(prog,Dpos,PIE.dom,st);
prog = lpi_eq_cdopvar(prog,Nneg + Dneg,'symmetric');                        % L174
prog = lpi_eq_cdopvar(prog,Npos - Dpos,'symmetric');                        % L175

prog = lpi_ineq(prog,gam - cx_h2_trace(Wm));                                % L179-180
end
