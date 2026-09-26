function prog = cx_H2_norm_2D_c(PIE,st,gam)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PROG = CX_H2_NORM_2D_C(PIE,ST,GAM) builds, with copvar/cdopvar only, the
% LPI of PIETOOLS_H2_norm_2D_c (2-D, controllability gramian, coercive) at
% the FIXED value GAM, and returns it unsolved:
%
%   W >= 0 (+ eppos I),  A*W*T' + T*W*A' + B1*B1' <= 0,
%   gam - trace((C1*W*C1').R00) >= 0.
%
% GAM is the executive's decision variable gam, the SQUARED bound (stock
% options.h2 branch: gam = options.h2^2, L136; output sqrt(gam), L212).
% For structure comparison only in this suite: 2-D programs are not solved.
%
% Stock lines mirrored (executives/2D/PIETOOLS_H2_norm_2D_c.m):
%   L58-65    operators, Tw == 0 check
%   L67-122   settings (settings_2d, eppos, LF/eq degrees, options and
%             psatz terms)                                  -> cx_h2_set2d
%   L125      lpiprogram
%   L147-155  W = poslpivar_2d(LF_deg,LF_opts) [+ psatz]    -> cx_h2_pos2d
%   L158-163  W += opvar2d(blkdiag(eppos(1)I,0,0,eppos(4)I)) -> plus of its
%             opvar2d2copvar image
%   L170      Dop = (A*W)*T' + T*(W*A') + B1*B1'
%   L185-196  slack at eq_deg [+ psatz], lpi_eq_2d 'symmetric' -> cx_h2_pos2d,
%             lpi_eq_cdopvar
%   L199-203  trace of (C1*W*C1').R00, gam - trace >= 0      -> cx_h2_trace,
%             stock scalar lpi_ineq
% NOT mirrored: L130-134 (gam as objective variable: gam is fixed);
% L180-182 (use_sosineq, error). States with R, L2[x] or L2[y] parts are
% refused (cx_h2_pos2d): the settings translation covers L2[x,y] only, and
% opvar2d2copvar has no zero-fill for empty rows (map_container_caps).
%
% Initial coding MMP, 09/25/2026
% MMP, 09/26/2026: opvar2d2copvar has zero-fill for empty rows since
%                  09/26/2026; the L2[x,y]-only scope remains for the settings
%                  translation.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if nargin<3,    error('cx_H2_norm_2D_c:gam','The container LPI is posed at fixed gamma.'),  end
if ~isa(PIE,'pie_struct'),  error('cx_H2_norm_2D_c:pie','PIE must be a pie_struct.'),  end
PIE = initialize(PIE);                                                      % L55
if ~(PIE.Tw==0)                                                             % L63-65
    error('cx_H2_norm_2D_c:Tw','H2 norm LPI cannot be solved with disturbances at the boundary.')
end
[X,eppos] = cx_h2_set2d(st);                                                % L67-122

Tm = opvar2d2copvar(PIE.T);     Am = opvar2d2copvar(PIE.A);
Bw = opvar2d2copvar(PIE.B1);    Cz = opvar2d2copvar(PIE.C1);

prog = lpiprogram(PIE.vars(:,1),PIE.vars(:,2),PIE.dom);                     % L125
[sp,dm] = cx_space_list(Tm,'out');
[prog,Wm] = cx_h2_pos2d(prog,sp,dm,PIE.dom,X.LF,X.LF_psatz);                % L147-155
if ~all(eppos==0)                                                           % L158-163
    np = PIE.T.dim(:,1);
    Ip = blkdiag(eppos(1)*eye(np(1)),zeros(np(2)),zeros(np(3)),eppos(4)*eye(np(4)));
    Wm = Wm + opvar2d2copvar(opvar2d(Ip,PIE.T.dim,PIE.dom,PIE.vars));
end

Dm = (Am*Wm)*Tm' + Tm*(Wm*Am') + Bw*Bw';                                    % L170

[spd,dmd] = cx_space_list(Dm,'out');
[prog,Nm] = cx_h2_pos2d(prog,spd,dmd,PIE.dom,X.eq,X.eq_psatz);              % L185-195
prog = lpi_eq_cdopvar(prog,Nm + Dm,'symmetric');                            % L196

prog = lpi_ineq(prog,gam - cx_h2_trace((Cz*Wm)*Cz'));                       % L199-203
end
