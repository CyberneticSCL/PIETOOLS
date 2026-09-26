function prog = cx_H2_norm_c_coercive(PIE,st,gam)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PROG = CX_H2_NORM_C_COERCIVE(PIE,ST,GAM) builds, with copvar/cdopvar
% only, the LPI of PIETOOLS_H2_norm_c_coercive (1-D, controllability
% gramian, coercive) at the FIXED value GAM, and returns it unsolved:
%
%   W >= 0,   A*W*T' + T*W*A' + B1*B1' <= 0,   gam - trace((C1*W*C1').P) >= 0.
%
% GAM is the executive's decision variable gam, which bounds the SQUARED
% H2 norm (the stock prints sqrt(gam), L171) - the quantity cx_stock_pinned
% pins. Mirrors the stock options.h2 branch (L107), with gam = GAM.
%
% Stock lines mirrored (executives/PIETOOLS_H2_norm_c_coercive.m):
%   L61-68    operators, Tw == 0 check;   L71-90 settings;   L95 lpiprogram
%   L118-125  W = poslpivar(dd1,options1) [+ dd12 psatz]    -> cx_h2_lf
%   L132      Dop = (A*W)*T' + T*(W*A') + B1*B1'            -> same products
%   L146-154  slack at dd2 [+ dd3 psatz], lpi_eq 'symmetric' -> cx_h2_slack
%   L157-164  trace of (C1*W*C1').P, gam - trace >= 0       -> cx_h2_trace,
%             stock scalar lpi_ineq
% NOT mirrored: L101-105 (gam decision variable, gam >= 0, objective);
% L141-143 (sosineq_on: no container lpi_ineq, error).
% When (C1*W*C1') has no R^n part the stock options.h2 branch passes a
% double to lpi_ineq, which errors (lpi_ineq.m:74-75); here the trace is a
% constant 0 dpvar and the constraint is gam >= 0 (cx_h2_trace).
%
% Initial coding MMP, 09/25/2026
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if nargin<3,    error('cx_H2_norm_c_coercive:gam','The container LPI is posed at fixed gamma.'),  end
if ~isa(PIE,'pie_struct'),  error('cx_H2_norm_c_coercive:pie','PIE must be a pie_struct.'),  end
PIE = initialize(PIE);
if PIE.dim==2,  error('cx_H2_norm_c_coercive:dim','2-D: use cx_H2_norm_2D_c.'),  end
if ~(PIE.Tw==0)                                                             % L66-68
    error('cx_H2_norm_c_coercive:Tw','H2 norm LPI cannot be solved with disturbances at the boundary.')
end
if st.sosineq_on                                                            % L82-83, L141-143
    error('cx_H2_norm_c_coercive:sosineq','sosineq_on: no container lpi_ineq (library gap).')
end

Tm = opvar2copvar(PIE.T);   Am = opvar2copvar(PIE.A);
Bw = opvar2copvar(PIE.B1);  Cz = opvar2copvar(PIE.C1);

prog = lpiprogram(PIE.vars(:,1),PIE.vars(:,2),PIE.dom);                     % L95
[sp,dm] = cx_space_list(Tm,'out');
[prog,Wm] = cx_h2_lf(prog,dm,sp,PIE.dom,st);                                % L118-125

Dm = (Am*Wm)*Tm' + Tm*(Wm*Am') + Bw*Bw';                                    % L132

[prog,Nm] = cx_h2_slack(prog,Dm,PIE.dom,st);                                % L146-153
prog = lpi_eq_cdopvar(prog,Nm + Dm,'symmetric');                            % L154

% L157-164: tempObj = Czop*Wop*Czop' evaluates left to right.
prog = lpi_ineq(prog,gam - cx_h2_trace((Cz*Wm)*Cz'));
end
