function prog = cx_H2_norm_o_coercive(PIE,st,gam)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PROG = CX_H2_NORM_O_COERCIVE(PIE,ST,GAM) builds, with copvar/cdopvar
% only, the LPI of PIETOOLS_H2_norm_o_coercive (1-D, observability gramian,
% coercive) at the FIXED value GAM, and returns it unsolved:
%
%   W >= 0,   A'*W*T + T'*W*A + C1'*C1 <= 0,   gam - trace((B1'*W*B1).P) >= 0.
%
% GAM is the executive's decision variable gam, which bounds the SQUARED
% H2 norm (the stock prints sqrt(gam), L172) - the quantity cx_stock_pinned
% pins. Mirrors the stock options.h2 branch (L106), with gam = GAM.
%
% Stock lines mirrored (executives/PIETOOLS_H2_norm_o_coercive.m):
%   L60-67    operators, Tw == 0 check;   L69-88 settings;   L94 lpiprogram
%   L117-124  W = poslpivar(dd1,options1) [+ dd12 psatz]    -> cx_h2_lf
%   L131      Dop = (A'*W)*T + T'*(W*A) + C1'*C1            -> same products
%   L145-153  slack at dd2 [+ dd3 psatz], lpi_eq 'symmetric' -> cx_h2_slack
%   L156-163  trace of (B1'*W*B1).P, gam - trace >= 0       -> cx_h2_trace,
%             stock scalar lpi_ineq
% NOT mirrored: L100-104 (gam decision variable, gam >= 0, objective);
% L140-142 (sosineq_on: no container lpi_ineq, error).
% With a distributed w, B1'*W*B1 has no R^n part: the stock options.h2
% branch then passes a double to lpi_ineq, which errors (lpi_ineq.m:74-75);
% here the trace is a constant 0 dpvar and the constraint is gam >= 0.
%
% Initial coding MMP, 09/25/2026
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if nargin<3,    error('cx_H2_norm_o_coercive:gam','The container LPI is posed at fixed gamma.'),  end
if ~isa(PIE,'pie_struct'),  error('cx_H2_norm_o_coercive:pie','PIE must be a pie_struct.'),  end
PIE = initialize(PIE);
if PIE.dim==2,  error('cx_H2_norm_o_coercive:dim','2-D: use cx_H2_norm_2D_o.'),  end
if ~(PIE.Tw==0)                                                             % L65-67
    error('cx_H2_norm_o_coercive:Tw','H2 norm LPI cannot be solved with disturbances at the boundary.')
end
if st.sosineq_on                                                            % L80-81, L140-142
    error('cx_H2_norm_o_coercive:sosineq','sosineq_on: no container lpi_ineq (library gap).')
end

Tm = opvar2copvar(PIE.T);   Am = opvar2copvar(PIE.A);
Bw = opvar2copvar(PIE.B1);  Cz = opvar2copvar(PIE.C1);

prog = lpiprogram(PIE.vars(:,1),PIE.vars(:,2),PIE.dom);                     % L94
[sp,dm] = cx_space_list(Tm,'out');
[prog,Wm] = cx_h2_lf(prog,dm,sp,PIE.dom,st);                                % L117-124

Dm = (Am'*Wm)*Tm + Tm'*(Wm*Am) + Cz'*Cz;                                    % L131

[prog,Nm] = cx_h2_slack(prog,Dm,PIE.dom,st);                                % L145-152
prog = lpi_eq_cdopvar(prog,Nm + Dm,'symmetric');                            % L153

% L156-163: tempObj = Bwop'*Wop*Bwop evaluates left to right.
prog = lpi_ineq(prog,gam - cx_h2_trace((Bw'*Wm)*Bw));
end
