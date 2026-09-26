function prog = cx_H2_control(PIE,st,gam)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PROG = CX_H2_CONTROL(PIE,ST,GAM) builds, with copvar/cdopvar only, the LPI
% of PIETOOLS_H2_control (1-D H2-optimal full-state feedback synthesis) at
% the FIXED bound GAM, and returns it unsolved:
%
%   P >= 0 (+ eppos margin), Z free (state -> u), W free (z -> z),
%   Dneg = [-gam*Iw, B1'; B1, T*P*A' + A*P*T' + B2*Z*T' + T*Z'*B2'] <= 0,
%   Dpos = [W, C1*P + D12*Z; (C1*P + D12*Z)', P] >= 0,
%   gam - trace(W.P) >= 0.
%
% GAM is the executive's decision variable gam, the closed-loop H2 bound
% (stock output gam, L209) - the quantity cx_stock_pinned pins. gam enters
% both -gam*Iw and the trace constraint, so the program cannot be recast as
% min trace(W).
%
% Stock lines mirrored (executives/PIETOOLS_H2_control.m):
%   L65-90    operators; 2-D, Tw, Tu, D11 checks
%   L99-118   settings read exactly as the stock
%   L122-124  lpiprogram, nx1, nx2
%   L142-148  P = poslpivar([nx1;nx2],dd1,options1) [+ dd12] -> cx_h2_lf
%   L150-151  P.P += eppos*I, P.R.R0 += eppos2*I   -> plus of the multiplier
%             blkdiag(eppos*I,eppos2*I) (containers have no component fields)
%   L153-154  Z = lpivar([B2.dim(:,2),T.dim(:,1)],ddZ)   -> lpivar_cdopvar,
%             rectangular, legacy [d1 d2 d3] degrees
%   L156-157  W = lpivar([C1.dim(:,1),C1.dim(:,1)],ddZ)  -> lpivar_cdopvar
%   L160-166  Iw, Dneg, Dp12, Dpos; products associate left to right as in
%             MATLAB; D12 restated on Z's registry (cx_on_registry), since
%             container mtimes demands identical registries
%   L178-191  slacks at dd2 [+ dd3 psatz], lpi_eq 'symmetric' -> cx_h2_slack
%   L167,194  gam - trace(W.P) >= 0 via the stock scalar lpi_ineq
%             -> cx_h2_trace
% NOT mirrored: L129-130, L195, L201 (gam decision variable, gam >= 0,
% objective: gam is fixed); L172-175 (sosineq_on, no container lpi_ineq,
% error); L206-214 (getsol / getController: no container getsol).
%
% Initial coding MMP, 09/25/2026
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if nargin<3,    error('cx_H2_control:gam','The container LPI is posed at fixed gamma.'),  end
if ~isa(PIE,'pie_struct'),  error('cx_H2_control:pie','PIE must be a pie_struct.'),  end
PIE = initialize(PIE);
if PIE.dim==2,          error('cx_H2_control:dim','2-D PIEs are not supported (L76-78).'),  end
if ~(PIE.Tw==0),        error('cx_H2_control:Tw','Disturbances at the boundary (L80-82).'),  end
if ~(PIE.Tu==0),        error('cx_H2_control:Tu','Inputs at the boundary (L84-86).'),  end
if ~(PIE.D11==0),       error('cx_H2_control:D11','Feedthrough D11 (L88-90).'),  end
if st.sosineq_on                                                            % L110-111, L172-175
    error('cx_H2_control:sosineq','sosineq_on: no container lpi_ineq (library gap).')
end

nx1 = PIE.A.dim(1,1);   nx2 = PIE.A.dim(2,1);                               % L123-124
Tm  = opvar2copvar(PIE.T);   Am  = opvar2copvar(PIE.A);
Bw  = opvar2copvar(PIE.B1);  B2m = opvar2copvar(PIE.B2);
Cz  = opvar2copvar(PIE.C1);  D12m = opvar2copvar(PIE.D12);
Iw  = opvar2copvar(mat2opvar(eye(size(PIE.B1,2)),PIE.B1.dim(:,2),PIE.vars,PIE.dom));    % L160

prog = lpiprogram(PIE.vars(:,1),PIE.vars(:,2),PIE.dom);                     % L122
[spx,dmx] = cx_space_list(Am,'out');       % [nx1;nx2] with zero-dim dropped
[prog,Pm] = cx_h2_lf(prog,dmx,spx,PIE.dom,st);                              % L142-148
% L150-151 as one multiplier: P + blkdiag(eppos*I_nx1, eppos2*I_nx2).
Pm = Pm + opvar2copvar(mat2opvar(blkdiag(st.eppos*eye(nx1),st.eppos2*eye(nx2)), ...
                                 [nx1;nx2],PIE.vars,PIE.dom));

% L153-154: Z maps the state (T's output spaces) to u (B2's input spaces).
[spu,dmu] = cx_space_list(B2m,'in');
[prog,Zm] = lpivar_cdopvar(prog,struct('out',dmu,'in',dmx), ...
                           struct('out',{spu},'in',{spx}),PIE.dom,st.ddZ);
% L156-157: W free and square on the z spaces; positive only through Dpos.
[spz,dmz] = cx_space_list(Cz,'out');
[prog,Wm] = lpivar_cdopvar(prog,dmz,spz,PIE.dom,st.ddZ);

Dneg = [-(gam*Iw),  Bw';                                                    % L162-163
        Bw,         (Tm*Pm)*Am' + (Am*Pm)*Tm' + (B2m*Zm)*Tm' + (Tm*Zm')*B2m'];
Dp12 = Cz*Pm + cx_on_registry(D12m,Zm)*Zm;                                  % L164
Dpos = [Wm,     Dp12;                                                       % L165-166
        Dp12',  Pm];

[prog,Nneg] = cx_h2_slack(prog,Dneg,PIE.dom,st);                            % L178-189
[prog,Npos] = cx_h2_slack(prog,Dpos,PIE.dom,st);
prog = lpi_eq_cdopvar(prog,Nneg + Dneg,'symmetric');                        % L190
prog = lpi_eq_cdopvar(prog,Npos - Dpos,'symmetric');                        % L191

prog = lpi_ineq(prog,gam - cx_h2_trace(Wm));                                % L167, L194
end
