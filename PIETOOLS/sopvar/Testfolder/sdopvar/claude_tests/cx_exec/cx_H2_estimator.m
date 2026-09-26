function prog = cx_H2_estimator(PIE,st,gam)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PROG = CX_H2_ESTIMATOR(PIE,ST,GAM) builds, with copvar/cdopvar only, the
% LPI of PIETOOLS_H2_estimator (1-D H2-optimal estimator synthesis) at the
% FIXED bound GAM, and returns it unsolved:
%
%   P >= 0 (+ eppos margin), Z free (y -> state), W free (w -> w),
%   Dneg = [-gam*Iz, C1; C1', T'*P*A + A'*P*T + T'*Z*C2 + C2'*Z'*T] <= 0,
%   Dpos = [W, -(B1'*P + D21'*Z'); -(B1'*P + D21'*Z')', P] >= 0,
%   gam - trace(W.P) >= 0.
%
% GAM is the executive's decision variable gam, the error H2 bound (stock
% output gam, L203) - the quantity cx_stock_pinned pins.
%
% Stock lines mirrored (executives/PIETOOLS_H2_estimator.m):
%   L66-87    operators; 2-D, Tw, D11 checks
%   L96-115   settings read exactly as the stock
%   L119-121  lpiprogram, nx1, nx2
%   L134-140  P = poslpivar(T.dim(:,1),dd1,options1) [+ dd12] -> cx_h2_lf
%   L142-143  P.P += eppos*I, P.R.R0 += eppos2*I  -> plus of the multiplier
%   L145-146  Z = lpivar([T.dim(:,1),C2.dim(:,1)],ddZ)  -> lpivar_cdopvar
%   L148-149  W = lpivar([B1.dim(:,2),B1.dim(:,2)],ddZ) -> lpivar_cdopvar on
%             the w spaces (an L2 operator when w is distributed)
%   L154-160  Iz, Dneg, D12, Dpos; left-to-right products; D21' restated on
%             the registry of Z' (cx_on_registry)
%   L174-187  slacks at dd2 [+ dd3 psatz], lpi_eq 'symmetric' -> cx_h2_slack
%   L161,190  gam - trace(W.P) >= 0 via the stock scalar lpi_ineq
%             -> cx_h2_trace (constant 0 when w has no R^n part: the stock
%             trace of a 0x0 dpvar is likewise a 1x1 dpvar)
% NOT mirrored: L126-127, L191, L195 (gam decision variable, gam >= 0,
% objective); L168-171 (sosineq_on, error); L200-208 (getsol/getObserver).
%
% Initial coding MMP, 09/25/2026
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if nargin<3,    error('cx_H2_estimator:gam','The container LPI is posed at fixed gamma.'),  end
if ~isa(PIE,'pie_struct'),  error('cx_H2_estimator:pie','PIE must be a pie_struct.'),  end
PIE = initialize(PIE);
if PIE.dim==2,          error('cx_H2_estimator:dim','2-D PIEs are not supported (L77-79).'),  end
if ~(PIE.Tw==0),        error('cx_H2_estimator:Tw','Disturbances at the boundary (L81-83).'),  end
if ~(PIE.D11==0),       error('cx_H2_estimator:D11','Feedthrough D11 (L85-87).'),  end
if st.sosineq_on                                                            % L107-108, L168-171
    error('cx_H2_estimator:sosineq','sosineq_on: no container lpi_ineq (library gap).')
end

nx1 = PIE.A.dim(1,1);   nx2 = PIE.A.dim(2,1);                               % L120-121
Tm  = opvar2copvar(PIE.T);   Am  = opvar2copvar(PIE.A);
Bw  = opvar2copvar(PIE.B1);  C2m = opvar2copvar(PIE.C2);
Cz  = opvar2copvar(PIE.C1);  D21m = opvar2copvar(PIE.D21);
Iz  = opvar2copvar(mat2opvar(eye(size(PIE.C1,1)),PIE.C1.dim(:,1),PIE.vars,PIE.dom));    % L154

prog = lpiprogram(PIE.vars(:,1),PIE.vars(:,2),PIE.dom);                     % L119
[spx,dmx] = cx_space_list(Tm,'out');       % T.dim(:,1) with zero-dim dropped
[prog,Pm] = cx_h2_lf(prog,dmx,spx,PIE.dom,st);                              % L134-140
Pm = Pm + opvar2copvar(mat2opvar(blkdiag(st.eppos*eye(nx1),st.eppos2*eye(nx2)), ...
                                 [nx1;nx2],PIE.vars,PIE.dom));              % L142-143

% L145-146: Z maps y (C2's output spaces) to the state.
[spy,dmy] = cx_space_list(C2m,'out');
[prog,Zm] = lpivar_cdopvar(prog,struct('out',dmx,'in',dmy), ...
                           struct('out',{spx},'in',{spy}),PIE.dom,st.ddZ);
% L148-149: W free and square on the w spaces.
[spw,dmw] = cx_space_list(Bw,'in');
[prog,Wm] = lpivar_cdopvar(prog,dmw,spw,PIE.dom,st.ddZ);

Dneg = [-(gam*Iz),  Cz;                                                     % L156-157
        Cz',        (Tm'*Pm)*Am + (Am'*Pm)*Tm + (Tm'*Zm)*C2m + (C2m'*Zm')*Tm];
D12 = Bw'*Pm + cx_on_registry(D21m',Zm')*Zm';                               % L158
Dpos = [Wm,       -D12;                                                     % L159-160
        -D12',    Pm];

[prog,Nneg] = cx_h2_slack(prog,Dneg,PIE.dom,st);                            % L174-185
[prog,Npos] = cx_h2_slack(prog,Dpos,PIE.dom,st);
prog = lpi_eq_cdopvar(prog,Nneg + Dneg,'symmetric');                        % L186
prog = lpi_eq_cdopvar(prog,Npos - Dpos,'symmetric');                        % L187

prog = lpi_ineq(prog,gam - cx_h2_trace(Wm));                                % L161, L190
end
