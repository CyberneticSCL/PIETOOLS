function prog = cx_H2_norm_2D_o(PIE,st,gam)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PROG = CX_H2_NORM_2D_O(PIE,ST,GAM) builds, with copvar/cdopvar only, the
% LPI of PIETOOLS_H2_norm_2D_o (2-D, observability gramian, coercive) at
% the FIXED value GAM, and returns it unsolved:
%
%   W >= 0 (+ eppos I),  A'*W*T + T'*W*A + C1'*C1 <= 0,
%   gam - trace((B1'*W*B1).R00) >= 0.
%
% GAM is the executive's decision variable gam, the SQUARED bound (stock
% options.h2 branch: gam = options.h2^2, L138; output sqrt(gam), L214).
% For structure comparison only in this suite: 2-D programs are not solved.
%
% Stock lines mirrored (executives/2D/PIETOOLS_H2_norm_2D_o.m):
%   L58-65    operators, Tw == 0 check;   L67-124 settings -> cx_h2_set2d
%   L127      lpiprogram
%   L149-157  W = poslpivar_2d(LF_deg,LF_opts) [+ psatz]    -> cx_h2_pos2d
%   L160-165  W += opvar2d(blkdiag(eppos(1)I,0,0,eppos(4)I)) -> plus
%   L172      Dop = (A'*W)*T + T'*(W*A) + C1'*C1
%   L187-198  slack at eq_deg [+ psatz], lpi_eq_2d 'symmetric'
%   L201-205  trace of (B1'*W*B1).R00, gam - trace >= 0      -> cx_h2_trace
%             (constant 0 when w has no R part: with a distributed w the
%             stock options.h2 branch then passes a double to lpi_ineq,
%             which errors, lpi_ineq.m:74-75; the objective form does not)
% NOT mirrored: L132-136 (gam as objective variable); L182-184
% (use_sosineq, error); non-L2[x,y] states (see cx_H2_norm_2D_c).
%
% Initial coding MMP, 09/25/2026
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if nargin<3,    error('cx_H2_norm_2D_o:gam','The container LPI is posed at fixed gamma.'),  end
if ~isa(PIE,'pie_struct'),  error('cx_H2_norm_2D_o:pie','PIE must be a pie_struct.'),  end
PIE = initialize(PIE);                                                      % L55
if ~(PIE.Tw==0)                                                             % L63-65
    error('cx_H2_norm_2D_o:Tw','H2 norm LPI cannot be solved with disturbances at the boundary.')
end
[X,eppos] = cx_h2_set2d(st);                                                % L67-124

Tm = opvar2d2copvar(PIE.T);     Am = opvar2d2copvar(PIE.A);
Bw = opvar2d2copvar(PIE.B1);    Cz = opvar2d2copvar(PIE.C1);

prog = lpiprogram(PIE.vars(:,1),PIE.vars(:,2),PIE.dom);                     % L127
[sp,dm] = cx_space_list(Tm,'out');
[prog,Wm] = cx_h2_pos2d(prog,sp,dm,PIE.dom,X.LF,X.LF_psatz);                % L149-157
if ~all(eppos==0)                                                           % L160-165
    np = PIE.T.dim(:,1);
    Ip = blkdiag(eppos(1)*eye(np(1)),zeros(np(2)),zeros(np(3)),eppos(4)*eye(np(4)));
    Wm = Wm + opvar2d2copvar(opvar2d(Ip,PIE.T.dim,PIE.dom,PIE.vars));
end

Dm = (Am'*Wm)*Tm + Tm'*(Wm*Am) + Cz'*Cz;                                    % L172

[spd,dmd] = cx_space_list(Dm,'out');
[prog,Nm] = cx_h2_pos2d(prog,spd,dmd,PIE.dom,X.eq,X.eq_psatz);              % L187-197
prog = lpi_eq_cdopvar(prog,Nm + Dm,'symmetric');                            % L198

prog = lpi_ineq(prog,gam - cx_h2_trace((Bw'*Wm)*Bw));                       % L201-205
end
