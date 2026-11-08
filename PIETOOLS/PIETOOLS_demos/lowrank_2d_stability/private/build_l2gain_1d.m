function [prog,H] = build_l2gain_1d(PIE,st,dual)                            % CC, 09/23/2026
% [PROG,H] = BUILD_L2GAIN_1D(PIE,ST,DUAL) assembles, without solving, the KYP
% LPI of PIETOOLS_Hinf_gain_sop, form 'Q' (DUAL false) or 'Qd' (DUAL true),
% with the stock negativity slack, for the L2 gain gam from w to z of the PIE
%
%       T x_t = A x + B1 w,      z = C1 x + D11 w:
%
%   'Q'   R >= 0,  T'*Q = R,
%         [-gam I, D11', B1'*Q; D11, -gam I, C1; Q'*B1, C1', A'*Q + Q'*A] + N = 0;
%   'Qd'  R >= 0,  T*Q = R,
%         [-gam I, D11, C1*Q; D11', -gam I, B1'; Q'*C1', B1, Q'*A' + A*Q] + N = 0;
%
% with R the positive storage operator, Q a free operator and N >= 0 the stock
% slack (degrees dd2, dd3, psatz flag override2). gam is a decision variable
% with gam >= 0 and the objective min gam, or the number ST.gamfix, which
% makes the LPI a feasibility test at that gain.
%
% INPUT
% - PIE:   'pie_struct' object of a 1-D PIE with 'copvar' (or 'opvar')
%          operators T, A, B1, C1, D11 and Tw = 0;
% - st:    'struct' of LPI settings as returned by lpisettings_sop: dd1,
%          dd12, options1, options12, override1 (storage), dd2, dd3,
%          options2, options3, override2 (slack); sosineq_on, if present,
%          must be 0; gamfix (optional) 'double', a fixed gain;
% - dual:  (optional) logical, true for form 'Qd' (default false).
% OUTPUT
% - prog:  'struct', the unsolved LPI program;
% - H:     'struct' with fields Top, Aop, Bwop, Czop, Dzwop ('copvar', T, A,
%          B1, C1, D11), Iw, Iz ('copvar' identities on w and z), Rop
%          ('cdopvar', R), Qop ('cdopvar', Q), Dop ('cdopvar', the KYP
%          operator), Deop ('cdopvar', the slack N), gam ('dpvar', or
%          'double' with ST.gamfix), dual, st, PIE.
%
% CHANGE LOG
% CC, 09/23/2026: initial coding, on opvar/poslpivar/lpivar/lpi_eq/lpi_ineq
%   (executives PIETOOLS_Hinf_gain and PIETOOLS_Hinf_gain_dual).
% CC, 09/24/2026: ST.gamfix, a numeric gain substituted for the decision
%   variable.
% CC, 11/08/2026: rewritten on the container LPI layer, the setup and the
%   'Q'/'Qd' assemble of PIETOOLS_Hinf_gain_sop with the stock slack branch;
%   opvar, poslpivar, lpivar, lpi_eq and lpi_ineq left the path (f83486c8,
%   433980d7). The opvar body of 09/23/2026 and 09/24/2026 is deleted.

% BEGIN CC, 11/08/2026: container build (replaces the opvar body, 09/23/2026
% and 09/24/2026)
if nargin<3 || isempty(dual), dual = false; end
% A disturbance at the boundary makes the executive pose form 'P' instead,
% a different LPI with a different residual: refuse rather than mis-adapt
if ~(PIE.Tw==0)
    error('build_l2gain_1d:Tw', ...
        ['PIE.Tw is nonzero (disturbance at the boundary): PIETOOLS_Hinf_gain_sop ' ...
         'poses form ''P'' there, a different LPI than this builder.']);
end
% The operator gate needs the equality route: sosineq_on poses no slack N
if isfield(st,'sosineq_on') && st.sosineq_on
    error('build_l2gain_1d:sosineq', ...
          'sosineq_on=1 has no slack operator; the operator gate needs the equality route.');
end
% Container PIE, completed settings, storage/slack builders, plant operators
PIE = exec_pie_sop(PIE);
st  = exec_sop_settings(st);
tl  = exec_tools_sop(PIE,st);    st.sop = tl.sop;
Tm  = PIE.T;     Am = PIE.A;
Bw  = PIE.B1;    Cz = PIE.C1;    Dzw = PIE.D11;
Iw  = tl.ident(PIE.B1,2);    Iz = tl.ident(PIE.C1,1);
% gam: a decision variable with gam >= 0 and min gam, or the fixed number
prog = lpiprogram(PIE.vars(:,1),PIE.vars(:,2),PIE.dom);
if isfield(st,'gamfix') && ~isempty(st.gamfix)
    gamv = st.gamfix;
else
    [prog,gamv] = lpidecvar(prog,'gam');
    prog = lpi_ineq_sop(prog,gamv);
    prog = lpisetobj(prog,gamv);
end
% Storage R >= 0 (no margin in the 'Q' forms), Q free, coupling T'Q = R (or TQ = R)
[prog,Rm] = tl.storage(prog,Tm,'out');
[sp,dm] = copvar_space_list(Tm,'out');
[prog,Qm] = lpivar_cdopvar(prog,dm,sp,tl.dom,get_lpivar_degs_sop(Rm));
if ~dual
    prog = lpi_eq_sop(prog,Tm'*Qm-Rm);
    Km = [-(gamv*Iw),   Dzw',       Bw'*Qm;
           Dzw,         -(gamv*Iz), Cz;
           Qm'*Bw,      Cz',        Am'*Qm + Qm'*Am];
else
    prog = lpi_eq_sop(prog,Tm*Qm-Rm);
    Km = [-(gamv*Iz),   Dzw,        Cz*Qm;
           Dzw',        -(gamv*Iw), Bw';
           Qm'*Cz',     Bw,         Qm'*Am' + Am*Qm];
end
% Km <= 0 through the stock slack: Km + N = 0 with N >= 0 of dd2/dd3/override2
[prog,Nm] = tl.slack(prog,Km);
prog = lpi_eq_sop(prog,Nm+Km,'symmetric');
% Handles for the operator gate (pielr_lpi l2g_resid)
H.Top = Tm;     H.Aop = Am;     H.Bwop = Bw;    H.Czop = Cz;    H.Dzwop = Dzw;
H.Iw = Iw;      H.Iz = Iz;      H.Rop = Rm;     H.Qop = Qm;
H.Dop = Km;     H.Deop = Nm;    H.gam = gamv;
H.dual = dual;  H.st = st;      H.PIE = PIE;
% END CC, 11/08/2026
end
