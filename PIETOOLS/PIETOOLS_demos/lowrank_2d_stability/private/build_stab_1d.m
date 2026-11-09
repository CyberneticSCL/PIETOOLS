function [prog,H] = build_stab_1d(PIE,st,dual)                              % CC, 09/23/2026
% [PROG,H] = BUILD_STAB_1D(PIE,ST,DUAL) assembles, without solving, the
% direct-form stability LPI of PIETOOLS_stability_sop, form 'P' (DUAL false)
% or 'Pd' (DUAL true), with the stock negativity slack:
%
%   'P'   P >= eppos I,   T'*P*A + A'*P*T + c*epneg*T'*P*T + N = 0,   N >= 0;
%   'Pd'  P >= eppos I,   T*P*A' + A*P*T' + c*epneg*T*P*T' + N = 0,   N >= 0;
%
% with c = 2 in every dimension, P the positive storage operator and N the
% slack. The slack is the stock one (degrees dd2, dd3 and
% the psatz flag override2), never the sized slack of LPI_INEQ_SOP, so the
% settings decide the program.
%
% INPUT
% - PIE:   'pie_struct' object of a 1-D or 2-D PIE with 'copvar' (or 'opvar',
%          'opvar2d') operators T and A;
% - st:    'struct' of LPI settings as returned by lpisettings_sop: dd1,
%          dd12, options1, options12, override1 (storage), dd2, dd3,
%          options2, options3, override2 (slack), eppos, eppos2, epneg, and in
%          2-D settings_2d; sosineq_on, if present, must be 0;
% - dual:  (optional) logical, true for form 'Pd' (default false).
% OUTPUT
% - prog:  'struct', the unsolved LPI program;
% - H:     'struct' with fields Top, Aop ('copvar', T and A), Pop ('cdopvar',
%          P including the margin), Dop ('cdopvar', the negativity
%          operator), Deop ('cdopvar', the slack N), eppos ('double' vector of
%          the margin of each space), epneg, epfac ('double', the factor c),
%          dual, st, PIE.
%
% CHANGE LOG
% CC, 09/23/2026: initial coding, on opvar/poslpivar/lpi_eq (executive
%   PIETOOLS_PDEstability and its dual).
% CC, 11/08/2026: rewritten on the container LPI layer, the body of the
%   'P'/'Pd' assemble of PIETOOLS_stability_sop with the stock slack branch;
%   opvar, poslpivar and lpi_eq left the path (f83486c8). The opvar body of
%   09/23/2026 is deleted; it could no longer run.
% MMP, 11/09/2026: c = 2 in 1-D too, as PIETOOLS_stability_sop poses the
%   decay term in every dimension since 11/09/2026 (maintainer: the
%   executives apply eppos and epneg the same way in every dimension).

% BEGIN CC, 11/08/2026: container build (replaces the opvar body, 09/23/2026)
if nargin<3 || isempty(dual), dual = false; end
% The operator gate needs the equality route: sosineq_on poses no slack N
if isfield(st,'sosineq_on') && st.sosineq_on
    error('build_stab_1d:sosineq', ...
          'sosineq_on=1 has no slack operator; the operator gate needs the equality route.');
end
% Container PIE, completed settings, and the storage/slack builders of this
% dimension, exactly as PIETOOLS_stability_sop prepares them
PIE = exec_pie_sop(PIE);
st  = exec_sop_settings(st);
tl  = exec_tools_sop(PIE,st);    st.sop = tl.sop;
eppos = tl.eppos;   epneg = tl.epneg;
if ~isempty(st.sop.eppos_stab),  eppos(:) = st.sop.eppos_stab;  end
% epfac = 1;  if tl.nv>=2,  epfac = 2;  end   % 2 epneg in 2-D, as the 2-D LPI % MMP, 11/09/2026 (was)
epfac = 2;          % 2 epneg in every dimension, as PIETOOLS_stability_sop % MMP, 11/09/2026
Tm = PIE.T;     Am = PIE.A;
% Storage P on the output side of T ('P') or its input side ('Pd'), plus margin
prog = lpiprogram(PIE.vars(:,1),PIE.vars(:,2),PIE.dom);
if ~dual
    [prog,Pm] = tl.storage(prog,Tm,'out');
    if ~all(eppos==0),  Pm = Pm + tl.eye(1,eppos);  end
    Dm = Tm'*Pm*Am + Am'*Pm*Tm;
    if epneg~=0,  Dm = Dm + epfac*epneg*(Tm'*Pm*Tm);  end
else
    [prog,Pm] = tl.storage(prog,Tm,'in');
    if ~all(eppos==0),  Pm = Pm + tl.eye(2,eppos);  end
    Dm = Tm*Pm*Am' + Am*Pm*Tm';
    if epneg~=0,  Dm = Dm + epfac*epneg*(Tm*Pm*Tm');  end
end
% Dm <= 0 through the stock slack: Dm + N = 0 with N >= 0 of dd2/dd3/override2
[prog,Nm] = tl.slack(prog,Dm);
prog = lpi_eq_sop(prog,Dm+Nm,'symmetric');
% Handles for the operator gate (pielr_lpi stab_resid)
H.Top = Tm;     H.Aop = Am;     H.Pop = Pm;
H.Dop = Dm;     H.Deop = Nm;
H.eppos = eppos;    H.epneg = epneg;    H.epfac = epfac;
H.dual = dual;  H.st = st;      H.PIE = PIE;
% END CC, 11/08/2026
end
