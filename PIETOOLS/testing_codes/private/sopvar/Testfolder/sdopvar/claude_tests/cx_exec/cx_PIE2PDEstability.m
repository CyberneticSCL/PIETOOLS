function [prog,ops] = cx_PIE2PDEstability(PIE,st)                           % MMP, 09/29/2026
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PROG = CX_PIE2PDESTABILITY(PIE,ST) builds, with copvar/cdopvar only, the
% LPI of executives/PIETOOLS_PIE2PDEstability.m (Q form, Cor. 7.1, the
% default 'stability' executive):
%
%   P = P1 (+P2) + eppos2*T'*T,     T'*Q - P = 0   (Q free),
%   D = A'*Q + Q'*A + epneg*P = -N,  N = N1 (+N2) >= 0,
%
% returned UNSOLVED. Mirrors lines 60-68 (2-D hand-off: the stock dispatches
% to the DIRECT-form PIETOOLS_stability_2D, so this goes to cx_stability_2D),
% 81-99 (settings), 105 (program), 112-120 (P), 123-125 (Q and T'Q = P),
% 137 (D), 148-163 (equality). Not mirrored: line 150 (sosineq_on=1, no
% container lpi_ineq, errors) and line 170 (no getsol for cdopvar).
%
% Container specifics:
% - Qdeg = get_lpivar_degs(Pop,Top) (line 123) has no container overload;
%   cx_stability_lpivar_degs reads the same maxima off the cdopvar P.
% - Q from lpivar_cdopvar with lpivar's [d1 d2 d3] convention, which
%   reproduces lpivar exactly in 1-D (lpivar_cdopvar.m:47-50, 429-430).
% - P lives on T's INPUT spaces (T'*T, T'*Q: in -> in).
%
% Initial coding MMP, 09/25/2026
% MMP, 09/29/2026: Second output OPS = struct of the program's operators
%                  T, A (copvar) and P, Q, D, N (cdopvar), so a solved
%                  program can be checked at operator level
%                  (lpi_programming_sopvar/tests/test_getsol_sop). Empty in
%                  2-D, whose program cx_stability_2D builds. One-output
%                  callers are unaffected. The signature was
%function prog = cx_PIE2PDEstability(PIE,st)                                % MMP, 09/29/2026 (was)
% MMP, 10/06/2026: P1(+P2) and N1(+N2) pairs -> library
%                  poslpivar_settings_sop ('lf', 'slack'),
%                  cx_stability_lpivar_degs -> get_lpivar_degs_sop,
%                  cx_space_list -> copvar_space_list, so the transcription
%                  runs on the library translators (lpi_programming_sopvar);
%                  programs unchanged (126-program bit-identity check).
%                  'cx_stability_lpivar_degs' above now refers to the
%                  library routine get_lpivar_degs_sop.
% MMP, 11/07/2026: T and A are used directly: PIEs hold 'copvar'
%                  operators, on which opvar2copvar errors.
% MMP, 11/08/2026: Default settings from lpisettings_sop('heavy') (was the
%                  stock preset settings_PIETOOLS_heavy, whose folder leaves
%                  the path): the fields read below (dd1, dd12, options1,
%                  options12, override1, sosineq_on) are the same; epneg,
%                  eppos2 and sos_opts are set below as before.
% MMP, 11/09/2026: eppos2 read as EXEC_SOP_SETTINGS sets it: the L2 margin
%                  when ST has one, else the one margin eppos. The
%                  settings_sopvar tiers set eppos only since 11/09/2026, so
%                  an ST from lpisettings_sop stopped at st.eppos2
%                  (e2e_stability_1d, test_getsol_sop part ii). The value is
%                  unchanged for every tier: eppos2 was eppos, 1e-4. The
%                  decay term epneg*P of the stock form is kept (the
%                  container executive poses 2*epneg*P since 11/09/2026).
% MMP, 11/12/2026: The settings are native (phase 2 of the knob reduction,
%                  maintainer 11/12/2026): the stock-format degree fields,
%                  read here and never used, are no longer read; P comes
%                  from lpi.pos and N from sop.stock (POSLPIVAR_SETTINGS_SOP),
%                  which CX_SETTINGS and the default ST now carry. The same
%                  program for every tier.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if ~isa(PIE,'pie_struct')
    error('The PIE should be a ''pie_struct''.')
end
PIE = initialize(PIE);
ops = struct();     % 2-D: operators not returned (built in cx_stability_2D) % MMP, 09/29/2026
if PIE.dim==2                                       % stock lines 60-68
    if nargin<2,    prog = cx_stability_2D(PIE);
    else,           prog = cx_stability_2D(PIE,st);
    end
    return
end
if nargin<2                                         % stock lines 73-79
%   st = settings_PIETOOLS_heavy;                                           % MMP, 11/08/2026 (was)
    st = lpisettings_sop('heavy');                                          % MMP, 11/08/2026
    [~,st.sop.stock] = lpisettings_stock_sop('heavy');  % the stock slack   % MMP, 11/12/2026
    st.sos_opts.simplify = 1;
    st.eppos = 1e-4;    st.eppos2 = 1e-6;   st.epneg = 0;
end
% dd1 = st.dd1;   dd12 = st.dd12;     options1 = st.options1;     options12 = st.options12; % MMP, 11/12/2026 (was)
% override1 = st.override1;   epneg = st.epneg;   eppos2 = st.eppos2;      % MMP, 11/09/2026 (was)
% override1 = st.override1;   epneg = st.epneg;                             % MMP, 11/12/2026 (was); MMP, 11/09/2026
% The degree fields above were read and never used (P and N come from       % MMP, 11/12/2026
% POSLPIVAR_SETTINGS_SOP, which reads lpi.pos and sop.stock since           % MMP, 11/12/2026
% 11/12/2026, when they left the settings)                                  % MMP, 11/12/2026
epneg = st.epneg;                                                           % MMP, 11/12/2026
% The L2 margin, else the one margin of the settings_sopvar tiers           % MMP, 11/09/2026
if isfield(st,'eppos2') && ~isempty(st.eppos2),  eppos2 = st.eppos2;  else,  eppos2 = st.eppos;  end % MMP, 11/09/2026
% if st.sosineq_on                                    % stock lines 148-150 % MMP, 11/12/2026 (was)
%     error('cx_stability:gap',['sosineq_on=1 needs lpi_ineq, which has no container '... % MMP, 11/12/2026 (was)
%           'counterpart; only the equality path (the default) is transcribed.']) % MMP, 11/12/2026 (was)
% end                                                                       % MMP, 11/12/2026 (was)
% override2 = st.override2;   options2 = st.options2;     options3 = st.options3; % MMP, 11/12/2026 (was)
% dd2 = st.dd2;   dd3 = st.dd3;                                             % MMP, 11/12/2026 (was)
% (sosineq_on left the settings: the transcription poses the equality path) % MMP, 11/12/2026

% Tm = opvar2copvar(PIE.T);   Am = opvar2copvar(PIE.A);                     % MMP, 11/07/2026 (was)
Tm = PIE.T;     Am = PIE.A;                 % 'copvar' since initialize     % MMP, 11/07/2026
prog = lpiprogram(PIE.vars(:,1),PIE.vars(:,2),PIE.dom);        % line 105

% P = P1 (+P2) + eppos2*T'*T, lines 112-120.
% [prog,Pm] = cx_stability_pos(prog,Tm,'in',PIE.dom,dd1,options1);          % MMP, 10/06/2026 (was)
% if override1~=1                                                           % MMP, 10/06/2026 (was)
%     [prog,P2m] = cx_stability_pos(prog,Tm,'in',PIE.dom,dd12,options12);   % MMP, 10/06/2026 (was)
%     Pm = Pm + P2m;                                                        % MMP, 10/06/2026 (was)
% end                                                                       % MMP, 10/06/2026 (was)
[prog,Pm] = poslpivar_settings_sop(prog,Tm,st,'lf','in',PIE.dom);           % MMP, 10/06/2026
Pm = Pm + eppos2*Tm'*Tm;                            % stock association (eppos2*T')*T

% Q free with T'*Q = P, lines 123-125. Q has T's shape (Top.dim).
% Qdeg = cx_stability_lpivar_degs(Pm);                                      % MMP, 10/06/2026 (was)
% [spO,dmO] = cx_space_list(Tm,'out');    [spI,dmI] = cx_space_list(Tm,'in'); % MMP, 10/06/2026 (was)
Qdeg = get_lpivar_degs_sop(Pm);                                             % MMP, 10/06/2026
[spO,dmO] = copvar_space_list(Tm,'out');    [spI,dmI] = copvar_space_list(Tm,'in'); % MMP, 10/06/2026
[prog,Qm] = lpivar_cdopvar(prog,struct('out',dmO,'in',dmI), ...
                           struct('out',{spO},'in',{spI}),PIE.dom,Qdeg);
prog = lpi_eq_cdopvar(prog,Tm'*Qm - Pm);            % NOT symmetric (line 125)

% D, line 137.
Dm = Am'*Qm + Qm'*Am;
if epneg~=0
    Dm = Dm + epneg*Pm;
end

% D = -N, lines 154-162.
% [prog,Nm] = cx_stability_pos(prog,Dm,'out',PIE.dom,dd2,options2);         % MMP, 10/06/2026 (was)
% if override2~=1                                                           % MMP, 10/06/2026 (was)
%     [prog,N2m] = cx_stability_pos(prog,Dm,'out',PIE.dom,dd3,options3);    % MMP, 10/06/2026 (was)
%     Nm = Nm + N2m;                                                        % MMP, 10/06/2026 (was)
% end                                                                       % MMP, 10/06/2026 (was)
[prog,Nm] = poslpivar_settings_sop(prog,Dm,st,'slack','out',PIE.dom);       % MMP, 10/06/2026
prog = lpi_eq_cdopvar(prog,Dm+Nm,'symmetric');
% Operators the two lpi_eq relations constrain, for operator-level checks   % MMP, 09/29/2026
ops = struct('T',Tm,'A',Am,'P',Pm,'Q',Qm,'D',Dm,'N',Nm);                    % MMP, 09/29/2026
end
