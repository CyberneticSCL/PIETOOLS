function [X,eppos] = cx_h2_set2d(st)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% [X,EPPOS] = CX_H2_SET2D(ST) reads the 2-D coercive H2 executives'
% settings as they do (PIETOOLS_H2_norm_2D_c.m:67-122, identical in
% PIETOOLS_H2_norm_2D_o.m) and translates the positive-operator parts for
% poscopvar through settings2possopvar:
%
%   X.LF, X.LF_psatz   the gramian W and its enabled psatz terms
%   X.eq, X.eq_psatz   the negativity slack and its enabled psatz terms
%   EPPOS              settings_2d.eppos, default [1e-4;1e-6;1e-6;1e-6]
%
% Mirrored: settings = settings.settings_2d (L69); eppos default (L85-91);
% the psatz slack options OR-ed with eq_opts before use (L190-191), applied
% here to the settings so the translation sees what poslpivar_2d sees.
% Refused, not approximated: use_sosineq (L180-182, no container lpi_ineq),
% and anything settings2possopvar reports as lossy (psatz = 2 terms, mixed
% sep regimes): dropping them would build a different cone.
%
% Initial coding MMP, 09/25/2026
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

s = st.settings_2d;                                                         % L69
if isfield(s,'eppos'),  eppos = s.eppos;    else,   eppos = [1e-4;1e-6;1e-6;1e-6];  end
if s.use_sosineq
    error('cx_h2_set2d:sosineq','use_sosineq: no container lpi_ineq (library gap).')
end
for j = 1:numel(s.eq_use_psatz)                                             % L190-191
    if s.eq_use_psatz(j)~=0
        s.eq_opts_psatz{j}.exclude = s.eq_opts_psatz{j}.exclude | s.eq_opts.exclude;
        s.eq_opts_psatz{j}.sep = s.eq_opts_psatz{j}.sep | s.eq_opts.sep;
    end
end
X = settings2possopvar(s);
if ~isempty(X.report.lossy)
    error('cx_h2_set2d:lossy','settings not representable: %s',strjoin(X.report.lossy,' | '))
end
end
