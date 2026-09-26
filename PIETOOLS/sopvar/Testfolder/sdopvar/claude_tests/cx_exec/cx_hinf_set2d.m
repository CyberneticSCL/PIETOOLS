function S = cx_hinf_set2d(st,eppos_default)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% S = CX_HINF_SET2D(ST,EPPOS_DEFAULT) reads the 2-D LPI settings exactly as
% the 2-D H-infinity gain executives do (PIETOOLS_Hinf_gain_2D.m:72-128,
% _non_coercive:87-126, Hinf_gain_dual_2D:77-133):
%   - an lpisettings struct without is2D is replaced by its settings_2d,
%     keeping sos_opts;
%   - eppos defaults to EPPOS_DEFAULT when absent ([1e-4;1e-6;1e-6;1e-6]
%     in Hinf_gain_2D, [1e-4;0;0;1e-6] in Hinf_gain_dual_2D) and a scalar
%     is expanded to 4 entries;
%   - LF_deg/LF_opts, the enabled LF psatz terms, use_sosineq, eq_deg/eq_opts
%     and the enabled eq psatz terms (extract_psatz_deg/opts, local to each
%     executive, reproduced below).
% use_bisect/bisect_opts are not read: the caller fixes gamma and bisects
% outside, and no shipped settings file sets them.
%
% Initial coding MMP, 09/25/2026
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if ~isfield(st,'is2D') || ~st.is2D
    so = st.sos_opts;   st = st.settings_2d;    st.sos_opts = so;
end
S = struct();
if ~isfield(st,'eppos'),    S.eppos = eppos_default;
else,                       S.eppos = st.eppos;
end
if numel(S.eppos)==1,       S.eppos = S.eppos*ones(4,1);    end
S.LF_deg = st.LF_deg;       S.LF_opts = st.LF_opts;
S.LF_use_psatz  = st.LF_use_psatz;
S.LF_deg_psatz  = pick(st.LF_deg_psatz,S.LF_use_psatz);
S.LF_opts_psatz = pick(st.LF_opts_psatz,S.LF_use_psatz);
S.use_sosineq = st.use_sosineq;
if ~S.use_sosineq
    S.eq_opts = st.eq_opts;     S.eq_deg = st.eq_deg;
    S.eq_use_psatz  = st.eq_use_psatz;
    S.eq_deg_psatz  = pick(st.eq_deg_psatz,S.eq_use_psatz);
    S.eq_opts_psatz = pick(st.eq_opts_psatz,S.eq_use_psatz);
end
end

function out = pick(in,use)
% extract_psatz_deg / extract_psatz_opts (identical bodies in the executive).
out = {};
if all(use==0),     return,     end
if isa(in,'struct'),                out = repmat({in},1,numel(use));
elseif numel(in)==1,                out = repmat(in(1),1,numel(use));
elseif numel(in)==numel(use),       out = in;
elseif numel(in)>=max(use),         out = in(use);
else,   error('cx_hinf_set2d:psatz','A deg/opts entry is needed for each use_psatz term.')
end
end
