function S = settings_2d_sop(st)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% S = SETTINGS_2D_SOP(ST) the 2-D settings the stock 2-D executives read,
% from an lpisettings struct ST (its field settings_2d) or from a 2-D
% settings struct itself (is2D set), with the defaults the stock executives
% apply filled in:
%   eppos         4 x 1 (a scalar is spread), default [1e-4;1e-6;1e-6;1e-6];
%   epneg         default 0;
%   LF_deg_psatz, LF_opts_psatz, eq_deg_psatz, eq_opts_psatz
%                 cells with one entry per nonzero use_psatz term, as the
%                 stock extract_psatz_deg / extract_psatz_opts return them
%                 ({} when no term is used);
%   sos_opts      ST's.
% The other fields (LF_deg, LF_opts, LF_use_psatz, eq_deg, eq_opts,
% eq_use_psatz, use_sosineq, Zop_deg, Zop_opts, ...) are passed through.
%
% Initial coding MMP, 10/08/2026
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if isfield(st,'is2D') && st.is2D
    S = st;
elseif isfield(st,'settings_2d')
    S = st.settings_2d;
    if isfield(st,'sos_opts'),  S.sos_opts = st.sos_opts;   end
else
    error('settings_2d_sop:st','ST should be an lpisettings struct (field settings_2d) or a 2-D settings struct.');
end
S.eppos_default = ~isfield(S,'eppos') || isempty(S.eppos);   % the dual gain executive defaults differently
if S.eppos_default,     S.eppos = [1e-4;1e-6;1e-6;1e-6];    end
if numel(S.eppos)==1,   S.eppos = S.eppos*ones(4,1);    end
S.eppos = S.eppos(:);
if ~isfield(S,'epneg') || isempty(S.epneg),     S.epneg = 0;    end
S.LF_deg_psatz  = pick(S.LF_deg_psatz, S.LF_use_psatz);
S.LF_opts_psatz = pick(S.LF_opts_psatz,S.LF_use_psatz);
if isfield(S,'use_sosineq') && ~S.use_sosineq
    S.eq_deg_psatz  = pick(S.eq_deg_psatz, S.eq_use_psatz);
    S.eq_opts_psatz = pick(S.eq_opts_psatz,S.eq_use_psatz);
end
end


function out = pick(in,use)
% one entry per use_psatz term, as the stock extract_psatz_deg/opts
out = {};
if all(use==0),     return,     end
if isa(in,'struct'),                out = repmat({in},1,numel(use));
elseif numel(in)==1,                out = repmat(in(1),1,numel(use));
elseif numel(in)==numel(use),       out = in;
elseif numel(in)>=max(use),         out = in(use);
else,   error('settings_2d_sop:psatz','A deg/opts entry is needed for each use_psatz term.');
end
end
