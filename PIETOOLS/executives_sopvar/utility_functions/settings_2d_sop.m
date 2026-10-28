function S = settings_2d_sop(st)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% S = SETTINGS_2D_SOP(ST) returns the 2-D LPI settings held by ST, with
% the margins eppos and epneg set to their defaults when absent, and the
% degrees and options of the Positivstellensatz terms laid out as one cell
% entry per entry of the corresponding use_psatz field.
%
% INPUT
% - st:   'struct', either LPI settings with the field settings_2d and,
%         optionally, sos_opts (as returned by lpisettings or
%         lpisettings_sop), or 2-D settings with is2D true (as returned by
%         settings_PIETOOLS_light_2D); the 2-D settings have the fields
%         LF_deg_psatz, LF_opts_psatz and LF_use_psatz and, when
%         use_sosineq is 0, eq_deg_psatz, eq_opts_psatz and eq_use_psatz.
%
% OUTPUT
% - S:    'struct', the 2-D settings with every field kept, except
%           eppos          4x1 double, one margin for each of the spaces R,
%                          L2[s1], L2[s2] and L2[s1,s2] (a scalar is
%                          repeated; default [1e-4;1e-6;1e-6;1e-6]);
%           eppos_default  logical, true when the settings carry no eppos;
%           epneg          double (default 0);
%           LF_deg_psatz, LF_opts_psatz   cell, one degree / options
%                          'struct' per entry of LF_use_psatz ({} when every
%                          entry is 0): a 'struct' or a one-entry cell is
%                          repeated, a cell with numel(LF_use_psatz)
%                          entries is kept, and otherwise a cell with at
%                          least max(LF_use_psatz) entries is indexed by
%                          LF_use_psatz;
%           eq_deg_psatz, eq_opts_psatz   the same for eq_use_psatz, when
%                          use_sosineq is 0;
%           sos_opts       ST.sos_opts, when ST holds settings_2d and
%                          sos_opts.
%
% NOTES
% ST with neither is2D true nor the field settings_2d is an error, as is a
% psatz cell that fits none of the cases above.
%
% For support, contact M. Peet, Arizona State University at mpeet@asu.edu

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PIETOOLS - settings_2d_sop
%
% Copyright (C) 2026 PIETOOLS Team
%
% This program is free software; you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation; either version 2 of the License, or
% (at your option) any later version.
%
% This program is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
%
% You should have received a copy of the GNU General Public License
% along with this program; if not, write to the Free Software
% Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% If you modify this code, document all changes carefully and include date
% authorship, and a brief description of modifications
%
% Initial coding MMP, 10/08/2026
% MMP, 10/28/2026: Header to the header protocol (what the function does,
%                inputs, outputs), the licence block added, step comments
%                in the body and the header of the local function; code
%                unchanged. The deleted header text was the comparison with
%                the stock 2-D executives.

% The 2-D settings: ST itself, or ST.settings_2d with the solver options
% of ST
if isfield(st,'is2D') && st.is2D
    S = st;
elseif isfield(st,'settings_2d')
    S = st.settings_2d;
    if isfield(st,'sos_opts'),  S.sos_opts = st.sos_opts;   end
else
    error('settings_2d_sop:st','ST should be an lpisettings struct (field settings_2d) or a 2-D settings struct.');
end
% Margins: the default when absent, a scalar repeated on the four spaces, a
% column; the decay rate 0 when absent
S.eppos_default = ~isfield(S,'eppos') || isempty(S.eppos);   % the dual gain executive defaults differently
if S.eppos_default,     S.eppos = [1e-4;1e-6;1e-6;1e-6];    end
if numel(S.eppos)==1,   S.eppos = S.eppos*ones(4,1);    end
S.eppos = S.eppos(:);
if ~isfield(S,'epneg') || isempty(S.epneg),     S.epneg = 0;    end
% One psatz degree / options entry per use_psatz entry: the terms of the
% storage, and those of the equality slack when sosineq is not used
S.LF_deg_psatz  = pick(S.LF_deg_psatz, S.LF_use_psatz);
S.LF_opts_psatz = pick(S.LF_opts_psatz,S.LF_use_psatz);
if isfield(S,'use_sosineq') && ~S.use_sosineq
    S.eq_deg_psatz  = pick(S.eq_deg_psatz, S.eq_use_psatz);
    S.eq_opts_psatz = pick(S.eq_opts_psatz,S.eq_use_psatz);
end
end


function out = pick(in,use)
% OUT = PICK(IN,USE) returns one entry of IN for each entry of USE, or {}
% when every entry of USE is 0.
% INPUT  in:  'struct', or cell of degree or options entries;
%        use: double array, a use_psatz field.
% OUTPUT out: cell: IN repeated for a 'struct' or a one-entry cell, IN
%             itself with numel(USE) entries, else IN(USE) (an error when
%             IN has fewer than max(USE) entries).
% The cases of the stock extract_psatz_deg / extract_psatz_opts
out = {};
if all(use==0),     return,     end
if isa(in,'struct'),                out = repmat({in},1,numel(use));
elseif numel(in)==1,                out = repmat(in(1),1,numel(use));
elseif numel(in)==numel(use),       out = in;
elseif numel(in)>=max(use),         out = in(use);
else,   error('settings_2d_sop:psatz','A deg/opts entry is needed for each use_psatz term.');
end
end
