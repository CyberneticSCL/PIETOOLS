function st = exec_sop_settings(st)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ST = EXEC_SOP_SETTINGS(ST) completes the LPI settings ST with the options
% of the container executives: each field of ST.sop listed below that is
% missing or empty is set to its default, and ST.sos_opts.solver, when
% missing or empty, is set to 'mosek' if mosekopt is on the MATLAB path and
% to 'sedumi' otherwise. A scalar ST.eppos and ST.epneg set the margin and
% decay-rate fields of every dimension that ST leaves unset.
%
% INPUT
% - st:   'struct' of LPI settings, as returned by lpisettings_sop, with
%         the optional fields
%           eppos     double, the margin of the positive operator on every
%                     space in every dimension: a scalar sets eppos2 (1-D,
%                     L2), settings_2d.eppos (2-D) and settings_nd.eppos (3
%                     or more variables) where unset;
%           epneg     double, the decay rate in every dimension: sets
%                     settings_2d.epneg where unset (3 or more variables
%                     read epneg itself);
%           sos_opts  'struct' of solver options for lpisolve;
%           sop       'struct' of options of the container executives, each
%                     field optional:
%             slack       char, 'like' (default) or 'stock': each
%                         negativity constraint -D >= 0 is imposed by
%                         LPI_INEQ_SOP with a slack sized from the positive
%                         operator of the LPI ('like'), or as D + N == 0
%                         with N the positive slack declared from the
%                         degrees of the settings ('stock');
%             dw_Q        double (default 1), the raise of the weight
%                         degree of the 'like' slack (opts.dw of
%                         LPI_INEQ_SOP) in the non-coercive forms;
%             dw_P        double (default 2), the same in the other LPIs;
%             max_raise   double (default 2), the largest number of raises
%                         of the degree loop LPI_SOLVE_LOOP_SOP, which
%                         rebuilds the LPI with the slack weight raised by
%                         one while the solver returns 'unknown' or
%                         'infeasible';
%             margin      double eps >= 0 (default 0), certify
%                         -D + eps*e_D >= 0 in place of -D >= 0, with
%                         e_D = A'*A and A the lift of the plain term of
%                         the slack (opts.margin of LPI_INEQ_SOP);
%             witness     logical (default true), compute for a 1-D PIE the
%                         gain, decay rate or H2 norm of its Legendre
%                         discretization (PIE_WITNESS_SOP);
%             N_cheb      double (default 24), the degree of that
%                         discretization;
%             nfreq       double (default 240), its number of frequency
%                         samples;
%             keepdual    logical (default true), record the equality rows
%                         of the program (prog.sopeq) and assemble the dual
%                         kernel of the negativity constraint
%                         (LPIGETDUAL_SOP);
%             verbose     logical (default true), print each pass and the
%                         result;
%             dw_feas     double, a scalar or one entry per dimension
%                         ([1-D, 2-D, ...], the last entry serving every
%                         higher dimension), the weight raise of the
%                         stability and well-posedness LPIs in place of dw_P
%                         and dw_Q (default []: dw_P and dw_Q);
%             max_raise_feas  double, a scalar or one entry per dimension,
%                         as dw_feas, max_raise of the stability and
%                         well-posedness LPIs (default []: max_raise);
%             cap_delta   double, a scalar or one entry per dimension, as
%                         dw_feas: when set, the lift basis of every slack
%                         LPI_INEQ_SOP declares is bounded per direction,
%                         theta_k + s_k <= Jmax_k(T) + cap_delta (Jmax of
%                         GET_LIFT_DEGS on the PIE's T; default []: no
%                         bound, every program unchanged);
%             eppos_stab  double, a scalar or one entry per space, the
%                         margin of the positive operator of the stability
%                         and well-posedness LPIs in place of eppos, eppos2
%                         (1-D) and settings_2d.eppos (2-D) (default []:
%                         those fields).
%
% OUTPUT
% - st:   'struct', ST with the fields of ST.sop above set (its other
%         fields kept), ST.sos_opts.solver set, and the margin and
%         decay-rate fields of each dimension set from eppos and epneg
%         where unset.
%
% NOTES
% ST other than a 'struct', or ST.sop.slack other than 'like' and 'stock',
% is an error.
%
% For support, contact M. Peet, Arizona State University at mpeet@asu.edu

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PIETOOLS - exec_sop_settings
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
% MMP, 10/16/2026: Moved from lpi_programming_sopvar/ to
%                executives_sopvar/utility_functions/: read by the container
%                executives only.
% MMP, 10/28/2026: Fields dw_feas, max_raise_feas, eppos_stab, read by the
%                feasibility executives (PIETOOLS_stability_sop,
%                PIETOOLS_well_posedness_sop) in place of dw_P / dw_Q,
%                max_raise and eppos, and set by the tiers of settings_sopvar
%                (LPISETTINGS_SOP): feasibility LPIs certify at a lower
%                weight than the metric LPIs that share dw_P. Their default
%                [] keeps every program of the existing settings unchanged.
%                Header to the header protocol, the licence block added, step
%                comments in the body; code unchanged.
% MMP, 10/29/2026: max_raise_feas documented as one entry per dimension, as
%                dw_feas (the last entry serves every higher dimension, entry
%                min(nv,end) in PIETOOLS_stability_sop). Field cap_delta
%                (default [], no change): per-direction bounds
%                theta_k + s_k <= Jmax_k(T) + cap_delta on the slack bases
%                (opts.cap of LPI_INEQ_SOP, EXEC_TOOLS_SOP capfor), the bound
%                read from T that made the 3-D heat benchmark solvable.
% MMP, 11/08/2026: Comment pass (maintainer, 11/08/2026): log condensed, help
%                names lpisettings_sop; code unchanged.
% MMP, 11/09/2026: One margin and one decay rate for every dimension
%                (maintainer: "I see no reason to treat the eppos and epneg
%                options differently for 1D,2D and ND"): a scalar eppos and
%                epneg of ST set eppos2 (1-D), settings_2d.eppos and epneg
%                (2-D) and settings_nd.eppos (3 or more variables) where ST
%                leaves them unset, so the tiers of settings_sopvar set each
%                once; a per-dimension field that is set is kept. Here and
%                not in EXEC_TOOLS_SOP, which reads those fields: every
%                caller of EXEC_TOOLS_SOP calls this function first.

% Defaults of the fields of ST.sop
% d = struct('slack','like','dw_Q',1,'dw_P',2,'max_raise',2,'margin',0,'witness',true, ...
%            'N_cheb',24,'nfreq',240,'keepdual',true,'verbose',true);       % MMP, 10/28/2026 (was)
% d = struct('slack','like','dw_Q',1,'dw_P',2,'max_raise',2,'margin',0,'witness',true, ...
%            'N_cheb',24,'nfreq',240,'keepdual',true,'verbose',true, ...
%            'dw_feas',[],'max_raise_feas',[],'eppos_stab',[]);             % MMP, 10/29/2026 (was)
d = struct('slack','like','dw_Q',1,'dw_P',2,'max_raise',2,'margin',0,'witness',true, ...
           'N_cheb',24,'nfreq',240,'keepdual',true,'verbose',true, ...
           'dw_feas',[],'max_raise_feas',[],'eppos_stab',[],'cap_delta',[]); % MMP, 10/29/2026
% Check ST; set each missing or empty field of ST.sop to its default
if ~isstruct(st),   error('exec_sop_settings:st','ST should be a settings struct (lpisettings).');   end
if ~isfield(st,'sop') || isempty(st.sop),   st.sop = struct();     end
fn = fieldnames(d);
for i = 1:numel(fn)
    if ~isfield(st.sop,fn{i}) || isempty(st.sop.(fn{i})),   st.sop.(fn{i}) = d.(fn{i});    end
end
% Check the slack mode
if ~any(strcmp(st.sop.slack,{'like','stock'}))
    error('exec_sop_settings:slack','sop.slack should be ''like'' or ''stock''.');
end
% Solver: MOSEK when mosekopt is on the path, else SeDuMi
if ~isfield(st,'sos_opts') || isempty(st.sos_opts),    st.sos_opts = struct();    end
if ~isfield(st.sos_opts,'solver') || isempty(st.sos_opts.solver)
    if ~isempty(which('mosekopt')),     st.sos_opts.solver = 'mosek';
    else,                               st.sos_opts.solver = 'sedumi';
    end
end
% if st.sop.keepdual,     st.sos_opts.simplify = 0;   end                   % MMP, 10/08/2026 (was)
% BEGIN MMP, 11/09/2026: one margin and one decay rate for every dimension.
% EXEC_TOOLS_SOP reads the margin from eppos and eppos2 in 1-D, from
% settings_2d.eppos in 2-D (SETTINGS_2D_SOP) and settings_nd.eppos in 3 or
% more variables, and the decay rate from epneg, settings_2d.epneg and
% settings_nd.epneg (default epneg); a scalar eppos and epneg fill those
% left unset. A non-scalar eppos is the 1-D R^n margin of the stock
% settings and fills nothing.
if isfield(st,'eppos') && isscalar(st.eppos)
    if ~isfield(st,'eppos2') || isempty(st.eppos2),     st.eppos2 = st.eppos;   end
    if isfield(st,'settings_2d') && isstruct(st.settings_2d) && ...
            (~isfield(st.settings_2d,'eppos') || isempty(st.settings_2d.eppos))
        st.settings_2d.eppos = st.eppos;
    end
    if ~isfield(st,'settings_nd') || ~isstruct(st.settings_nd),  st.settings_nd = struct();  end
    if ~isfield(st.settings_nd,'eppos') || isempty(st.settings_nd.eppos)
        st.settings_nd.eppos = st.eppos;
    end
end
if isfield(st,'epneg') && isscalar(st.epneg) && isfield(st,'settings_2d') && ...
        isstruct(st.settings_2d) && (~isfield(st.settings_2d,'epneg') || isempty(st.settings_2d.epneg))
    st.settings_2d.epneg = st.epneg;
end
% END MMP, 11/09/2026
end
