function st = exec_sop_settings(st)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ST = EXEC_SOP_SETTINGS(ST) completes the LPI settings ST with the options
% of the executives: each field of ST.sop listed below that is missing or
% empty is set to its default, and ST.sos_opts.solver, when missing or
% empty, is set to 'mosek' if mosekopt is on the MATLAB path and to
% 'sedumi' otherwise. A scalar ST.eppos and ST.epneg set the margin and
% decay-rate fields of every dimension that ST leaves unset.
%
% INPUT
% - st:   'struct' of LPI settings, with the optional fields
%           eppos     double, the margin of the positive operator on every
%                     space in every dimension: a scalar sets eppos2 (1-D,
%                     L2), settings_2d.eppos (2-D) and settings_nd.eppos (3
%                     or more variables) where unset;
%           epneg     double, the decay rate in every dimension: sets
%                     settings_2d.epneg where unset (3 or more variables
%                     read epneg itself);
%           sos_opts  'struct' of solver options;
%           sop       'struct' of options of the executives, each field
%                     optional:
%             slack       char, 'like' (default) or 'stock': each
%                         negativity constraint -D >= 0 is imposed with a
%                         slack sized from the positive operator of the LPI
%                         ('like'), or as D + N == 0 with N the positive
%                         slack declared from the degrees of sop.stock
%                         ('stock');
%             dw_Q        double (default 1), the raise of the weight
%                         degree of the 'like' slack in the non-coercive
%                         forms;
%             dw_P        double (default 2), the same in the other LPIs;
%             max_raise   double (default 2), the largest number of raises
%                         of the degree loop, which rebuilds the LPI with
%                         the slack weight raised by one while the solver
%                         returns 'unknown' or 'infeasible';
%             margin      double eps >= 0 (default 0), certify
%                         -D + eps*e_D >= 0 in place of -D >= 0, with
%                         e_D = A'*A and A the lift of the plain term of
%                         the slack;
%             witness     logical (default true), compute for a 1-D PIE the
%                         gain, decay rate or H2 norm of its Legendre
%                         discretization;
%             N_cheb      double (default 24), the degree of that
%                         discretization;
%             nfreq       double (default 240), its number of frequency
%                         samples;
%             keepdual    logical (default true), record the equality rows
%                         of the program (prog.sopeq) and assemble the dual
%                         kernel of the negativity constraint;
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
%                         is bounded per direction, theta_k + s_k <=
%                         Jmax_k(T) + cap_delta, with Jmax_k(T) the largest
%                         total degree in direction k of a kernel monomial
%                         of the PIE operator T (default []: no bound);
%             eppos_stab  double, a scalar or one entry per space, the
%                         margin of the positive operator of the stability
%                         and well-posedness LPIs in place of eppos, eppos2
%                         (1-D) and settings_2d.eppos (2-D) (default []:
%                         those fields);
%             deg_Q       1x3 'cell' (default []: lpi.var.deg), the degrees
%                         of Q in the stability form 'Q';
%             max_raise_Q double (default []: max_raise_feas), the raises
%                         of the stability form 'Q';
%             raise_Q     char, 'slack' (default), 'deg' or 'both': what a
%                         raise of the stability form 'Q' adds, the weight
%                         of the slacks, the first entry of deg_Q, or both;
%           lpi       'struct' with fields pos and var, the positive and
%                     free operators of the LPIs (default those of the
%                     'light' tier).
%
% OUTPUT
% - st:   'struct', ST with the fields of ST.sop above set (its other
%         fields kept), ST.sos_opts.solver set, the margin and
%         decay-rate fields of each dimension set from eppos and epneg
%         where unset, and lpi.pos and lpi.var set where missing.
%
% NOTES
% ST other than a 'struct', ST.sop.slack other than 'like' and 'stock',
% ST.sop.raise_Q other than 'slack', 'deg' and 'both', or a stock-format
% field (dd1, dd12, options1, options12, override1, ddZ, dd2, dd3,
% options2, options3, override2, sosineq_on; settings_nd.LF_int, LF_mult,
% Zop_deg) is an error that names its replacement. A storage row lpi.pos
% (field len_com) that also carries mult, int, joint or sep is an error:
% the row is read alone.
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

% Defaults of the fields of ST.sop
d = struct('slack','like','dw_Q',1,'dw_P',2,'max_raise',2,'margin',0,'witness',true, ...
           'N_cheb',24,'nfreq',240,'keepdual',true,'verbose',true, ...
           'dw_feas',[],'max_raise_feas',[],'eppos_stab',[],'cap_delta',[], ...
           'deg_Q',[],'max_raise_Q',[],'raise_Q','slack');
% Check ST
if ~isstruct(st),   error('exec_sop_settings:st','ST should be a settings struct (lpisettings).');   end
% Stock-format fields: an error naming the replacement
old = {'dd1','dd12','options1','options12','override1','ddZ','dd2','dd3', ...
       'options2','options3','override2','sosineq_on'};
new = [repmat({'use lpi.pos'},1,5),{'use lpi.var.deg{1}'}, ...
       repmat({'use sop.stock (LPISETTINGS_STOCK)'},1,5), ...
       {'remove it: the container executives pose no sosineq'}];
for k = find(isfield(st,old))
    error('exec_sop_settings:retired',['settings.%s was retired;'...
          ' %s.'],old{k},new{k});
end
if isfield(st,'settings_nd') && isstruct(st.settings_nd) && ...
        any(isfield(st.settings_nd,{'LF_int','LF_mult','Zop_deg'}))
    error('exec_sop_settings:retired',['settings_nd.LF_int, LF_mult and Zop_deg were '...
          'retired; use entry 3 of lpi.pos and lpi.var.'])
end
% lpi.pos and lpi.var, the operators of the LPIs: those of 'light' where
% missing
if ~isfield(st,'lpi') || ~all(isfield(st.lpi,{'pos','var'}))
    d0 = lpisettings();
    if ~isfield(st,'lpi'),          st.lpi = struct();      end
    if ~isfield(st.lpi,'pos'),      st.lpi.pos = d0.lpi.pos;    end
    if ~isfield(st.lpi,'var'),      st.lpi.var = d0.lpi.var;    end
end
% A storage row (len_com) with stock fields: the row is read alone, so
% mult, int, joint, sep would be ignored without notice
if isfield(st.lpi.pos,'len_com') && any(isfield(st.lpi.pos,{'mult','int','joint','sep'}))
    error('exec_sop_settings:row',['lpi.pos is a storage row (len_com): its fields ' ...
          'mult, int, joint, sep are not read; set len_inc, len_com, len_sos and cap ' ...
          '(see lpisettings).'])
end
% Each missing or empty field of ST.sop: its default
if ~isfield(st,'sop') || isempty(st.sop),   st.sop = struct();     end
fn = fieldnames(d);
for i = 1:numel(fn)
    if ~isfield(st.sop,fn{i}) || isempty(st.sop.(fn{i})),   st.sop.(fn{i}) = d.(fn{i});    end
end
% Check the slack mode
if ~any(strcmp(st.sop.slack,{'like','stock'}))
    error('exec_sop_settings:slack','sop.slack should be ''like'' or ''stock''.');
end
% The raise target of the stability form 'Q'
if ~any(strcmp(st.sop.raise_Q,{'slack','deg','both'}))
    error('exec_sop_settings:raise_Q', ...
          'sop.raise_Q should be ''slack'', ''deg'' or ''both''.');
end
% Solver: MOSEK when mosekopt is on the path, else SeDuMi
if ~isfield(st,'sos_opts') || isempty(st.sos_opts),    st.sos_opts = struct();    end
if ~isfield(st.sos_opts,'solver') || isempty(st.sos_opts.solver)
    if ~isempty(which('mosekopt')),     st.sos_opts.solver = 'mosek';
    else,                               st.sos_opts.solver = 'sedumi';
    end
end
% One margin and one decay rate for every dimension: a scalar eppos fills
% the margin fields left unset, eppos2 (1-D), settings_2d.eppos (2-D) and
% settings_nd.eppos (3 or more variables), a scalar epneg the field
% settings_2d.epneg (3 or more variables read epneg itself); a non-scalar
% eppos, a 1-D R^n margin, fills nothing.
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
end
