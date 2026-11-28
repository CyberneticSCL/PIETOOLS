function st = settings_PIETOOLS_veryheavy()
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ST = SETTINGS_PIETOOLS_VERYHEAVY() returns the LPI settings of the
% 'veryheavy' tier:
% the degrees of the basis operators of the positive operator P of an LPI
% declared without degrees (the storage row lpi.pos), the slack options of
% an operator inequality (lpi.ineq), the degrees of the free operators
% (lpi.var), the options of the executives (sop), the margin eppos of P,
% the decay rate epneg and the solver options. Entry k of a 1 x 3 cell, or
% column k of a row, applies to a space or an inequality in k variables,
% the last entry to every higher k.
%
% INPUT
% none.
%
% OUTPUT
% - st:   'struct' with fields
%           lpi       'struct' with fields pos (len_R, len_inc, len_com,
%                     len_sos, patterns, cap, len_cross, terms, offsets),
%                     ineq (dD, dw, wmin, wcross, psatz, Dmax, cap_inc,
%                     psatz_offset) and var (field deg);
%           sop       'struct' with fields slack, dw_feas, dw_P, dw_Q,
%                     max_raise_feas, max_raise, cap_delta, deg_Q,
%                     max_raise_Q, raise_Q, margin;
%           eppos     'double', the margin of P on every space, in every
%                     dimension;
%           epneg     'double', the decay rate, in every dimension;
%           sos_opts  'struct' with fields simplify (0 or 1) and solver.
%         The body gives the value of each field and what it does.
%
% For support, contact M. Peet, Arizona State University at mpeet@asu.edu

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PIETOOLS - settings_PIETOOLS_veryheavy
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
% Initial coding MMP, 10/28/2026

% Veryheavy tier. Entry k (cell) or column k (row): a space or inequality in
% k variables, the last entry for every higher k.
% 1. The basis operators of P declared without degrees: lengths of the
% monomial vectors, 0 = absent; letters a the multiplier, b the upper
% integral (the full integral where no space carries c), c the lower
% integral; s the variables of the space, theta those of the registry.
st.lpi.pos = struct( ...
   'len_R',    1, ...                                   % R^q: the identity
   'len_inc',  {{4,       [0 0],         [0 0 0]}}, ...             % s-length of c
   'len_com',  {{4,       [3 3],         [2 2 2]}}, ...             % s-length of b
   'len_sos',  {{[4;4;4], [2 2;3 3;0 0], [1 1 1;2 2 2;0 0 0]}}, ... % theta-lengths of a;b;c
   'patterns', {{'all',   'separable',   'uniform'}}, ...           % letter patterns
   'cap',      {{Inf,     [Inf Inf],     [Inf Inf Inf]}}, ...       % theta_k+s_k caps
   'len_cross', {{1,      1,             1}}, ...        % non-owned, 1: none
   'terms',    {{'none',  'none',        'none'}}, ...              % no Psatz term
   'offsets',  {{[],      [],            []}});                     % no term offsets
% 2. The slack of an operator inequality; degrees.
st.lpi.ineq = struct( ...
   'dD',           1, ...         % lift raise: D = Dmin + dD
   'dw',           1, ...         % weight raise, calls without dw
   'wmin',         [2 2 0], ...   % weight floor, support read from the target
   'wcross',       0, ...         % weight, non-owned dirs; 0: none
   'psatz',        {{'auto','faces1','auto'}}, ... % term sets; 2-D 'faces1'
   'Dmax',         false, ...     % lift at the support ceiling; off
   'cap_inc',      [], ...        % triangular cap; 'auto' = Dmax-1; off
   'psatz_offset', []);           % read with a row of codes only
% 3. The free operators; degrees.
st.lpi.var = struct('deg', {{[6 6 6], 2, 1}});  % [d1 d2 d3]: mult, out, in d1; int [d3 d2]
% 4. The executives.
st.sop = struct( ...
   'slack',          'like', ...      % slack sized from the storage
   'dw_feas',        [1 2 2], ...     % weight raise: stability, w.-p.
   'dw_P',           2, ...           % weight raise: gains, synthesis
   'dw_Q',           1, ...           % weight raise: the Q forms
   'max_raise_feas', [3 3 0], ...     % raises: stability, w.-p.
   'max_raise',      [2 2 0], ...     % raises: every other LPI
   'cap_delta',      [NaN NaN 2], ... % Jmax(T) + cap_delta; NaN: none
   'deg_Q',          {{[6 6 6], 2, 1}}, ... % Q, direct stability form
   'max_raise_Q',    [3 3 0], ...     % its raises
   'raise_Q',        'slack', ...     % a raise adds: slack weight
   'margin',         0);           % eps of -D + eps*e_D >= 0
% 5. Tolerances and solver.
st.eppos = 1e-4;                      % P >= eppos*I on every space
st.epneg = 0;                         % decay rate
st.sos_opts = struct( ...
   'simplify', 0, ...                 % 1: simplify the SDP
   'solver',   '');                   % '': MOSEK if found, else SeDuMi
end
