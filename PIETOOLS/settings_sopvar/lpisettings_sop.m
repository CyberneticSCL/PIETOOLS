function st = lpisettings_sop(type,derivative_strictness,simplify,solver)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ST = LPISETTINGS_SOP(TYPE,DERIVATIVE_STRICTNESS,SIMPLIFY,SOLVER) returns
% the settings structure of the container executives (PIETOOLS_*_sop) at
% the tier TYPE, the structure SETTINGS_SOPVAR_<TYPE> returns, with the
% decay rate, the simplification and the solver set by the other inputs.
% The tiers set the storage row lpi.pos (POSCOPVAR with DEG omitted;        % MMP, 11/17/2026
% settings report 11/16/2026, Sec. 6), lengths per direction, 0 = absent:   % MMP, 11/17/2026
%
%   TYPE       1-D d1,d2,d3(a,b,c)  2-D d2,d3(a,b),Pi   3-D d2,d3(a,b),Pi   % MMP, 11/17/2026
%   extreme    0,0,[2 0 0]          0,[1 0],aa          0,[1 0],aaa         % MMP, 11/17/2026
%   stripped   2,2,[2 2 2]          2,[1 2],uniform     1,[1 1],uniform     % MMP, 11/17/2026
%   light      2,2,[2 2 2]          2,[1 2],separable   2,[1 2],uniform     % MMP, 11/17/2026
%   heavy      2,2,[3 2 2]          2,[2 2],separable   2,[1 2],uniform     % MMP, 11/17/2026
%   veryheavy  4,4,[4 4 4]          3,[2 3],separable   2,[1 2],uniform     % MMP, 11/17/2026
%   custom     as light                                                     % MMP, 11/17/2026
%
% and the container options as follows:                                     % MMP, 11/17/2026
%
%   TYPE       dw_feas  dw_P     max_raise_feas  max_raise  cap_delta       % MMP, 11/17/2026
%   extreme    [0 1 1]  1        [1 1 0]         0          [NaN NaN 0]     % MMP, 11/17/2026
%   stripped   [0 1 1]  1        [1 1 0]         0          [NaN NaN 0]     % MMP, 11/17/2026
%   light      [0 1 1]  [2 2 1]  [2 1 0]         [1 1 0]    [NaN NaN 0]     % MMP, 11/17/2026
%   heavy      [1 2 2]  2        [2 2 0]         [1 1 0]    [NaN NaN 1]     % MMP, 11/17/2026
%   veryheavy  [1 2 2]  2        [3 3 0]         [2 2 0]    [NaN NaN 2]     % MMP, 11/17/2026
%   custom     as light                                                     % MMP, 11/17/2026
%
% A vector holds the entries of 1, 2, and 3 or more variables; NaN in       % MMP, 11/17/2026
% cap_delta is no bound; dw_Q is 1 in every tier. The 1-D entry has the     % MMP, 11/17/2026
% patterns 'all' (extreme: a alone), the 2-D and 3-D entries d1 = 0. The    % MMP, 11/17/2026
% output options keepdual and witness, which the tiers do not carry, are    % MMP, 11/17/2026
% set false here at extreme and stripped, and N_cheb, nfreq 32, 480 at      % MMP, 11/17/2026
% veryheavy (the tier values until 11/17/2026); EXEC_SOP_SETTINGS supplies  % MMP, 11/17/2026
% the others (true, true, 24, 240). The 'stock' slack, which the tiers      % MMP, 11/17/2026
% held until 11/12/2026, comes from LPISETTINGS_STOCK_SOP.                  % MMP, 11/17/2026
%
% INPUT
% - type:     (optional) char, 'extreme', 'stripped', 'light' (default),
%             'heavy', 'veryheavy' or 'custom';
% - derivative_strictness: (optional) double >= 0, the decay rate epneg of
%             the stability LPIs, in every dimension (default 0);
% - simplify: (optional) char, '' (default) or 'psimplify': simplify the
%             SDP before it is solved;
% - solver:   (optional) char, '' (default: MOSEK when installed, else
%             SeDuMi), 'mosek', 'sedumi', 'sdpt3', 'sdpnalplus' or 'cuadmm'
%             (the GPU executable of CUADMM_SEDUMI; the executives return
%             status 'residual', never certified).
%
% OUTPUT
% - st:       'struct' with fields (settings report 11/16/2026, Sec. 5):    % MMP, 11/17/2026
%               lpi        'struct' with fields                             % MMP, 11/17/2026
%                 pos      'struct', the storage row read by POSCOPVAR      % MMP, 11/17/2026
%                          with DEG omitted: len_R (1, the identity on      % MMP, 11/17/2026
%                          R^q), and 1x3 cells, entry k for a space in k    % MMP, 11/17/2026
%                          variables (the last for more): len_inc (d1),     % MMP, 11/17/2026
%                          len_com (d2), len_sos (d3, rows a, b, c),        % MMP, 11/17/2026
%                          patterns (Pi), cap (Inf: none), len_cross (1:    % MMP, 11/17/2026
%                          none), terms ('none'), offsets ([]);             % MMP, 11/17/2026
%                 ineq     the options of LPI_INEQ_SOP: dD 1 and dw 1,      % MMP, 11/17/2026
%                          added to the lift and weight degrees read from   % MMP, 11/17/2026
%                          the support; wmin, the least weight degree when  % MMP, 11/17/2026
%                          the support is read from the target itself,      % MMP, 11/17/2026
%                          [2 2 0] (0 at stripped, extreme); wcross 0, no   % MMP, 11/17/2026
%                          weight in a direction a space does not own;      % MMP, 11/17/2026
%                          psatz 'auto' and psatz_offset [], the rule;      % MMP, 11/17/2026
%                 var      'struct' with field deg, a 1x3 'cell' of degree  % MMP, 11/17/2026
%                          specifications of the free operator              % MMP, 11/17/2026
%                          (LPIVAR_CDOPVAR); light {[2 2 2],1,1}, heavy     % MMP, 11/17/2026
%                          {[4 2 2],2,1}, veryheavy {[6 6 6],2,1},          % MMP, 11/17/2026
%                          extreme {[2 2 0],1,1};                           % MMP, 11/17/2026
%               sop        'struct' of container options:                   % MMP, 11/17/2026
%                 slack          char, 'like': each negativity slack sized  % MMP, 11/17/2026
%                                from the positive operator of the LPI;     % MMP, 11/17/2026
%                 dw_feas, dw_P, dw_Q, max_raise_feas, max_raise,           % MMP, 11/17/2026
%                 cap_delta      double, a scalar or one entry per          % MMP, 11/17/2026
%                                dimension, as in the table above           % MMP, 11/17/2026
%                                (EXEC_SOP_SETTINGS);                       % MMP, 11/17/2026
%                 margin         double, 0: no regularizer on the certified % MMP, 11/17/2026
%                                inequality;                                % MMP, 11/17/2026
%                 keepdual, witness  logical, false (extreme, stripped      % MMP, 11/17/2026
%                                only): no dual kernel, no 1-D witness;     % MMP, 11/17/2026
%                 N_cheb, nfreq  double, 32, 480 (veryheavy only): the      % MMP, 11/17/2026
%                                Legendre degree and number of frequencies  % MMP, 11/17/2026
%                                of the witness;                            % MMP, 11/17/2026
%               eppos      double, 1e-4: the margin of the storage on       % MMP, 11/17/2026
%                          every space, in every dimension and every LPI    % MMP, 11/17/2026
%                          that imposes one (EXEC_SOP_SETTINGS);            % MMP, 11/17/2026
%               epneg      double, DERIVATIVE_STRICTNESS, the decay rate in % MMP, 11/17/2026
%                          every dimension;                                 % MMP, 11/17/2026
%               sos_opts   'struct' with fields simplify (0 or 1) and       % MMP, 11/17/2026
%                          solver ('' unless SOLVER is given).              % MMP, 11/17/2026
%
% For support, contact M. Peet, Arizona State University at mpeet@asu.edu

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PIETOOLS - lpisettings_sop
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
% Initial coding MMP, 10/28/2026. Tier values from the Posopvar, New
%                Posopvar and cuADMM sessions (10/28/2026) and the
%                measurements in README.md.
% MMP, 10/29/2026: Help text: light max_raise_feas [2 1] (one raise in 2-D,
%                settings_sopvar_light), and the per-dimension reading of
%                dw_feas and max_raise_feas. Code unchanged.
% MMP, 10/29/2026: Help text: the field cap_delta of EXEC_SOP_SETTINGS
%                (empty in every tier). Code unchanged.
% MMP, 10/29/2026: Help text: the 3-D entries of every tier (dw_feas, dw_P,
%                max_raise_feas, max_raise, cap_delta) and settings_nd; the
%                tiers now set cap_delta, so 'empty in every tier' in the
%                entry above no longer holds. Code unchanged.
% MMP, 11/07/2026: solver 'cuadmm' admitted. The arguments block rejected
%                it, although sossolve has the option since 10/30/2026 and
%                the executives classify its exit ('residual').
% MMP, 11/08/2026: Tier 'extreme' (SETTINGS_SOPVAR_EXTREME), the counterpart
%                of the stock extreme level, so that every level name of
%                the stock LPISETTINGS has a tier; help text updated.
% MMP, 11/09/2026: The decay rate is set once, in epneg (maintainer: "I see
%                no reason to treat the eppos and epneg options differently
%                for 1D,2D and ND"); EXEC_SOP_SETTINGS sets settings_2d.epneg
%                from it. The tiers set one margin eppos, 1e-4; help text
%                updated (eppos_stab, eppos2 and the 2-D and N-D margins
%                removed).
% MMP, 11/10/2026: Help text: the field lpi, the options of LPI_INEQ_SOP
%                that every tier now sets (st.lpi.ineq). Code unchanged.
% MMP, 11/12/2026: Help text: the tiers carry the operators in native form,
%                lpi.pos and lpi.var, in place of the stock-format 1-D fields
%                and settings_nd; the 'stock' slack left for
%                LPISETTINGS_STOCK_SOP (phase 2 of the knob reduction). Code
%                unchanged.
% MMP, 11/17/2026: The tiers take the structure of the settings report
%                11/16/2026, Sec. 5 and 6 (work package 5): help text for
%                the storage row lpi.pos and the fields that left
%                (settings_2d; lpi.pos mult, int, joint, sep; sop
%                keepdual, witness, N_cheb, nfreq). Code: keepdual and
%                witness false at extreme and stripped, N_cheb 32 and
%                nfreq 480 at veryheavy, the values the tier files carried
%                until 11/17/2026, so the output of the executives at each
%                tier is unchanged.

% Admissible values and defaults of the inputs
arguments
%   type {mustBeMember(type,{'stripped','light','heavy','veryheavy','custom'})} = 'light'; % MMP, 11/08/2026 (was)
    type {mustBeMember(type,{'extreme','stripped','light','heavy','veryheavy','custom'})} = 'light'; % MMP, 11/08/2026
    derivative_strictness {mustBeNonnegative} = 0;
    simplify {mustBeMember(simplify,{'','psimplify'})} = '';
%   solver {mustBeMember(solver,{'','sedumi','mosek','sdpt3','sdpnalplus'})} = ''; % MMP, 11/07/2026 (was)
    solver {mustBeMember(solver,{'','sedumi','mosek','sdpt3','sdpnalplus','cuadmm'})} = ''; % MMP, 11/07/2026
end
% The settings of the tier
switch type
    case 'extreme',     st = settings_sopvar_extreme();                     % MMP, 11/08/2026
    case 'stripped',    st = settings_sopvar_stripped();
    case 'light',       st = settings_sopvar_light();
    case 'heavy',       st = settings_sopvar_heavy();
    case 'veryheavy',   st = settings_sopvar_veryheavy();
    case 'custom',      st = settings_sopvar_custom();
end
% Output options the tiers no longer carry (settings report, Sec. 5):       % MMP, 11/17/2026
% their former tier values where they differ from the defaults of           % MMP, 11/17/2026
% EXEC_SOP_SETTINGS (keepdual, witness true; N_cheb 24, nfreq 240)          % MMP, 11/17/2026
switch type                                                                 % MMP, 11/17/2026
    case {'extreme','stripped'},   st.sop.keepdual = false;   st.sop.witness = false; % MMP, 11/17/2026
    case 'veryheavy',              st.sop.N_cheb = 32;        st.sop.nfreq = 480; % MMP, 11/17/2026
end                                                                         % MMP, 11/17/2026
% Decay rate, simplification and solver from the inputs
% st.epneg = derivative_strictness;   st.settings_2d.epneg = derivative_strictness; % MMP, 11/09/2026 (was)
st.epneg = derivative_strictness;                                           % MMP, 11/09/2026
st.sos_opts.simplify = strcmp(simplify,'psimplify');
if ~isempty(solver),    st.sos_opts.solver = solver;    end
end
