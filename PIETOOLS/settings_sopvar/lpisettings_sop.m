function st = lpisettings_sop(type,derivative_strictness,simplify,solver)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ST = LPISETTINGS_SOP(TYPE,DERIVATIVE_STRICTNESS,SIMPLIFY,SOLVER) returns
% the settings structure of the container executives (PIETOOLS_*_sop) at
% the tier TYPE, the structure SETTINGS_SOPVAR_<TYPE> returns, with the
% decay rate, the simplification and the solver set by the other inputs.
% The tiers set the storage degrees from a stock level and the container
% options as follows:
%
%   TYPE       storage 1-D / 2-D   dw_feas  dw_P  dw_Q  max_raise_feas /  keepdual,
%                                                       max_raise         witness
%   stripped   light / light       [0 1]    1     1     1 / 0             false
%   light      light / light       [0 1]    2     1     2 / 1             true
%   heavy      heavy / stripped    [1 2]    2     1     2 / 1             true
%   veryheavy  veryheavy / heavy   [1 2]    2     1     3 / 2             true
%   custom     as light
%
% INPUT
% - type:     (optional) char, 'stripped', 'light' (default), 'heavy',
%             'veryheavy' or 'custom';
% - derivative_strictness: (optional) double >= 0, the decay rate epneg of
%             the stability LPIs, in 1-D and 2-D (default 0);
% - simplify: (optional) char, '' (default) or 'psimplify': simplify the
%             SDP before it is solved;
% - solver:   (optional) char, '' (default: MOSEK when installed, else
%             SeDuMi), 'mosek', 'sedumi', 'sdpt3' or 'sdpnalplus'.
%
% OUTPUT
% - st:       'struct' with fields
%               dd1, dd12, options1, options12, override1, ddZ, dd2, dd3,
%               options2, options3, override2, sosineq_on
%                          the 1-D degrees and options of the storage level's
%                          settings_PIETOOLS_<level>;
%               settings_2d  'struct', settings_PIETOOLS_<level>_2D of the 2-D
%                          storage level, with eppos (4x1 double, 1e-2) and
%                          epneg (DERIVATIVE_STRICTNESS);
%               eppos, eppos2  double, 1e-4: the 1-D margins of the positive
%                          operator in the gain, norm and synthesis LPIs;
%               epneg      double, DERIVATIVE_STRICTNESS;
%               sos_opts   'struct' with field simplify (0 or 1) and, when
%                          SOLVER is given, solver;
%               sop        'struct' of container options:
%                 slack          char, 'like': each negativity slack sized
%                                from the positive operator of the LPI;
%                 dw_feas        1x2 double, the degree raise of that slack
%                                in the stability and well-posedness LPIs,
%                                1-D and 2-D;
%                 dw_P, dw_Q     double, the same in the coercive (and
%                                synthesis) LPIs and in the Q forms;
%                 max_raise_feas, max_raise  double, how many times the
%                                program is rebuilt with that degree one
%                                higher after the solver returns UNKNOWN or
%                                INFEASIBLE (stability and well-posedness;
%                                the other LPIs);
%                 eppos_stab     double, 1e-2: the margin of the positive
%                                operator in the stability and well-posedness
%                                LPIs, on every space;
%                 margin         double, 0: no regularizer on the certified
%                                inequality;
%                 keepdual       logical: assemble the dual kernel of the
%                                negativity constraint;
%                 witness        logical: compute the 1-D numerical
%                                counterpart of the certificate;
%                 N_cheb, nfreq  double, its Legendre degree (24; 32 at
%                                veryheavy) and number of frequencies (240;
%                                480 at veryheavy).
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

% Admissible values and defaults of the inputs
arguments
    type {mustBeMember(type,{'stripped','light','heavy','veryheavy','custom'})} = 'light';
    derivative_strictness {mustBeNonnegative} = 0;
    simplify {mustBeMember(simplify,{'','psimplify'})} = '';
    solver {mustBeMember(solver,{'','sedumi','mosek','sdpt3','sdpnalplus'})} = '';
end
% The settings of the tier
switch type
    case 'stripped',    st = settings_sopvar_stripped();
    case 'light',       st = settings_sopvar_light();
    case 'heavy',       st = settings_sopvar_heavy();
    case 'veryheavy',   st = settings_sopvar_veryheavy();
    case 'custom',      st = settings_sopvar_custom();
end
% Decay rate, simplification and solver from the inputs
st.epneg = derivative_strictness;   st.settings_2d.epneg = derivative_strictness;
st.sos_opts.simplify = strcmp(simplify,'psimplify');
if ~isempty(solver),    st.sos_opts.solver = solver;    end
end
