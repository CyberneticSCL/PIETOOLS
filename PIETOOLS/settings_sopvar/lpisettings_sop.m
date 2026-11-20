function st = lpisettings_sop(type,derivative_strictness,simplify,solver)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ST = LPISETTINGS_SOP(TYPE,DERIVATIVE_STRICTNESS,SIMPLIFY,SOLVER) returns
% the LPI settings structure of the tier TYPE, with the decay rate, the
% simplification and the solver set by the other inputs. The tiers set the
% storage row lpi.pos, the lengths of the monomial vectors of the basis
% operators of a positive operator declared without degrees, per
% direction, 0 = absent: d1 the s-length of the lower-integral operators
% (letter c), d2 that of the upper- or full-integral operators (letter b),
% d3(a,b,c) the theta-lengths of the letters a (multiplier), b and c, and
% Pi the pattern set of the letters:
%
%   TYPE       1-D d1,d2,d3(a,b,c)  2-D d2,d3(a,b),Pi   3-D d2,d3(a,b),Pi
%   extreme    0,0,[2 0 0]          0,[1 0],aa          0,[1 0],aaa
%   stripped   2,2,[2 2 2]          2,[1 2],uniform     1,[1 1],uniform
%   light      2,2,[2 2 2]          2,[1 2],separable   2,[1 2],uniform
%   heavy      2,2,[3 2 2]          2,[2 2],separable   2,[1 2],uniform
%   veryheavy  4,4,[4 4 4]          3,[2 3],separable   2,[1 2],uniform
%   custom     as light
%
% and the options of the executives as follows:
%
%   TYPE       dw_feas  dw_P     max_raise_feas  max_raise  cap_delta
%   extreme    [0 1 1]  1        [1 1 0]         0          [NaN NaN 0]
%   stripped   [0 1 1]  1        [1 1 0]         0          [NaN NaN 0]
%   light      [0 1 1]  [2 2 1]  [2 1 0]         [1 1 0]    [NaN NaN 0]
%   heavy      [1 2 2]  2        [2 2 0]         [1 1 0]    [NaN NaN 1]
%   veryheavy  [1 2 2]  2        [3 3 0]         [2 2 0]    [NaN NaN 2]
%   custom     as light
%
% A vector holds the entries of 1, 2, and 3 or more variables; NaN in
% cap_delta is no bound; dw_Q is 1 in every tier. The 1-D entry has the
% patterns 'all' (extreme: a alone), the 2-D and 3-D entries d1 = 0. The
% output options keepdual and witness, which the tiers do not carry, are
% set false here at extreme and stripped, and N_cheb, nfreq to 32, 480 at
% veryheavy; unset, they default to true, true, 24 and 240.
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
%             (a first-order GPU solver that decides no feasibility: its
%             solve is judged by the residual, never certified).
%
% OUTPUT
% - st:       'struct' with fields
%               lpi        'struct' with fields
%                 pos      'struct', the storage row: len_R (1, the
%                          identity on R^q), and 1x3 cells, entry k for a
%                          space in k variables (the last for more):
%                          len_inc (d1), len_com (d2), len_sos (d3, rows
%                          a, b, c), patterns (Pi), cap (Inf: none),
%                          len_cross (1: none), terms ('none'), offsets
%                          ([]);
%                 ineq     the slack options of an operator inequality: dD
%                          1 and dw 1, added to the lift and weight degrees
%                          read from the support; wmin, the least weight
%                          degree when the support is read from the target
%                          itself, [2 2 0] (0 at stripped, extreme); wcross
%                          0, no weight in a direction a space does not
%                          own; psatz {'auto','faces1','auto'}, the term
%                          sets, and psatz_offset []; Dmax false and
%                          cap_inc [] (off);
%                 var      'struct' with field deg, a 1x3 'cell' of degree
%                          specifications of the free operators; light
%                          {[2 2 2],1,1}, heavy {[4 2 2],2,1}, veryheavy
%                          {[6 6 6],2,1}, extreme {[2 2 0],1,1};
%               sop        'struct' of the options of the executives:
%                 slack          char, 'like': each negativity slack sized
%                                from the positive operator of the LPI;
%                 dw_feas, dw_P, dw_Q, max_raise_feas, max_raise,
%                 cap_delta      double, a scalar or one entry per
%                                dimension, as in the table above;
%                 deg_Q          1x3 'cell', the degrees of Q in the
%                                stability form 'Q' (light {[4 2 2],1,1});
%                 max_raise_Q    double, its raises (as max_raise_feas);
%                 raise_Q        char, 'slack': a raise adds slack weight;
%                 margin         double, 0: no regularizer on the certified
%                                inequality;
%                 keepdual, witness  logical, false (extreme, stripped
%                                only): no dual kernel, no 1-D witness;
%                 N_cheb, nfreq  double, 32, 480 (veryheavy only): the
%                                Legendre degree and number of frequencies
%                                of the witness;
%               eppos      double, 1e-4: the margin of the storage on
%                          every space, in every dimension and every LPI
%                          that imposes one;
%               epneg      double, DERIVATIVE_STRICTNESS, the decay rate in
%                          every dimension;
%               sos_opts   'struct' with fields simplify (0 or 1) and
%                          solver ('' unless SOLVER is given).
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
% Initial coding MMP, 10/28/2026

% Admissible values and defaults of the inputs
arguments
    type {mustBeMember(type,{'extreme','stripped','light','heavy','veryheavy','custom'})} = 'light';
    derivative_strictness {mustBeNonnegative} = 0;
    simplify {mustBeMember(simplify,{'','psimplify'})} = '';
    solver {mustBeMember(solver,{'','sedumi','mosek','sdpt3','sdpnalplus','cuadmm'})} = '';
end
% The settings of the tier
switch type
    case 'extreme',     st = settings_sopvar_extreme();
    case 'stripped',    st = settings_sopvar_stripped();
    case 'light',       st = settings_sopvar_light();
    case 'heavy',       st = settings_sopvar_heavy();
    case 'veryheavy',   st = settings_sopvar_veryheavy();
    case 'custom',      st = settings_sopvar_custom();
end
% Output options the tiers do not carry: keepdual and witness false at
% extreme and stripped, N_cheb 32 and nfreq 480 at veryheavy (unset, they
% default to true, true, 24 and 240)
switch type
    case {'extreme','stripped'},   st.sop.keepdual = false;   st.sop.witness = false;
    case 'veryheavy',              st.sop.N_cheb = 32;        st.sop.nfreq = 480;
end
% Decay rate, simplification and solver from the inputs
st.epneg = derivative_strictness;
st.sos_opts.simplify = strcmp(simplify,'psimplify');
if ~isempty(solver),    st.sos_opts.solver = solver;    end
end
