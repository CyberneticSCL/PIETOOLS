function st = settings_sopvar_custom()
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ST = SETTINGS_SOPVAR_CUSTOM() returns the settings structure of the
% 'custom' tier of LPISETTINGS_SOP, a template to edit: as shipped, the
% settings of the 'light' tier, with the range of every field stated in
% the body.
%
% INPUT
% none.
%
% OUTPUT
% - st:   'struct' with fields
%           dd1, dd12, options1, options12, override1, ddZ, dd2, dd3,
%           options2, options3, override2, sosineq_on
%                    the 1-D degrees and options of settings_PIETOOLS_light;
%           settings_2d  'struct', settings_PIETOOLS_light_2D with eppos =
%                    1e-2*ones(4,1) and epneg = 0;
%           eppos, eppos2  double, 1e-4: the 1-D margins of the positive
%                    operator in the gain, norm and synthesis LPIs;
%           epneg    double, 0: the decay rate;
%           sos_opts 'struct' with simplify = 0;
%           sop      'struct' of container options: slack 'like', dw_feas
%                    [0 1], dw_P 2, dw_Q 1, max_raise_feas 2, max_raise 1,
%                    eppos_stab 1e-2, margin 0, keepdual true, witness
%                    true, N_cheb 24, nfreq 240, verbose true, each
%                    described where it is set.
%
% For support, contact M. Peet, Arizona State University at mpeet@asu.edu

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PIETOOLS - settings_sopvar_custom
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

% Storage degrees: any stock level, 1-D (settings_PIETOOLS_<level>: dd1,
% dd12, options1, options12, ddZ) and 2-D (settings_PIETOOLS_<level>_2D:
% LF_deg, LF_opts, Zop_deg, ...), or those fields edited.
st = settings_PIETOOLS_light();
st.settings_2d = settings_PIETOOLS_light_2D();
% Margins: eppos (R^n), eppos2 (L2) in 1-D, one per space of R x L2[s1] x
% L2[s2] x L2[s1,s2] in 2-D, for the gain, norm and synthesis LPIs.
st.eppos = 1e-4;    st.eppos2 = 1e-4;   st.epneg = 0;
st.settings_2d.eppos = 1e-2*ones(4,1);  st.settings_2d.epneg = 0;
% sop.slack = 'stock' needs the equality slack of the stock files: dd2,
% dd3 (1-D, sosineq_on 0) and eq_deg, eq_opts (2-D, use_sosineq 0). The
% stock light 2-D file sets use_sosineq 1 and declares no eq fields; take
% another 2-D level there. The default 'like' reads neither.
% Solver: sos_opts.solver 'mosek' or 'sedumi' (unset: MOSEK when installed,
% else SeDuMi); simplify 0 or 1.
st.sos_opts = struct('simplify',0);
% Container knobs (EXEC_SOP_SETTINGS).
st.sop = struct( ...
    'slack','like', ...         'like' | 'stock' (the stock slack, for comparison)
    'dw_feas',[0 1], ...        weight raise, stability and well-posedness, [1-D 2-D]
    'dw_P',2, ...               weight raise, coercive gains and norms, synthesis
    'dw_Q',1, ...               weight raise, the Q (non-coercive) forms
    'max_raise_feas',2, ...     raises on UNKNOWN / INFEASIBLE, stability, well-posedness
    'max_raise',1, ...          raises on UNKNOWN / INFEASIBLE, every other LPI
    'eppos_stab',1e-2, ...      storage margin of stability and well-posedness, every space
    'margin',0, ...             eps > 0 certifies -D + eps*e_D >= 0 instead of -D >= 0
    'keepdual',true, ...        assemble the dual kernel of the negativity constraint
    'witness',true, ...         1-D numerical counterpart: gain, rate or H2 norm
    'N_cheb',24, ...            its Legendre order
    'nfreq',240, ...            its frequency samples
    'verbose',true);
end
