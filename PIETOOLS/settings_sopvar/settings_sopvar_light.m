function st = settings_sopvar_light()
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ST = SETTINGS_SOPVAR_LIGHT() returns the settings structure of the
% 'light' tier of LPISETTINGS_SOP: the storage degrees of the stock light
% level in 1-D and 2-D, and the container options below.
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
%           settings_nd  'struct', the storage in 3 or more variables
%                    (EXEC_TOOLS_SOP): LF_int 1, LF_mult 0, eppos 1e-2;
%           sop      'struct' of container options, a vector holding one
%                    entry per dimension (1-D, 2-D, 3-D and higher):
%             slack 'like'       each negativity slack sized from the
%                                positive operator of the LPI;
%             dw_feas [0 1 1]    degree raise of that slack in stability
%                                and well-posedness;
%             dw_P [2 2 1], dw_Q 1   the same in the coercive and synthesis
%                                LPIs, and in the Q forms;
%             max_raise_feas [2 1 0], max_raise [1 1 0]   rebuilds with that
%                                degree one higher after UNKNOWN or
%                                INFEASIBLE, in stability and
%                                well-posedness, and in the other LPIs;
%             cap_delta [NaN NaN 0]  no bound in 1-D and 2-D; in 3-D and
%                                higher the slack bases bounded per
%                                direction by Jmax(T) (EXEC_SOP_SETTINGS);
%             eppos_stab 1e-2    margin of the positive operator in stability
%                                and well-posedness, on every space;
%             margin 0           no regularizer on the certified inequality;
%             keepdual true      assemble the dual kernel;
%             witness true, N_cheb 24, nfreq 240   compute the 1-D numerical
%                                counterpart, Legendre degree 24, 240
%                                frequencies.
%
% For support, contact M. Peet, Arizona State University at mpeet@asu.edu

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PIETOOLS - settings_sopvar_light
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
% MMP, 10/29/2026: max_raise_feas [2 1] (was 2): one raise in 2-D. In the
%                10/28 example-library baseline at this tier, the second
%                2-D raise ran only on PDE 36 and 38 (2D_Telegraph_Eq,
%                2D_Wave_Eq; infeasible on both paths at every pass), took
%                2,803 s of the tier's 3,819 s of solve time, and certified
%                nothing; the first 2-D raise certified PDE 34
%                (2D_KISS_Model). The 1-D raises are kept: 80 s, 35 passes
%                certified.
% MMP, 10/29/2026: 3-D entries (maintainer: every tier defined for 3-D, light
%                to run multi-state 3-D problems on this workstation):
%                settings_nd, the lean N-D storage of EXEC_TOOLS_SOP; a
%                third entry of dw_feas, dw_P, max_raise_feas, max_raise and
%                cap_delta for 3 or more variables: weight 1, no raise, the
%                per-direction bound Jmax(T) (cap_delta 0). Measured 10/29 on
%                3-D heat: m 28,956 (one field, stability) to 31,493 (an ODE
%                plus a field, H-infinity), the size of the 3-D benchmark
%                solved here; two coupled fields give m 166,572 (stability),
%                past MOSEK on this workstation. In 2-D this setting gives
%                gains 36-39% above exact; cap_delta 1 gives 0.13-0.15%.

% Storage degrees: the stock light levels, read by poslpivar_settings_sop
% (1-D) and poslpivar_settings_2d_sop (2-D).
st = settings_PIETOOLS_light();
st.settings_2d = settings_PIETOOLS_light_2D();
% Margins of the positive operator in the gain, norm and synthesis LPIs;
% stability and well-posedness take sop.eppos_stab instead.
st.eppos = 1e-4;    st.eppos2 = 1e-4;   st.epneg = 0;
st.settings_2d.eppos = 1e-2*ones(4,1);  st.settings_2d.epneg = 0;
st.sos_opts = struct('simplify',0);
% Container knobs (EXEC_SOP_SETTINGS).
% st.sop = struct('slack','like','dw_feas',[0 1],'dw_P',2,'dw_Q',1, ...
%                 'max_raise_feas',2,'max_raise',1,'eppos_stab',1e-2,'margin',0, ...
%                 'keepdual',true,'witness',true,'N_cheb',24,'nfreq',240);  % MMP, 10/29/2026 (was)
% max_raise_feas per dimension [1-D 2-D]: the second 2-D raise is dropped
% st.sop = struct('slack','like','dw_feas',[0 1],'dw_P',2,'dw_Q',1, ...
%                 'max_raise_feas',[2 1],'max_raise',1,'eppos_stab',1e-2,'margin',0, ...
%                 'keepdual',true,'witness',true,'N_cheb',24,'nfreq',240);  % MMP, 10/29/2026 (was)
% Storage in 3 or more variables (EXEC_TOOLS_SOP): the lean rule            % MMP, 10/29/2026
st.settings_nd = struct('LF_int',1,'LF_mult',0,'eppos',1e-2);               % MMP, 10/29/2026
% Entries [1-D 2-D 3-D] where a vector: 3-D at weight 1, no raise, bounded  % MMP, 10/29/2026
st.sop = struct('slack','like','dw_feas',[0 1 1],'dw_P',[2 2 1],'dw_Q',1, ...
                'max_raise_feas',[2 1 0],'max_raise',[1 1 0],'cap_delta',[NaN NaN 0], ...
                'eppos_stab',1e-2,'margin',0, ...
                'keepdual',true,'witness',true,'N_cheb',24,'nfreq',240);    % MMP, 10/29/2026
end
