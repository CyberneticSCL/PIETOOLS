function st = settings_sopvar_extreme()
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ST = SETTINGS_SOPVAR_EXTREME() returns the settings structure of the
% 'extreme' tier of LPISETTINGS_SOP: the storage degrees of the stock
% extreme level in 1-D and 2-D, the lowest of the stock levels, and the
% container options below, those of the stripped tier.
%
% INPUT
% none.
%
% OUTPUT
% - st:   'struct' with fields
%           dd1, dd12, options1, options12, override1, ddZ, dd2, dd3,
%           options2, options3, override2, sosineq_on
%                    the 1-D degrees and options of the extreme level;
%           settings_2d  'struct', the 2-D extreme level with eppos =
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
%             dw_P 1, dw_Q 1     the same in the coercive and synthesis
%                                LPIs, and in the Q forms;
%             max_raise_feas [1 1 0], max_raise 0   rebuilds with that
%                                degree one higher after UNKNOWN or
%                                INFEASIBLE, in stability and
%                                well-posedness, and in the other LPIs;
%             cap_delta [NaN NaN 0]  no bound in 1-D and 2-D; in 3-D and
%                                higher the slack bases bounded per
%                                direction by Jmax(T) (EXEC_SOP_SETTINGS);
%             eppos_stab 1e-2    margin of the positive operator in stability
%                                and well-posedness, on every space;
%             margin 0           no regularizer on the certified inequality;
%             keepdual false     no dual kernel;
%             witness false      no numerical counterpart (N_cheb 24, nfreq
%                                240 unused).
%
% For support, contact M. Peet, Arizona State University at mpeet@asu.edu

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PIETOOLS - settings_sopvar_extreme
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
% Initial coding MMP, 11/08/2026: the container counterpart of the stock
%                'extreme' level, which had no tier (maintainer,
%                11/08/2026): the stock extreme storage degrees
%                (STORAGE_EXTREME, STORAGE_EXTREME_2D in
%                settings_sopvar/private) with the margins, the N-D storage
%                and the container options of the stripped tier.

% Storage degrees: the stock extreme levels (see the header).
st = storage_extreme();
st.settings_2d = storage_extreme_2D();
% Margins of the positive operator in the gain, norm and synthesis LPIs;
% stability and well-posedness take sop.eppos_stab instead.
st.eppos = 1e-4;    st.eppos2 = 1e-4;   st.epneg = 0;
st.settings_2d.eppos = 1e-2*ones(4,1);  st.settings_2d.epneg = 0;
st.sos_opts = struct('simplify',0);
% Storage in 3 or more variables (EXEC_TOOLS_SOP): the lean rule of every tier
st.settings_nd = struct('LF_int',1,'LF_mult',0,'eppos',1e-2);
% Container knobs (EXEC_SOP_SETTINGS) of the stripped tier, entries
% [1-D 2-D 3-D] where a vector: no dual kernel, no witness, no raise in the
% gain, norm and synthesis LPIs.
st.sop = struct('slack','like','dw_feas',[0 1 1],'dw_P',1,'dw_Q',1, ...
                'max_raise_feas',[1 1 0],'max_raise',0,'cap_delta',[NaN NaN 0], ...
                'eppos_stab',1e-2,'margin',0, ...
                'keepdual',false,'witness',false,'N_cheb',24,'nfreq',240);
end
