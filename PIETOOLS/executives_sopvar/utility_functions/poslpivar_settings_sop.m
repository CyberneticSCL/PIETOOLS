function [prog,P] = poslpivar_settings_sop(prog,n,settings,role,varargin)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% [PROG,P] = POSLPIVAR_SETTINGS_SOP(PROG,N,SETTINGS,ROLE[,SIDE[,DOM]])
% declares the positive operator the 1-D executives build from their LPI
% settings, for either operator family (see 'poslpivar_sop' for N/X, SIDE
% and DOM):
%
%   ROLE 'lf'     the Lyapunov / storage operator
%       [prog,P1] = poslpivar(prog,Top.dim,dd1,options1);
%       if override1~=1, [prog,P2] = poslpivar(prog,Top.dim,dd12,options12);
%                        P = P1+P2; end
%   ROLE 'slack'  the positive slack of the equality branch
%       [prog,P1] = poslpivar(prog,Dop.dim,dd2,options2);
%       if override2~=1, [prog,P2] = poslpivar(prog,Dop.dim,dd3,options3);
%                        P = P1+P2; end
%
% as written in every 1-D executive (e.g. PIETOOLS_Hinf_gain: the
% 'poslpivar(prog,Top.dim,dd1,options1)' block, and the 'dd2,options2'
% block before its 'lpi_eq(prog,Deop+Dop,'symmetric')'). It imposes
% nothing: the caller writes the equality, with its sign.
%
% INPUT
% - prog, n:  as for 'poslpivar_sop';
% - settings: LPI settings struct (lpisettings / settings_PIETOOLS_*), with
%             the dd*, options* and override* fields of ROLE;
% - role:     'lf' or 'slack'.
%
% OUTPUT
% - prog:     the program with the Gram matrices declared, P1's first;
% - P:        P1, or P1 + P2 (in that operand order).
%
% NOTES
% ROLE 'slack' is the equality branch only; settings.sosineq_on = 1 selects
% lpi_ineq in the executives, which has no container form here, and is an
% error. Do not use it where an executive interleaves declarations (e.g.
% PIETOOLS_well_posedness declares P1, R1, P2, R2): the decision variables
% would come out permuted.
% Cost: one or two 'poslpivar_sop' calls; with override ~= 1 one operator
% sum, whose decision-list merge is O(q log q) in the q new variables.
%
% See also POSLPIVAR_SOP, LPISETTINGS, LPI_EQ_SOP.
%
% For support, contact M. Peet, Arizona State University at mpeet@asu.edu

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PIETOOLS - poslpivar_settings_sop
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
% Initial coding MMP, 10/06/2026. Replaces the test-folder cx_hinf_lf,
%                cx_h2_lf, cx_hinf_slack and cx_h2_slack (cx_exec, MMP
%                09/25/2026) and the P1(+P2) / N1(+N2) pairs written out in
%                the stability transcriptions. cx_hinf_slack imposed the
%                equality itself; returning N instead serves both forms.

switch role
    case 'lf'
        f = {'dd1','options1','override1','dd12','options12'};
    case 'slack'
        if isfield(settings,'sosineq_on') && settings.sosineq_on
            error('poslpivar_settings_sop:sosineq',['settings.sosineq_on = 1 '...
                  'selects lpi_ineq, which has no container form here.'])
        end
        f = {'dd2','options2','override2','dd3','options3'};
    otherwise
        error('poslpivar_settings_sop:role',"ROLE must be 'lf' or 'slack'.")
end
[prog,P] = poslpivar_sop(prog,n,settings.(f{1}),settings.(f{2}),varargin{:});
if settings.(f{3})~=1
    [prog,P2] = poslpivar_sop(prog,n,settings.(f{4}),settings.(f{5}),varargin{:});
    P = P + P2;
end
end
