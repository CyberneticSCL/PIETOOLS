function [prog,P] = poslpivar_settings_sop(prog,n,settings,role,varargin)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% [PROG,P] = POSLPIVAR_SETTINGS_SOP(PROG,X,SETTINGS,ROLE[,SIDE[,DOM]])
% declares the positive operator P = P1 + P2 that the 1-D LPI settings
% SETTINGS define for ROLE, on the SIDE spaces of the container X, each Pi
% the positive semidefinite operator Z'*(g*Q)*Z, Q >= 0, that
% POSLPIVAR_SOP(PROG,X,dd,options,SIDE,DOM) declares for the degrees dd
% and options of the settings:
%   ROLE 'lf'     P1 of settings.lpi.pos, entry 1: dd = {mult, [int int
%                 joint], [int int joint]}, options sep and, for int < 0,
%                 exclude of the integral bases; no P2;
%   ROLE 'slack'  P1 of dd2 and options2, P2 of dd3 and options3, P2
%                 omitted when override2 = 1, all in settings.sop.stock
%                 (LPISETTINGS_STOCK_SOP).
% No constraint other than P1 >= 0 and P2 >= 0 is added. An X that is not a
% 'copvar' or 'cdopvar' object, such as the dimension N of the retired
% POSLPIVAR, is an error (POSLPIVAR_SOP).
%
% INPUT
% - prog:      'struct', the LPI program to which the decision variables
%              are added;
% - X:         'copvar' or 'cdopvar' object whose SIDE spaces P maps, each
%              R^q or L2^m[s] in one spatial variable;
% - settings:  'struct' of LPI settings (LPISETTINGS_SOP) with the fields
%              of ROLE: lpi.pos ('lf': 'struct' with fields mult, int,
%              joint, one entry per dimension, and sep, as the tiers set
%              them) or sop.stock ('slack': 'struct' with fields dd2,
%              options2, override2, dd3, options3, from
%              LPISETTINGS_STOCK_SOP). Each dd is a 'cell' {d1,d2,d3} of monomial
%              degrees: d1 of the multiplier basis, d2 and d3 [degree in
%              the integration variable, in s, joint] of the lower and upper
%              integral bases. Each options is a 'struct' with fields psatz
%              (1: the weight (s-a)*(b-s)), exclude (1x4 logical: the R^q,
%              multiplier, lower and upper integral bases) and sep (1: one
%              integral over [a,b] in place of the lower and upper ones);
% - role:      'char', 'lf' or 'slack';
% - side:      (optional, X only) 'char', 'out' or 'in', the spaces of X;
%              omitted, X must have the same spaces and dimensions in and
%              out, and 'out' is used;
% - dom:       (optional, X only) the domain: 'double' 1x2 [a b], or
%              'struct' with fields vars ('cellstr') and dom (one row [a b]
%              per variable); default the variables and domain of X.
%
% OUTPUT
% - prog:      'struct', the program with the decision variables and Gram
%              constraints of P1 and then of P2 added;
% - P:         'cdopvar' object P1 + P2 (P1 when P2 is omitted).
%
% NOTES
% ROLE 'slack' without settings.sop.stock is an error, as is a ROLE other
% than 'lf' and 'slack'.
%
% See also POSLPIVAR_SOP, LPISETTINGS_SOP, LPI_EQ_SOP.
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
% MMP, 10/28/2026: Header to the header protocol, step comments in the body;
%                code unchanged.
% MMP, 11/08/2026: Comment pass (maintainer, 11/08/2026): the help no longer
%                offers the dimension N of the retired POSLPIVAR in place of
%                X (POSLPIVAR_SOP rejects it since 11/08/2026) and names
%                LPISETTINGS_SOP; code unchanged.
% MMP, 11/12/2026: The settings are native (maintainer, 11/12/2026: phase 2
%                of the knob reduction): ROLE 'lf' reads lpi.pos, entry 1,
%                and builds the stock-format dd and options that
%                POSLPIVAR_SOP takes from it, the same as the tiers' former
%                dd1 and options1, so the same P; P2 of dd12 is gone
%                (override1 was 1 in every tier). ROLE 'slack' reads the
%                stock slack from sop.stock (LPISETTINGS_STOCK_SOP), which
%                left the tiers; the check of sosineq_on, which left them
%                too, is removed.

% The settings fields of ROLE: degrees and options of P1, the override,
% degrees and options of P2
switch role
    case 'lf'
%       f = {'dd1','options1','override1','dd12','options12'};              % MMP, 11/12/2026 (was)
        % P1 of lpi.pos, entry 1: the multiplier basis of degree mult; the  % MMP, 11/12/2026
        % lower and upper integral bases of degree int in each variable and % MMP, 11/12/2026
        % joint in total, excluded when int < 0; no P2                      % MMP, 11/12/2026
        ps = settings.lpi.pos;  b = ps.int(1);  j = ps.joint(1);            % MMP, 11/12/2026
        settings = struct('dd1',{{ps.mult(1),[b b j],[b b j]}}, ...
            'options1',struct('sep',ps.sep,'exclude',[0 0 1 1]*(b<0)),'override1',1); % MMP, 11/12/2026
        f = {'dd1','options1','override1'};                                 % MMP, 11/12/2026
    case 'slack'
%       % sosineq_on (lpi_ineq in place of the equality) is not supported   % MMP, 11/12/2026 (was)
%       if isfield(settings,'sosineq_on') && settings.sosineq_on            % MMP, 11/12/2026 (was)
%           error('poslpivar_settings_sop:sosineq',['settings.sosineq_on = 1 '... % MMP, 11/12/2026 (was)
%                 'selects lpi_ineq, which has no container form here.'])   % MMP, 11/12/2026 (was)
%       end                                                                 % MMP, 11/12/2026 (was)
        % The stock slack of LPISETTINGS_STOCK_SOP                          % MMP, 11/12/2026
        if ~isfield(settings,'sop') || ~isfield(settings.sop,'stock')       % MMP, 11/12/2026
            error('poslpivar_settings_sop:stock',['The stock slack reads settings.sop.stock;'...
                  ' take the settings from LPISETTINGS_STOCK_SOP.'])        % MMP, 11/12/2026
        end                                                                 % MMP, 11/12/2026
        settings = settings.sop.stock;                                      % MMP, 11/12/2026
        f = {'dd2','options2','override2','dd3','options3'};
    otherwise
        error('poslpivar_settings_sop:role',"ROLE must be 'lf' or 'slack'.")
end
% P1; then, unless the override is 1, P2 and the sum P1 + P2
[prog,P] = poslpivar_sop(prog,n,settings.(f{1}),settings.(f{2}),varargin{:});
if settings.(f{3})~=1
    [prog,P2] = poslpivar_sop(prog,n,settings.(f{4}),settings.(f{5}),varargin{:});
    P = P + P2;
end
end
