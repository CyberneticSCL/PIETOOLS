function [prog,P] = poslpivar_settings_sop(prog,n,settings,role,varargin)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% [PROG,P] = POSLPIVAR_SETTINGS_SOP(PROG,X,SETTINGS,ROLE[,SIDE[,DOM]])
% declares the positive operator P = P1 + P2 that the 1-D LPI settings
% SETTINGS define for ROLE, on the SIDE spaces of the container X, each Pi
% the positive semidefinite operator Z'*(g*Q)*Z, Q >= 0, that
% POSLPIVAR_SOP(PROG,X,dd,options,SIDE,DOM) declares for the degrees dd
% and options of the settings:
%   ROLE 'lf'     P1 of dd1 and options1, P2 of dd12 and options12, P2
%                 omitted when override1 = 1;
%   ROLE 'slack'  P1 of dd2 and options2, P2 of dd3 and options3, P2
%                 omitted when override2 = 1.
% With the dimension N of an 'opvar' object in place of X, P1 and P2 are
% the 'dopvar' objects of POSLPIVAR(PROG,N,dd,options). No constraint other
% than P1 >= 0 and P2 >= 0 is added.
%
% INPUT
% - prog:      'struct', the LPI program to which the decision variables
%              are added;
% - X:         'copvar' or 'cdopvar' object whose SIDE spaces P maps, each
%              R^q or L2^m[s] in one spatial variable; or
%   N:         'double', the dimensions [n0;n1] or the 2x2 .dim of an
%              'opvar' object;
% - settings:  'struct' of 1-D LPI settings (LPISETTINGS,
%              settings_PIETOOLS_*) with the fields of ROLE: dd1, options1,
%              override1, dd12, options12 ('lf') or dd2, options2,
%              override2, dd3, options3 ('slack'), and (optional)
%              sosineq_on. Each dd is a 'cell' {d1,d2,d3} of monomial
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
% - P:         'cdopvar' object P1 + P2 (P1 when P2 is omitted); with N, a
%              'dopvar' object.
%
% NOTES
% ROLE 'slack' with settings.sosineq_on nonzero is an error, as is a ROLE
% other than 'lf' and 'slack'.
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
% MMP, 10/28/2026: Header to the header protocol (what the function does,
%                inputs, outputs), step comments in the body; code
%                unchanged. The deleted header text was the transcription
%                of the executives' code, references to the executives and
%                a cost note.

% The settings fields of ROLE: degrees and options of P1, the override,
% degrees and options of P2
switch role
    case 'lf'
        f = {'dd1','options1','override1','dd12','options12'};
    case 'slack'
        % sosineq_on (lpi_ineq in place of the equality) is not supported
        if isfield(settings,'sosineq_on') && settings.sosineq_on
            error('poslpivar_settings_sop:sosineq',['settings.sosineq_on = 1 '...
                  'selects lpi_ineq, which has no container form here.'])
        end
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
