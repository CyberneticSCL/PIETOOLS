function [prog,P] = poslpivar_settings(prog,n,settings,role,varargin)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% [PROG,P] = POSLPIVAR_SETTINGS(PROG,X,SETTINGS,ROLE[,SIDE[,DOM]])
% declares the positive operator P = P1 + P2 that the 1-D LPI settings
% SETTINGS define for ROLE, on the SIDE spaces of the container X, each Pi
% the positive semidefinite operator Z'*(g*Q)*Z, Q >= 0, with Z the
% identity on R^q and, on L2^m[s], the multiplier, lower- and
% upper-integral basis operators of the degrees dd that the options keep,
% g = (s-a)*(b-s) for options.psatz = 1 and g = 1 otherwise:
%   ROLE 'lf'     P1 of settings.lpi.pos, entry 1: dd = {mult, [int int
%                 joint], [int int joint]}, options sep and, for int < 0,
%                 exclude of the integral bases; no P2;
%   ROLE 'slack'  P1 of dd2 and options2, P2 of dd3 and options3, P2
%                 omitted when override2 = 1, all in settings.sop.stock.
% No constraint other than P1 >= 0 and P2 >= 0 is added. An X that is not an
% 'opvar' or 'dopvar' object is an error.
%
% INPUT
% - prog:      'struct', the LPI program to which the decision variables
%              are added;
% - X:         'opvar' or 'dopvar' object whose SIDE spaces P maps, each
%              R^q or L2^m[s] in one spatial variable;
% - settings:  'struct' of LPI settings with the fields of ROLE: lpi.pos
%              ('lf': 'struct' with fields mult, int, joint, one entry per
%              dimension, and sep) or sop.stock ('slack': 'struct' with
%              fields dd2, options2, override2, dd3, options3). Each dd is a
%              'cell' {d1,d2,d3} of monomial
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
% - P:         'dopvar' object P1 + P2 (P1 when P2 is omitted).
%
% NOTES
% ROLE 'slack' without settings.sop.stock is an error, as is a ROLE other
% than 'lf' and 'slack'.
%
% See also POSLPIVAR_LEGACY, LPISETTINGS, LPI_EQ.
%
% For support, contact M. Peet, Arizona State University at mpeet@asu.edu

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PIETOOLS - poslpivar_settings
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
% Initial coding MMP, 10/06/2026

% The settings fields of ROLE: degrees and options of P1, the override,
% degrees and options of P2
switch role
    case 'lf'
        % P1 of lpi.pos, entry 1: the multiplier basis of degree mult; the
        % lower and upper integral bases of degree int in each variable and
        % joint in total, excluded when int < 0; no P2
        ps = settings.lpi.pos;  b = ps.int(1);  j = ps.joint(1);
        settings = struct('dd1',{{ps.mult(1),[b b j],[b b j]}}, ...
            'options1',struct('sep',ps.sep,'exclude',[0 0 1 1]*(b<0)),'override1',1);
        f = {'dd1','options1','override1'};
    case 'slack'
        % The 'stock' slack, settings.sop.stock
        if ~isfield(settings,'sop') || ~isfield(settings.sop,'stock')
            error('poslpivar_settings:stock',['The stock slack reads settings.sop.stock;'...
                  ' take the settings from LPISETTINGS_STOCK.'])
        end
        settings = settings.sop.stock;
        f = {'dd2','options2','override2','dd3','options3'};
    otherwise
        error('poslpivar_settings:role',"ROLE must be 'lf' or 'slack'.")
end
% P1; then, unless the override is 1, P2 and the sum P1 + P2
[prog,P] = poslpivar_legacy(prog,n,settings.(f{1}),settings.(f{2}),varargin{:});
if settings.(f{3})~=1
    [prog,P2] = poslpivar_legacy(prog,n,settings.(f{4}),settings.(f{5}),varargin{:});
    P = P + P2;
end
end
