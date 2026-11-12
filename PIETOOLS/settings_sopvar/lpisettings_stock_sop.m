function [st,stock] = lpisettings_stock_sop(type,varargin)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% [ST,STOCK] = LPISETTINGS_STOCK_SOP(TYPE,...) returns the settings of the
% tier TYPE (LPISETTINGS_SOP) with the slack of the stock PIETOOLS executives,
% sop.slack = 'stock', for reproducing the baselines that compare the
% container executives with stock PIETOOLS. In one spatial variable the
% slack N = N1 + N2 is that of the stock level of the tier: N1 of dd2 and
% options2, N2 (the psatz term) of dd3 and options3, omitted when override2
% = 1; in two variables it is the stock slack of st.settings_2d; three or
% more variables have none.
%
% INPUT
% - type:   (optional) char, a tier of LPISETTINGS_SOP: 'extreme',
%           'stripped', 'light' (default), 'heavy', 'veryheavy' or 'custom';
% - ...:    (optional) the further inputs of LPISETTINGS_SOP: the decay
%           rate, the simplification and the solver.
%
% OUTPUT
% - st:     'struct', LPISETTINGS_SOP(TYPE) with sop.slack = 'stock' and
%           sop.stock = STOCK;
% - stock:  'struct', the slack alone, for a caller that builds it itself
%           (tl.slack of EXEC_TOOLS_SOP), with fields
%             override2  1 to omit N2;
%             options2, options3  'struct' with optional fields exclude
%                        (1x4: the R^n, multiplier, lower- and upper-
%                        integral bases), sep and psatz (POSLPIVAR_SOP);
%             dd2, dd3   {d1, d2, d3}: d1 the degree of the multiplier
%                        basis; d2, d3 those of the lower- and upper-
%                        integral bases, each [degree in the integration
%                        variable, in s, joint].
%
% NOTES
% The values are those the tiers held until 11/12/2026, the stock levels of
% M. Peet, S. Shivakumar and D. Jagt (2021-2022): stripped, light and custom
% the stock light level, heavy and veryheavy the stock heavy and veryheavy
% levels, extreme the stock extreme level.
%
% For support, contact M. Peet, Arizona State University at mpeet@asu.edu

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PIETOOLS - lpisettings_stock_sop
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
% Initial coding MMP, 11/12/2026: the stock slack leaves the settings tiers
%   (maintainer, 11/12/2026: "retire it, but create a utility function which
%   creates an ad hoc settings file which reproduces the baselines"); its
%   fields, top-level fields of every tier until 11/12/2026, are written out
%   here under sop.stock.

if nargin<1 || isempty(type),   type = 'light';     end
st = lpisettings_sop(type,varargin{:});
% The stock-format slack of each tier, as the tiers set it before 11/12/2026
switch type
    case {'stripped','light','custom'}          % the stock light level
        stock = struct('override2',0, ...
                       'options2',struct('exclude',[0 0 0 0]), ...
                       'dd2',{{2, [1 2 3], [1 2 3]}}, ...
                       'options3',struct('psatz',1,'exclude',[0 0 0 0]), ...
                       'dd3',{{1, [0 1 2], [0 1 2]}});
    case 'heavy'
        stock = struct('override2',0, ...
                       'options2',struct('psatz',0,'exclude',[0 0 0 0]), ...
                       'dd2',{{3, [2 2 3], [2 2 3]}}, ...
                       'options3',struct('psatz',1,'exclude',[0 0 0 0]), ...
                       'dd3',{{2, [1 1 2], [1 1 2]}});
    case 'veryheavy'
        stock = struct('override2',0, ...
                       'options2',struct('exclude',[0 0 0 0]), ...
                       'dd2',{{4, [5 5 8], [5 5 8]}}, ...
                       'options3',struct('psatz',1,'exclude',[0 0 0 0]), ...
                       'dd3',{{3, [4 4 7], [4 4 7]}});
    case 'extreme'                              % N2 omitted (override2 = 1)
        stock = struct('override2',1, ...
                       'options2',struct('exclude',[0 1 0 0],'sep',0), ...
                       'dd2',{{1, [1 0 1], [1 0 1]}}, ...
                       'options3',struct('psatz',1,'exclude',[0 1 0 1]), ...
                       'dd3',{{0, [0 -1 0], [-1 -1 0]}});
end
st.sop.slack = 'stock';
st.sop.stock = stock;

end
