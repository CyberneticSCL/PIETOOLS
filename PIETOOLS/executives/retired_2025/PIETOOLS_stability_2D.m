function varargout = PIETOOLS_stability_2D(varargin) %#ok<STOUT> the call only errors
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PIETOOLS_STABILITY_2D - an executive of PIETOOLS 2025, retired. A call
% stops with an error that gives the call of this release:
%
%   PIETOOLS 2025:
%     [prog,Pop,Qop,solve_val] = PIETOOLS_stability_2D(PIE,settings)
%   this release:
%     [prog,P,info] = PIETOOLS_stability(PIE,settings,'P')
%
% The form 'P' poses the LPI of PIETOOLS_stability_2D. The executive takes
% a PIE in any number of spatial variables; the solver status is
% info.status. See PIETOOLS_stability for the forms and the outputs, and
% LPISETTINGS for SETTINGS.
%
% PIETOOLS - PIETOOLS_stability_2D
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
% Initial coding MMP, 11/28/2026

% The LPI of this executive of PIETOOLS 2025 is posed by an executive of
% this release: stop with the call that poses it.
error('PIETOOLS_stability_2D:retired', ...
    ['PIETOOLS_stability_2D is an executive of PIETOOLS 2025. Call ' ...
     '[prog,P,info] = PIETOOLS_stability(PIE,settings,''P''): the ' ...
     'form ''P'' poses the same LPI. The executive takes a PIE in ' ...
     'any number of spatial variables; the solver status is ' ...
     'info.status.']);
end
