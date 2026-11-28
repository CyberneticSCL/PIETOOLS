function varargout = PIETOOLS_PIE2PDEstability_dual(varargin) %#ok<STOUT> the call only errors
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PIETOOLS_PIE2PDESTABILITY_DUAL - an executive of PIETOOLS 2025, retired.
% A call stops with an error that gives the call of this release:
%
%   PIETOOLS 2025:
%     [prog,P] = PIETOOLS_PIE2PDEstability_dual(PIE,settings)
%   this release:
%     [prog,P,info] = PIETOOLS_stability(PIE,settings,'Qd')
%
% The form 'Qd' poses the LPI of PIETOOLS_PIE2PDEstability_dual. See
% PIETOOLS_stability for the forms and the outputs, and LPISETTINGS for
% SETTINGS.
%
% PIETOOLS - PIETOOLS_PIE2PDEstability_dual
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
error('PIETOOLS_PIE2PDEstability_dual:retired', ...
    ['PIETOOLS_PIE2PDEstability_dual is an executive of PIETOOLS ' ...
     '2025. Call [prog,P,info] = ' ...
     'PIETOOLS_stability(PIE,settings,''Qd''): the form ''Qd'' poses ' ...
     'the same LPI.']);
end
