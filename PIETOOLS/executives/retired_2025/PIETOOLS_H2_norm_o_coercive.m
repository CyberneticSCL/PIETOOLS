function varargout = PIETOOLS_H2_norm_o_coercive(varargin) %#ok<STOUT> the call only errors
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PIETOOLS_H2_NORM_O_COERCIVE - an executive of PIETOOLS 2025, retired. A
% call stops with an error that gives the call of this release:
%
%   PIETOOLS 2025:
%     [prog,Wo,gam] = PIETOOLS_H2_norm_o_coercive(PIE,settings,options)
%   this release:
%     [prog,Wo,gam,~,~,info] = PIETOOLS_H2_norm(PIE,settings,'oco')
%
% The form 'oco' poses the LPI of PIETOOLS_H2_norm_o_coercive. The 2025
% input OPTIONS.h2, a fixed bound, has no counterpart: PIETOOLS_H2_norm
% minimizes the bound. See PIETOOLS_H2_norm for the forms and the outputs,
% and LPISETTINGS for SETTINGS.
%
% PIETOOLS - PIETOOLS_H2_norm_o_coercive
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
error('PIETOOLS_H2_norm_o_coercive:retired', ...
    ['PIETOOLS_H2_norm_o_coercive is an executive of PIETOOLS 2025. ' ...
     'Call [prog,Wo,gam,~,~,info] = ' ...
     'PIETOOLS_H2_norm(PIE,settings,''oco''): the form ''oco'' poses ' ...
     'the same LPI. The 2025 input OPTIONS.h2, a fixed bound, has no ' ...
     'counterpart: PIETOOLS_H2_norm minimizes the bound.']);
end
