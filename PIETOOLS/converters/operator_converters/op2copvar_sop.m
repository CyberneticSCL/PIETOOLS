function Pm = op2copvar_sop(Pop)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PM = OP2COPVAR_SOP(POP) returns the container POP as is: a 'copvar'       % MMP, 11/22/2026
% or 'cdopvar', the class of the operators of a PIE.                        % MMP, 11/22/2026
%
% INPUT
% - Pop:    'copvar' or 'cdopvar' object; other classes throw an error.     % MMP, 11/22/2026
%
% OUTPUT
% - Pm:     the same object.                                                % MMP, 11/22/2026

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PIETOOLS - op2copvar_sop
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
% For support, contact M. Peet, Arizona State University at mpeet@asu.edu
% Initial coding MMP, 10/08/2026
% MMP, 11/22/2026: the conversion of the retired classes opvar, dopvar and
%   opvar2d (maintainer, 11/22/2026) removed: a container is returned as
%   is, and any other class errors.

% Containers are returned as is                                             % MMP, 11/22/2026
if isa(Pop,'copvar') || isa(Pop,'cdopvar'),     Pm = Pop;   return,     end
% BEGIN MMP, 11/22/2026: the conversion of opvar and dopvar objects (one
% call) and of opvar2d objects (block by block, with the helper is_zero),
% MMP, 10/08/2026, removed
error('op2copvar_sop:class','POP should be a copvar or a cdopvar.');        % MMP, 11/22/2026
% END MMP, 11/22/2026
end
