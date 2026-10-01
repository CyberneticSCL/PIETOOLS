function check_reserved_names(vars)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% CHECK_RESERVED_NAMES(VARS) raises an error if a spatial variable name in
% VARS ends in '_int' or '_dum'.
%
% INPUT
% - vars:   cellstr of spatial variable names;
%
% NOTES
% The positive-operator constructors name the integration variable theta_d
% and the input dummy s'_d by appending '_int' and '_dum' to s_d, so the
% suffixes are reserved outright, not only where they would collide: the
% set of admissible names then does not depend on which other variables
% are present.
%
% Cost: one regexp over the registry, one name per spatial direction.
%
% See also PARSE_COPVAR_SPACES, COPQUADVAR, SOPQUADVAR, LPIVAR_CDOPVAR.
%
% For support, contact M. Peet, Arizona State University at mpeet@asu.edu

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PIETOOLS - check_reserved_names
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
% Initial coding MMP, 09/30/2026: one public form of the reserved-suffix
%                  check, written identically in 'copquadvar',
%                  'sopquadvar', 'lpivar_cdopvar' and 'spaces2meta_sop'.
%                  Body and message verbatim.

is_reserved = ~cellfun(@isempty,regexp(vars,'_(int|dum)$','once'));
if any(is_reserved)
    error("Spatial variable names may not end in '_int' or '_dum'; "...
          +"'"+string(vars{find(is_reserved,1)})+"' does.")
end

end
