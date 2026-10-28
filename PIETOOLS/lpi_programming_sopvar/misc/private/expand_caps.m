function caps = expand_caps(caps,n3,name)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% CAPS = EXPAND_CAPS(CAPS,N3,NAME) returns the degree bound CAPS as a
% 1 x N3 array: a scalar repeated N3 times, an array of N3 entries as a
% row.
%
% INPUT
% - caps:   'double', a scalar or an array of n3 entries;
% - n3:     'double', the number of variables;
% - name:   'char', the name of the degree, printed in the error message.
%
% OUTPUT
% - caps:   'double', 1 x n3.
%
% NOTES
% CAPS with other than 1 or n3 entries is an error.
%
% For support, contact M. Peet, Arizona State University at mpeet@asu.edu

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PIETOOLS - expand_caps
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
% MMP, 10/28/2026: Header to the header protocol (what the function does,
%                inputs, outputs), the licence block added, step comments
%                in the body; code unchanged. No help text deleted. The
%                file had no log: initial coding MP, 08/22/2026, as a local
%                function of 'possopvar' (the file since renamed
%                'sopquadvar'), moved to private/ by MMP, 09/21/2026.

% As a row; a scalar applies to every variable
caps = reshape(caps,1,[]);
if isscalar(caps)
    caps = repmat(caps,1,n3);
elseif numel(caps)~=n3
    error("The '"+string(name)+"' degree should be a scalar or have one "...
          +"element per spatial variable.")
end

end


