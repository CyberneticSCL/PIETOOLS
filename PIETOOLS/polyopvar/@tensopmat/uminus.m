function C_out = uminus(C_in)
% C_OUT = UMINUS(C_IN) takes a 'tensopmat' object C_IN representing a
% matrix of tensor-PI operator C and returns a 'tensopmat' object C_OUT 
% representing -C.
%
% INPUTS
% - C_in:   'tensopmat' object;
%
% OUTPUTS
% - C_out:  'tensopmat' object representing the negative complement of the
%           input operator;
%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PIETOOLS - uminus
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
% DJ, 04/21/2026: Initial coding

C_out = C_in;
for i=1:numel(C_out.ops)
    % Mutliply each of the tensor-PI operators in the matrix with -1
    if ~isempty(C_in.ops{i})
        C_out.ops{i} = -C_in.ops{i};
    end
end

end