function C_out = uminus(C_in)
% C_OUT = UMINUS(C_IN) takes a 'tensopvar' object C_IN representing a
% tensor-PI operator C and returns a 'tensopvar' object C_OUT representing 
% -C.
%
% INPUTS
% - C_in:   'tensopvar' object;
%
% OUTPUS
% - C_out:  'tensopvar' object representing the negative complement of the
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
% DJ, 01/26/2026: Initial coding

C_out = C_in;
for kk=1:numel(C_out.ops)
    if isa(C_out.ops{kk},'nopvar') || isa(C_out.ops{kk},'ndopvar')
        % Single PI operator, Top*x
        C_out.ops{kk} = -C_out.ops{kk};
    elseif isa(C_out.ops{kk},'cell')
        % Tensor product of PI operators, (T1*x)*(T2*x)
        % --> multiply only the first factor in each term with -1
        for ii=1:size(C_out.ops{kk},1)
            C_out.ops{kk}{ii,1} = -C_out.ops{kk}{ii,1};
        end
    elseif isa(C_out.ops{kk},'struct') && isfield(C_out.ops{kk},'params')
        % Functional operator
        C_out.ops{kk}.params = -C_out.ops{kk}.params;
    end
end

end