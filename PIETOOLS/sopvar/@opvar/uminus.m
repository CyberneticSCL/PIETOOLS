function P = uminus(A)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% P = UMINUS(A) or -A negates an 'opvar' blockwise, (-A)_{ij} = -A_{ij}.
%
% INPUTS
% - A:      'opvar' object;
%
% OUTPUTS
% - P:      'opvar' object equal to -A, on the same grid, spaces and
%           dimensions.
%
% NOTES
% Zero blocks stay [], and the metadata is carried over unchanged, since
% negation changes no space, dimension or basis. Each block is negated in
% one pass over its coefficients.
%
% Cost: O(nnz) over all blocks.
%
% See also MINUS, PLUS, OPVAR.
%
% For support, contact M. Peet, Arizona State University at mpeet@asu.edu

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PIETOOLS - uminus(opvar)
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
% Initial coding MMP, 09/25/2026

% Negate each nonzero block; the metadata is unchanged
Cc = A.C;
for ii = 1:numel(Cc)
    if ~isempty(Cc{ii})
        Cc{ii} = -Cc{ii};
    end
end
P = opvar(Cc,metadata(A));

end
