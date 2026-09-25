function P = uminus(A)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% P = UMINUS(A) or -A negates a 'cdopvar' blockwise, (-A)_{ij} = -A_{ij}.
%
% INPUTS
% - A:      'cdopvar' object;
%
% OUTPUTS
% - P:      'cdopvar' object equal to -A, on the same grid, spaces,
%           dimensions and decision variable list.
%
% NOTES
% Negating vec(C) = A + B'*d negates both A and B, so the decision variable
% list is unchanged and every block stays on it; the container invariant
% needs no reconciliation. Zero blocks stay []. 'sdopvar' blocks are negated
% by @sdopvar/uminus and 'sopvar' blocks by @sopvar/uminus.
%
% Cost: O(nnz) over all blocks, including the decision rows of each B.
%
% See also MINUS, PLUS, CDOPVAR.
%
% For support, contact M. Peet, Arizona State University at mpeet@asu.edu

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PIETOOLS - uminus(cdopvar)
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

Cc = A.C;
for ii = 1:numel(Cc)
    if ~isempty(Cc{ii})
        Cc{ii} = -Cc{ii};
    end
end
P = cdopvar(Cc,metadata(A));

end
