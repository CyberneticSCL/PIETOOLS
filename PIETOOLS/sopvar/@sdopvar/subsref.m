function [out,varargout] = subsref(obj,s)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% [out,varargout] = subsref(obj,s) overloads indexing for sdopvar objects,
% so that obj(I,J) returns the sub-operator mapping the input components J
% to the output components I, rather than indexing the MATLAB object array.
%
% INPUT
% obj:  'sdopvar' object
% s:    substruct specifying the indexing operation
%
% OUTPUT
% out:  for s(1).type '.', the requested property; for '()', an 'sdopvar'
%       object of dimension [numel(I),numel(J)]
% varargout: the remaining elements of a '.' comma list
%
% NOTES:
% Rows of the coefficient matrix C_gamma are ordered with the matrix row
% index outer and the ZL monomial index inner, and columns with the matrix
% column index outer and the ZR monomial index inner. Since C_gamma is
% stored vectorized as A + B'*d, selecting rows and columns of the operator
% is a gather of entries of vec(C_gamma), which is a gather of rows of A
% and of COLUMNS of B. The decision variables are
% the rows of B and are never touched, so the cost is independent of their
% number.
%
% The monomial bases are unaffected: a sub-block keeps every monomial, only
% fewer matrix rows and columns. The canonical multiplier form constrains
% which ZR monomials may carry content in a multiplier direction, and that
% is preserved by discarding whole (row,column) blocks, so a sub-block of a
% canonical object is canonical. P(I) of a one-row or one-column operator
% takes linear indices.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PIETOOLS - subsref (sdopvar)
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
%
% Initial coding MMP, 09/06/2026

switch s(1).type
    case '.'
        % Every property access, at any depth, is handled by the builtin.
        % Without this branch no property of an sdopvar is readable.
        if nargout~=1
            % A comma list: nargout is its length (0 for an empty list;
            % an explicit subsref(P,s) statement also has 0 and returns
            % the one value as 'ans').
            [varargout{1:nargout}] = builtin('subsref',obj,s);
            if ~isempty(varargout)
                out = varargout{1};  varargout(1) = [];
            end
            return
        end
        out = builtin('subsref',obj,s);
    case '()'
        % A comma-list brace after the slice cannot be served (one output).
        for k = 2:numel(s)
            if strcmp(s(k).type,'{}') && ~all(cellfun(@(x) isscalar(x) && (isnumeric(x) || islogical(x)),s(k).subs))
                error('sdopvar:parenThenList','Slice first, then read the list: Y = P(I,J); Y.field{...}.')
            end
        end
        out = op_slice(obj,s(1).subs{:});

        % Continue any further indexing, e.g. P(1,1).dims
        if numel(s)>1
            out = subsref(out,s(2:end));
        end
    case '{}'
        error('Not a valid indexing expression')
end

end
