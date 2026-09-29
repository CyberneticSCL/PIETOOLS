function [out,varargout] = subsref(obj,s)                                   % MMP, 09/29/2026
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% out = subsref(obj,s) overloads indexing for sdopvar objects, so that
% obj(I,J) returns the sub-operator mapping the input components J to the
% output components I, rather than indexing the MATLAB object array.
%
% INPUT
% obj:  sdopvar class object
% s:    substruct specifying the indexing operation
%
% OUTPUT
% out:  for s(1).type '.', the requested property; for '()', an sdopvar
%       object of dimension [numel(I),numel(J)]
%
% NOTES:
% Rows of the coefficient matrix C_gamma are ordered with the matrix row
% index outer and the ZL monomial index inner, and columns with the matrix
% column index outer and the ZR monomial index inner (sopvar document
% Sec. 4). Since C_gamma is stored vectorized as A + B'*d, selecting rows
% and columns of the operator is a gather of entries of vec(C_gamma), which
% is a gather of rows of A and of COLUMNS of B. The decision variables are
% the rows of B and are never touched, so the cost is independent of their
% number.
%
% The monomial bases are unaffected: a sub-block keeps every monomial, only
% fewer matrix rows and columns. The canonical multiplier form constrains
% which ZR monomials may carry content in a multiplier direction, and that
% is preserved by discarding whole (row,column) blocks, so a sub-block of a
% canonical object is canonical.
%
% For support, contact M. Peet, Arizona State University at mpeet@asu.edu

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Copyright (C)2026 PIETOOLS Team
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
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% If you modify this code, document all changes carefully and include date
% authorship, and a brief description of modifications
%
% Initial coding MMP, 09/06/2026
% MMP, 09/29/2026: '.' reads return every element of a comma list. With the
%                  one output 'out', {P.Zd{:}}, [P.params.A{:}] and
%                  f(P.ZL{:}) silently kept the first element only, and
%                  P.C{i,j}.Zd{:} through a container raised
%                  MATLAB:TooManyOutputs. A '.' read now returns nargout
%                  values, the builtin count set by
%                  'numArgumentsFromSubscript'; the one-value read and '()'
%                  run the same lines as before, so their outputs and cost
%                  are unchanged. The signature was
% function out = subsref(obj,s)                                             % MMP, 09/29/2026 (was)
% MMP, 09/29/2026: A slice followed by a comma-list brace, e.g.
%                  {P(I,J).Zd{:}}, now errors 'sdopvar:parenThenList'. A
%                  '()' read gets one output, so the list silently kept its
%                  first element.

switch s(1).type
    case '.'
        % Every property access, at any depth, is handled by the builtin.
        % Without this branch no property of an sdopvar is readable.
        if nargout~=1                                                       % MMP, 09/29/2026
            % A comma list: nargout is its length (0 for an empty list;     % MMP, 09/29/2026
            % an explicit subsref(P,s) statement also has 0 and returns     % MMP, 09/29/2026
            % the one value as 'ans', as before).                           % MMP, 09/29/2026
            [varargout{1:nargout}] = builtin('subsref',obj,s);              % MMP, 09/29/2026
            if ~isempty(varargout)                                          % MMP, 09/29/2026
                out = varargout{1};  varargout(1) = [];                     % MMP, 09/29/2026
            end                                                             % MMP, 09/29/2026
            return                                                          % MMP, 09/29/2026
        end                                                                 % MMP, 09/29/2026
        out = builtin('subsref',obj,s);
    case '()'
        if numel(s(1).subs)~=2
            error('An sdopvar must be indexed with two subscripts, as P(I,J).')
        end
        % A comma-list brace after the slice cannot be served (one output). % MMP, 09/29/2026
        for k = 2:numel(s)                                                  % MMP, 09/29/2026
            if strcmp(s(k).type,'{}') && ~all(cellfun(@(x) isscalar(x) && (isnumeric(x) || islogical(x)),s(k).subs)) % MMP, 09/29/2026
                error('sdopvar:parenThenList','Slice first, then read the list: Y = P(I,J); Y.field{...}.') % MMP, 09/29/2026
            end                                                             % MMP, 09/29/2026
        end                                                                 % MMP, 09/29/2026
        indr = s(1).subs{1};    % output components
        indc = s(1).subs{2};    % input components

        % Allow indices to be specified as e.g. (i,:)
        if strcmp(indr,':')
            indr = 1:obj.dims(1);
        end
        if strcmp(indc,':')
            indc = 1:obj.dims(2);
        end
        if islogical(indr),     indr = find(indr);      end
        if islogical(indc),     indc = find(indc);      end
        % Tested with 'any' rather than 'max' so that an empty index set
        % passes and gives an operator with a zero dimension, as it does for
        % an sopvar; 'max' of an empty index returns [], which makes the
        % surrounding '||' throw. Integrality is checked because a
        % fractional index survives the range test and then selects a
        % misaligned window of the coefficient matrix.
        if any(indr<1) || any(indr>obj.dims(1)) || any(mod(indr,1))
            error('row index into an sdopvar must be a positive integer no greater than %d',obj.dims(1))
        end
        if any(indc<1) || any(indc>obj.dims(2)) || any(mod(indc,1))
            error('column index into an sdopvar must be a positive integer no greater than %d',obj.dims(2))
        end

        % Number of monomials on each side, and the row count of the
        % unvectorized coefficient. The trailing 1 covers the case of no
        % spatial variables, where the cell is empty and the product is 1.
        NL = prod([cellfun(@numel,obj.ZL),1]);
        NR = prod([cellfun(@numel,obj.ZR),1]);
        nrow = obj.dims(1)*NL;

        % Rows and columns of C_gamma kept, then their positions in
        % vec(C_gamma). 'idx' is ordered column-major over the sub-block,
        % which is the order of the new vec.
        rows = reshape((indr(:).'-1)*NL + (1:NL).',[],1);
        cols = reshape((indc(:).'-1)*NR + (1:NR).',[],1);
        idx  = reshape(rows + (cols.'-1)*nrow,[],1);

        params = obj.params;
        for ii = 1:numel(params.A)
            params.A{ii} = params.A{ii}(idx);
            params.B{ii} = params.B{ii}(:,idx);
        end

        dims = [numel(indr),numel(indc)];
        out = sdopvar(params,obj.vars,obj.Zd,obj.ZL,obj.ZR,obj.dom,dims);

        % Continue any further indexing, e.g. P(1,1).dims. @sopvar/subsref
        % drops these subscripts and returns the sliced object instead.
        if numel(s)>1
            out = subsref(out,s(2:end));
        end
    case '{}'
        error('Not a valid indexing expression')
end

end
