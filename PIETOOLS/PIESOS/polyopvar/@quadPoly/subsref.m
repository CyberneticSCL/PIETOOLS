function out = subsref(obj, S)
%SUBSREF Indexing for quadPoly.
%
% Supports:
%   - obj(i,j)   : returns a quadPoly submatrix with same bases/variables
%   - obj(k)     : linear indexing into the m x n matrix of entries
%                  (column-major); for a vector k, the entries k as a
%                  vector, oriented as in MATLAB (obj(:) is a column)
%   - obj.field / obj.method(...) : delegated to builtin subsref
%
% Not supported:
%   - {} indexing
%
% INPUTS
% - obj:    'quadPoly' object;
% - S:      'struct' array of the indexing operations, as in
%           obj(i,j), obj(k) or obj.field;
%
% OUTPUTS
% - out:    'quadPoly' object of the selected entries, or the value of
%           the field or method.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PIETOOLS - subsref
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
% SS, DJ, 11/24/2026: Initial coding;

switch S(1).type
    case '.'
        % Delegate dot access (properties/methods) to MATLAB
        out = builtin('subsref', obj, S);
        return;

    case '{}'
        error('quadPoly:subsref:braceNotSupported', '{} indexing is not supported for quadPoly.');

    case '()'
        m = obj.dim(1);
        n = obj.dim(2);

        ds = prod(cellfun(@numel, obj.Zs));  % 1 if empty
        dt = prod(cellfun(@numel, obj.Zt));  % 1 if empty

        idx = S(1).subs;

        if numel(idx) == 1
            % Linear indexing into m×n (column-major)
            lin = idx{1};
            is_colon = ischar(lin) && strcmp(lin, ':');
            if is_colon
                lin = (1:(m*n))';
            elseif islogical(lin)
                lin = find(lin);
            end
            [I,J] = ind2sub([m n], lin(:));
            if numel(lin) > 1
                % The entries lin as a vector: a row for a row quadPoly, or
                % for a row lin into a matrix, else a column (always for :)
                as_row = ~is_colon && ((m==1 && n>1) || (n>1 && isrow(lin)) || (m==1 && n==1 && isrow(lin)));
                blocks = cell(numel(lin),1);
                for r = 1:numel(lin)
                    blocks{r} = obj.C((I(r)-1)*ds+(1:ds), (J(r)-1)*dt+(1:dt));
                end
                if as_row
                    out = quadPoly(horzcat(blocks{:}), obj.Zs, obj.Zt, [1 numel(lin)], obj.ns, obj.nt);
                else
                    out = quadPoly(vertcat(blocks{:}), obj.Zs, obj.Zt, [numel(lin) 1], obj.ns, obj.nt);
                end
                if numel(S) > 1
                    out = subsref(out, S(2:end));
                end
                return;
            end

        elseif numel(idx) == 2
            I = idx{1};
            J = idx{2};

            % Expand ':' explicitly
            if ischar(I) && strcmp(I, ':')
                I = 1:m;
            end
            if ischar(J) && strcmp(J, ':')
                J = 1:n;
            end

            % Expand logical indexing
            if islogical(I)
                I = find(I);
            end
            if islogical(J)
                J = find(J);
            end

            % Validate bounds for numeric indices
            if ~isempty(I) && (min(I) < 1 || max(I) > m)
                error('quadPoly:subsref:indexOutOfBounds', 'Row index out of bounds.');
            end
            if ~isempty(J) && (min(J) < 1 || max(J) > n)
                error('quadPoly:subsref:indexOutOfBounds', 'Column index out of bounds.');
            end

            I = I(:);
            J = J(:);
        else
            error('quadPoly:subsref:badSubscripts', 'Use obj(i,j) or obj(k).');
        end

        mi = numel(I);
        nj = numel(J);

        % Rows in C: for each matrix row i, take its ds rows
        % rowIdx is (mi*ds)×1
        if mi == 0
            rowIdx = zeros(0,1);
        else
            rowIdx = bsxfun(@plus, (I.'-1)*ds, (1:ds)'); % ds×mi
            rowIdx = rowIdx(:);
        end

        % Cols in C: for each matrix col j, take its dt cols
        % colIdx is (nj*dt)×1
        if nj == 0
            colIdx = zeros(0,1);
        else
            colIdx = bsxfun(@plus, (J.'-1)*dt, (1:dt)'); % dt×nj
            colIdx = colIdx(:);
        end

        Csub = obj.C(rowIdx, colIdx);

        % Construct result
        out = quadPoly(Csub, obj.Zs, obj.Zt, [mi nj], obj.ns, obj.nt);

        % Apply further subsref, if any
        if numel(S) > 1
            out = subsref(out, S(2:end));
        end
        return;


    otherwise
        error('quadPoly:subsref:unknownType', 'Unknown subsref type: %s', S(1).type);
end
end
