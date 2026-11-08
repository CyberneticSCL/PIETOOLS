function out = subsref(obj, S)
%SUBSREF Indexing for quadPoly.
%
% Supports:
%   - obj(i,j)   : returns a quadPoly submatrix with same bases/variables
%   - obj(k)     : linear indexing into the m×n matrix of entries (column-major)
%   - obj.field / obj.method(...) : delegated to builtin subsref
%
% Not supported:
%   - {} indexing

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
            if ischar(lin) && strcmp(lin, ':')
                lin = 1:(m*n);
            end
            [I,J] = ind2sub([m n], lin);
            I = I(:);
            J = J(:);

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
