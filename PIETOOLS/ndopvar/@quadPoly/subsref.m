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
        % Matrix indexing
        m = obj.dim(1);
        n = obj.dim(2);

        ds = prod(cellfun(@numel, obj.Zs));  % =1 if empty
        dt = prod(cellfun(@numel, obj.Zt));  % =1 if empty

        idx = S(1).subs;

        if numel(idx) == 1
            % Linear indexing into m×n (MATLAB column-major)
            lin = idx{1};
            [I,J] = ind2sub([m n], lin);

        elseif numel(idx) == 2
            I = idx{1};
            J = idx{2};

        else
            error('quadPoly:subsref:badSubscripts', 'Use obj(i,j) or obj(k).');
        end

        % Normalize ':' etc. by letting MATLAB do the indexing on 1:m,1:n
        I = (1:m); J = (1:n);  %#ok<NASGU>
        if numel(idx) == 1
            lin = idx{1};
            [I,J] = ind2sub([m n], lin);
        else
            I = idx{1};
            J = idx{2};
            I = (1:m); J = (1:n); %#ok<NASGU>
            I = idx{1}; J = idx{2};
        end

        % Expand ':' and logical indexing safely via MATLAB
        I = (1:m);
        J = (1:n);
        if numel(idx) == 1
            lin = idx{1};
            [I,J] = ind2sub([m n], lin);
        else
            I = (1:m); J = (1:n);
            I = I(idx{1});
            J = J(idx{2});
        end

        mi = numel(I);
        nj = numel(J);

        % Build row/col indices into coefficient matrix C
        % Rows correspond to blocks of size ds for each matrix row
        if mi == 0
            rowIdx = zeros(0,1);
        else
            rowIdx = bsxfun(@plus, (I(:).'-1)*ds, (1:ds)');  % ds × mi
            rowIdx = rowIdx(:);
        end

        % Columns correspond to blocks of size dt for each matrix column
        if nj == 0
            colIdx = zeros(0,1);
        else
            colIdx = bsxfun(@plus, (J(:).'-1)*dt, (1:dt)');  % dt × nj
            colIdx = colIdx(:);
        end

        Csub = obj.C(rowIdx, colIdx);

        out = quadPoly(Csub, obj.Zs, obj.Zt, [mi nj], obj.ns, obj.nt);

        % If there are further subsref operations, apply them
        if numel(S) > 1
            out = subsref(out, S(2:end));
        end
        return;

    otherwise
        error('quadPoly:subsref:unknownType', 'Unknown subsref type: %s', S(1).type);
end
end
