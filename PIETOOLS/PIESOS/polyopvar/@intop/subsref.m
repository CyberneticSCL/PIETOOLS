function out = subsref(obj,ref)
% OUT = SUBSREF(OBJ,REF) extracts the elements of the 'intop' object OBJ
% specified by REF. This function is automatically called when calling e.g.
% K(i,j) or K.params for an 'intop' object K.
%
% INPUTS
% - obj:    'intop' object of which to extract certain elements;
% - ref:    struct array with fields 'type' and 'subs' indicating which
%           elements/fields of obj to extract;
%
% OUTPUTS
% - out:    object representing the elements of 'obj' specified by 'ref'.
%           If ref.type='.', this will be out=obj.(ref.subs). If
%           ref.type='()', the output will be an intop object including the
%           elements in the indices specified by ref.subs;
%
% NOTES
% - For each of the parameters, the nonzeros in the selected rows and
%   columns are relabeled, and repeated indices copy the nonzeros.
% - Indexing with '()' accepts at most two indices (linear indexing, or
%   row and column indices), and indexing with '{}' is not supported.
%

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
% DJ, 10/15/2026: Initial coding;

if isempty(ref)
    out = obj;
    return
end

switch ref(1).type
    case '.'
        % Just use the built-in function
        out = builtin('subsref',obj,ref);

    case '()'
        % When calling C(ridcs,cidcs), extract elements ridcs, cidcs from
        % the matrix-valued operator;
        [m,n] = size(obj);
        if isscalar(ref(1).subs)
            % Linear indexing
            subs = ref(1).subs{1};
            if islogical(subs)
                if numel(subs)>numel(obj)
                    error("The logical indices contain a true value outside of the array bounds.")
                end
            elseif isa(subs,'double')
                if any(abs(round(subs))~=subs) || min(subs)<=0
                    error("Array indices must be positive integers or logical values.")
                elseif max(subs)>numel(obj)
                    error("Index exceeds the number of array elements.")
                end
            end
            % Determine the selected elements, with the output of the same shape
            % as for a 'polynomial' array: a column if the functional is a column,
            % and of the shape of the indices otherwise
            if ischar(subs) && strcmp(subs,':')
                subs = (1:m*n)';
            elseif islogical(subs)
                subs = find(subs);
            end
            if n==1
                odim = [numel(subs),1];
            else
                odim = size(subs);
            end
            rsel = mod(subs(:)-1,m) + 1;
            csel = ceil(subs(:)/m);
        elseif numel(ref(1).subs)==2
            subs1 = ref(1).subs{1};
            subs2 = ref(1).subs{2};
            if islogical(subs1)
                if numel(subs1)>size(obj,1)
                    error("The logical indices in position 1 contain a true value outside of the array bounds.")
                end
            elseif isa(subs1,'double')
                if any(abs(round(subs1))~=subs1) || min(subs1)<=0
                    error("Array indices must be positive integers or logical values.")
                elseif max(subs1)>size(obj,1)
                    error("Index in position 1 exceeds the number of array elements.")
                end
            end
            if islogical(subs2)
                if numel(subs2)>size(obj,2)
                    error("The logical indices in position 2 contain a true value outside of the array bounds.")
                end
            elseif isa(subs2,'double')
                if any(abs(round(subs2))~=subs2) || min(subs2)<=0
                    error("Array indices must be positive integers or logical values.")
                elseif max(subs2)>size(obj,2)
                    error("Index in position 2 exceeds the number of array elements.")
                end
            end
            % Determine the selected rows and columns
            if ischar(subs1) && strcmp(subs1,':'),  subs1 = 1:m;            end
            if islogical(subs1),                    subs1 = find(subs1);    end
            if ischar(subs2) && strcmp(subs2,':'),  subs2 = 1:n;            end
            if islogical(subs2),                    subs2 = find(subs2);    end
            odim = [numel(subs1),numel(subs2)];
            rsel = repmat(subs1(:),numel(subs2),1);
            csel = repelem(subs2(:),numel(subs1),1);
        else
            error("Number of indices must not exceed 2.")
        end

        % The element (rsel(e),csel(e)) of the functional is element e of the
        % output, with the elements of the output in column-major order. Relabel
        % the nonzeros of the parameters, copying them for each element of the
        % output that they are mapped to
        nd = numel(obj.Zd);
        [~,~,nZ] = Zdims(obj.Z.in);
        src = rsel + (csel-1)*m;                % element of the functional for each output element
        cnt = full(sparse(src,1,1,m*n,1));      % number of output elements for each element of the functional
        [~,eord] = sort(src);                   % output elements ordered by element of the functional
        start = cumsum(cnt) - cnt;
        out = obj;
        for k=1:numel(obj.params)
            P = obj.params{k};
            [i,j,v] = find(P);
            i = i(:);   j = j(:);   v = v(:);   % row vectors if P has one row
            blk = ceil(i/(nd+1));       d = i - (blk-1)*(nd+1);
            jb = ceil(j/nZ);            z = j - (jb-1)*nZ;
            s = blk + (jb-1)*m;
            c = cnt(s);
            tt = repelem((1:numel(s))',c);      tt = tt(:);     % repelem of a scalar returns a row
            off = repelem(cumsum(c)-c,c);       off = off(:);
            e = eord(start(s(tt)) + (1:numel(tt))' - off);
            ro = mod(e-1,odim(1)) + 1;
            co = ceil(e/odim(1));
            out.params{k} = sparse((ro-1)*(nd+1)+d(tt),(co-1)*nZ+z(tt),v(tt),odim(1)*(nd+1),odim(2)*nZ);
        end

    case '{}'
        error("'{}' type indexing of 'intop' objects is not supported.")
end
end
