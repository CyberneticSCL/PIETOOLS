function s = d3_hash(v)
% S = D3_HASH(V), 11/27/2026: SHA-256 hex digest of the exact contents of
% V: its class and size, then for a sparse matrix the row indices, column
% indices and values of its nonzeros in storage order, for a full numeric,
% logical or char array its elements, for a cell its elements in order,
% for a struct its sorted field names and values. Two arrays have the same
% digest only if they are equal bit for bit (up to hash collisions).
md = java.security.MessageDigest.getInstance('SHA-256');
feed(md,v);
s = sprintf('%02x',typecast(md.digest(),'uint8'));
end

function feed(md,v)
% FEED(MD,V): V into the digest MD, recursively.
md.update(typecast(uint8([class(v),0]),'int8'));
md.update(typecast(double(size(v)),'int8'));
if issparse(v)
    % No nonzero (c = 0 of a feasibility SDP): class and size only. An
    % empty array makes MATLAB pick update(byte) and fail, 11/27/2026
    [i,j,x] = find(v);
    if isempty(i),  return,     end
    md.update(typecast(double(i(:)),'int8'));
    md.update(typecast(double(j(:)),'int8'));
    if islogical(x),    x = double(x);  end
    md.update(typecast(x(:),'int8'));
elseif isnumeric(v) || islogical(v) || ischar(v)
    if isempty(v),  return,     end
    if islogical(v) || ischar(v),   v = double(v);  end
    md.update(typecast(v(:),'int8'));
elseif iscell(v)
    for k = 1:numel(v),     feed(md,v{k});  end
elseif isstruct(v)
    f = sort(fieldnames(v));
    for k = 1:numel(v)
        for j = 1:numel(f)
            md.update(typecast(uint8([f{j},0]),'int8'));
            feed(md,v(k).(f{j}));
        end
    end
else
    % objects and anything else: by their class and size only
end
end
