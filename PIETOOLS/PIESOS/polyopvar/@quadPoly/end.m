function e = end(F,k,n)
% A function for indexing quadPoly objects. This will allow operations such
% as obj(end), obj(:,end), obj(1:end,:), etc.

m = F.dim(1); % find matrix dimensions
p = F.dim(2);

if n == 1  % if obj(end), then use linear indexing
    e = m*p;
else
    if k == 1  % if obj(end,...) then last row
        e = m;
    else       % if obj(..., end) then last column
        e = p;
    end
end
end
