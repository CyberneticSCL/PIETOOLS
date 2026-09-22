function deg_list = process_degrees(deg,nblk,n3)
% Expand the degree specification into one struct per basis operator, each
% with fields 'int', 'mult' (1 x n3 arrays) and 'joint' (scalar).

if iscell(deg)
    if numel(deg)~=nblk
        error("A cell 'deg' should have one entry per included basis operator.")
    end
    deg_list = cell(1,nblk);
    for i=1:nblk
        deg_list{i} = process_degrees_one(deg{i},n3);
    end
else
    deg_list = repmat({process_degrees_one(deg,n3)},1,nblk);
end

end


