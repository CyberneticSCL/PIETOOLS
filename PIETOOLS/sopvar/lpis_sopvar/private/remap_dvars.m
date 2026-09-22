function Bnew = remap_dvars(B,dvars_old,dmap,ndec)
% Move the rows of B from the ordering dvars_old to the global ordering that
% dmap indexes, inserting zero rows for decision variables absent from B.
% The lookup is a hash per name rather than a search over the whole global
% list, so the cost is in the number of names B actually carries.

loc = zeros(numel(dvars_old),1);
for i = 1:numel(dvars_old)
    if ~isKey(dmap,dvars_old{i})
        error("Internal error: unrecognized decision variable.")
    end
    loc(i) = dmap(dvars_old{i});
end

Bnew = sparse(ndec,size(B,2));
if ~isempty(loc)
    Bnew(loc,:) = B;
end

end


