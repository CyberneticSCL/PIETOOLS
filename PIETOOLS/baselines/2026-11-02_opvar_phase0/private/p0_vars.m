function v = p0_vars(vnames)
% V = P0_VARS(VNAMES) returns the numel(VNAMES) x 1 'polynomial' vector of the
% variables named in the cell VNAMES.
v = cellfun(@(n) polynomial({n}), vnames(:), 'UniformOutput', false);
v = vertcat(v{:});
end
