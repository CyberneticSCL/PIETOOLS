function caps = expand_caps(caps,n3,name)
% Expand a scalar or 1 x n3 degree bound into a 1 x n3 array.

caps = reshape(caps,1,[]);
if isscalar(caps)
    caps = repmat(caps,1,n3);
elseif numel(caps)~=n3
    error("The '"+string(name)+"' degree should be a scalar or have one "...
          +"element per spatial variable.")
end

end


