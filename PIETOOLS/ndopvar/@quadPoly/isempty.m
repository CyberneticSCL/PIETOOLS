function logval = isempty(obj)
% almost trivial function. may have its uses.
logval = isempty(obj.C) || isempty(obj.dim);
end