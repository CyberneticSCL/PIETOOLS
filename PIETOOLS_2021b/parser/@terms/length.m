function out = length(obj)
dim = obj.operator.dim;
out = sum(dim(:,1));
end