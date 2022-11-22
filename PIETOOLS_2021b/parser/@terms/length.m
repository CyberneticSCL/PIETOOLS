function out = length(obj)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% terms-class method that specifies number of equations in obj
% Input: 
% obj - terms class objects
% Output:
% out - positive integer len(obj)

dim = obj.operator.dim;
out = sum(dim(:,1));
end