function out = subsref(obj,s)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% terms-class method used in sub-indexing
% Input: 
% obj - terms class objects
% s - sub indexing structure
% Output:
% out - obj.s or obj(s)

switch s(1).type
    case '.'
        out = builtin('subsref',obj,s);
    case '()'
        indc = 1:obj.operator.dim(2,2);
        out = terms(obj.operator(s(1).subs{1},indc),obj.statevec);
    otherwise
        error('Not a valid indexing expression');
end
end