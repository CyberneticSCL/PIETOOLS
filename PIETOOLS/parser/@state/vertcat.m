function obj = vertcat(varargin)
if nargin==1
    obj = varargin{1};
else
    objA = varargin{1};
    objB = varargin{2};
    if ~isa(objA,'state')||~isa(objB,'state')
        error('Only state class objects can be vertically concatenated');
    end
    
    temptype = [objA.type;objB.type];
    tempvecLength = objA.vecLength+objB.vecLength;
    tempdifforder = [objA.diff_order; objB.diff_order];
    tempvar = [objA.var_indep; objB.var_indep];
    
    tempsegregation = [objA.segregation, objA.vecLength+1];
    obj = state(temptype,tempvecLength,tempdifforder,tempvar,tempsegregation);
    if nargin>2
        obj = vertcat(obj,varargin{3:end});
    end
end
end