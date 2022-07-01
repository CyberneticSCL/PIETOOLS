function obj = vertcat(varargin)
if nargin==1
    obj = varargin{1};
else
    objA = varargin{1};
    objB = varargin{2};
    if ~isa(objA,'state')||~isa(objB,'state')
        error('Only state class objects can be vertically concatenated');
    end
    
    obj(1,1) = objA;
    obj(2,1) = objB;
    if nargin>2
        obj = vertcat(obj,varargin{3:end});
    end
end
end