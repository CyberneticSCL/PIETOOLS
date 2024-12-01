function obj = set(obj,varargin)
if floor(nargin/2)==nargin/2
    error("set() function for state class object requires name-value pair format");
end
if length(obj)>=2
    error("set() function cannot be used on state objects obtained through vertical concatenation");
end
for i=1:2:nargin-1
switch varargin{i}
    case 'dim'
        tmpvar(1) = pvar('t');
        for k=2:varargin{i+1}+1
                tmpvar(k) = pvar(['s',num2str(k-1)]);
        end
        obj.prop{3} = tmpvar;
        obj.prop{4} = repmat({[0,1]},length(tmpvar)-1,1);
        obj.prop{5} = inf(length(tmpvar),1);
        obj.prop{6} = zeros(length(tmpvar),1);
    case 'type'
        obj.prop{1} = varargin{i+1};
    case 'len'
        obj.prop{2} = varargin{i+1};
    case 'var'
        obj.prop{3} = varargin{i+1};
    case 'dom'
        obj.prop{4} = varargin{i+1};
    case 'maxdiff'
        obj.prop{5} = varargin{i+1};
    case 'diffOrder'
        obj.prop{6} = varargin{i+1};
    otherwise
        error("Unknown property for state class object");        
end
end
end