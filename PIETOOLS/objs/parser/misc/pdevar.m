function obj = pdevar(varargin)
if nargout==0
    for i=1:nargin
        if ischar(varargin{i})
            assignin('caller', varargin{i}, signals(1,'pde'));
        else
            error("Invalid syntax for pdevar()");
        end
    end
else
    if nargin==0
        obj = signals(1,'pde');
        return;
    end
    if nargin>=1
        obj = signals(1,'pde',varargin{:});
    end
end
end