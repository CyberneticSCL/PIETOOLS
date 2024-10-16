function obj = inputvar(varargin)
if nargout==0
    for i=1:nargin
        if ischar(varargin{i})
            assignin('caller', varargin{i}, signals("in"));
        else
            error("Invalid syntax for inputvar()");
        end
    end
else
    if nargin==0
        obj = signals("in");
        return;
    end
    if nargin>=1
        obj = signals("in",varargin{:});
    end
end
end