function obj = odevar(varargin)
if nargout==0
    for i=1:nargin
        if ischar(varargin{i})
            assignin('caller', varargin{i}, signals("ode"));
        else
            error("Invalid syntax for odevar()");
        end
    end
else
    if nargin==0
        obj = signals("ode");
        return;
    end
    if nargin>=1
        obj = signals("ode",varargin{:});
    end
end
end