function obj = state(type,varargin)
if nargin==0
    obj = odevar();
    return;
end
if nargin>=1
    if strcmp(type,'ode')
        obj = odevar(varargin{:});
    elseif strcmp(type,'pde')
        obj = pdevar(varargin);
    elseif strcmp(type,'out')
        obj = outputvar(varargin);
    elseif strcmp(type,'in')
        obj = inputvar(varargin);
    else
        error("Invalid input for state() function");
    end
end
end