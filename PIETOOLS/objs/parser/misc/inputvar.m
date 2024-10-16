function obj = inputvar(varargin)
if nargout==0
    for i=1:nargin
        if ischar(varargin{i})
            assignin('caller', varargin{i}, signals(1,'in'));
        else
            error("Invalid syntax for inputvar()");
        end
    end
else
    if nargin==0
        obj = signals(1,'in');
        return;
    end
    if nargin>=1
        obj = signals(1,'in',varargin{:});
    end
end
end