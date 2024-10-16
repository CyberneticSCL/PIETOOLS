function obj = outputvar(varargin)
if nargout==0
    for i=1:nargin
        if ischar(varargin{i})
            assignin('caller', varargin{i}, signals(1,'out'));
        else
            error("Invalid syntax for outputvar()");
        end
    end
else
    if nargin==0
        obj = signals(1,'out');
        return;
    end
    if nargin>=1
        obj = signals(1,'out',varargin{:});
    end
end
end