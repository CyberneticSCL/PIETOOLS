function out = buildopvar(varargin)
opvar out;
for i=1:2:nargin-1
if ischar(varargin{i})
    switch varargin{i}
        case 'lim'
            lim = varargin{i+1};
        case 'kernel'
            R1R2 = varargin{i+1};
        case 'var'
            tmp = varargin{i+1};
            out.var1 = tmp(1);
            out.var2 = tmp(2);
        case 'dom'
            out.I = unique(varargin{i+1});
        case 'multiplier'
            out.R.R0 = varargin{i+1};
    end
end
end

if exist('lim','var')
if poly2double(lim(2))
    out.R.R2 = R1R2;
end
if poly2double(lim(1))
    out.R.R1 = R1R2;
end
end
end