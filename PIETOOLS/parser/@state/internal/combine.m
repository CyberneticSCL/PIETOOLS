function [out, varargout] = combine(varargin)
out = state();
% pre allocate for speed?

for i=1:nargin
    out = [out;varargin{i}];
end

[matP,outUnique] = findUnique(out);

for i=1:nargin
    varargout{i} = findPerm(outUnique,varargin{i})*matP;
end
end
% rearrange to put inputs first, outputs, ODEs, PDEs in that order
function [A,obj] = findUnique(obj)
end
function [permMat] = findPerm(objA, objB)
end