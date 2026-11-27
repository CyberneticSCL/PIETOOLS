function varargout = sedumi(varargin)
% Phase 3 capture shadow, 11/27/2026: logs SHA-256 digests of the SeDuMi
% inputs A, b, c, K, calls the original, and logs those of the outputs x,
% y. Records only.
persistent h
if isempty(h),  h = d3_orig('sedumi');  end
nin = {'A','b','c','K'};    nout = {'x','y'};
n = min(4,nargin);
d3_log('SEDin',nin(1:n),varargin(1:n),'');
[varargout{1:max(nargout,1)}] = h(varargin{:});
m = min(2,numel(varargout));
d3_log('SEDout',nout(1:m),varargout(1:m),'');
end
