function varargout = cuadmm_private(fn,varargin)                           % CC, 09/25/2026
% CUADMM_PRIVATE  Call one of this package's private/ functions from outside it.
%
%   S = cuadmm_private('sdpshape',prog)
%
% private/ is visible only to .m files in the package folder itself, and
% privacy is resolved against the CALLER'S FILE, not the working directory, so
% the scripts in experiments/ cannot reach the builders or helpers directly.
% Same device as pielr_private in the low-rank package. Copying the helpers
% into experiments/ instead is what that package warns against: its forked
% copy drifted away from the shipped one.
%
% A hook for the experiments, not API.
f = str2func(fn);
[varargout{1:nargout}] = f(varargin{:});
end
