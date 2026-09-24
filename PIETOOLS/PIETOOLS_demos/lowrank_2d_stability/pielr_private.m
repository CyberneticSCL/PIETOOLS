function varargout = pielr_private(fn,varargin)                            % CC, 09/23/2026
% PIELR_PRIVATE  Test hook: call one of this package's private/ functions.
%
%   [a,b,...] = pielr_private('build_stab_1d',PIE,st,false)
%
% WHY THIS EXISTS.  private/ is visible only to .m files sitting in the
% package folder, so a regression script outside it cannot reach the builders,
% the layout extractor or the gate -- and cd'ing into the folder does not help,
% because privacy is resolved against the CALLER's file location, not the
% working directory.  The alternatives are worse: moving the internals out of
% private/ exposes them as API, and copying them into the test folder is what
% produced the forked tests_1d core whose bm_lm2 and bm_resid have since
% drifted away from the shipped ones.
%
% This file is a test hook, not API.  Nothing in the package calls it.

f = str2func(fn);
[varargout{1:nargout}] = f(varargin{:});
end
