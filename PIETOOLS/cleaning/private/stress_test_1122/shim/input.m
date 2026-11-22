function r = input(prompt, varargin)
% R = INPUT(PROMPT, ...) answers input() in a batch run, where the
% built-in errors ("Support for user input is required"). Each call takes
% the next answer of the environment variable ST_INPUT, answers separated
% by '|', in order; the last answer repeats. With the 's' option the
% answer is returned as text; without it, the answer is evaluated ('' gives
% []), as the built-in does. The prompt and the answer are printed.
%
% Stress test of 11/22/2026 (inventory_B.md): put this folder first on the
% path with addpath, set ST_INPUT with setenv, then run the library or
% script. MATLAB warns that the name shadows a built-in; that is intended.
persistent k
if isempty(k),  k = 0;  end
k = k + 1;
a = strsplit(getenv('ST_INPUT'),'|','CollapseDelimiters',false);
r = a{min(k,numel(a))};
fprintf('%s[input shim answer: ''%s'']\n',strrep(prompt,'\n',newline),r);
if nargin < 2,  if isempty(r),  r = [];  else,  r = eval(r);  end,  end
end
