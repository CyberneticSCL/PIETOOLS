function d = cuadmm_shadow(name)                                            % CC, 09/25/2026
% CUADMM_SHADOW  Temporarily REPLACE a PIETOOLS function with an instrumented
% copy, for the experiments that need one.
%
%   cuadmm_shadow('lpisolve')    % capture stub: records the program, never solves
%   cuadmm_shadow('poslpivar')   % 1-D poslpivar with function-handle multipliers
%   cuadmm_shadow('off')         % remove every active shadow
%
% WHY THE SHADOWS ARE STORED AS .txt.  pietools_path_update adds the whole
% repository with genpath, so a file named lpisolve.m anywhere under PIETOOLS
% would replace the real lpisolve for every user. The copies therefore live in
% shadows/ as <name>.m.txt, which MATLAB never treats as a function, and are
% written out as real .m files only into a temporary folder, only while an
% experiment runs.
%
% A shadow cannot be a private/ function either. The executives call lpisolve
% from THEIR OWN folder, where a private copy is invisible, so the stub would
% never intercept anything. It has to sit ahead of the real one on the path.
%
% THE SHADOW IS ACTIVE UNTIL REMOVED. Anything calling lpisolve or poslpivar
% gets the shadow in the meantime; the lpisolve stub warns on every call because
% it answers with a zero solution and a "clean feasible" status. Every experiment
% that installs a shadow removes it at its end, but one that errors partway does
% not. cuadmm_path REMOVES any active shadow (with a warning) on every call --
% including from inside the harness runners -- so its single-resolution check
% never sees one; never run a runner while a shadow is meant to stay active.
switch lower(name)
    case 'off'
        p = strsplit(path,pathsep);
        sh = contains(p,'cuadmm_shadow_');
        if any(sh), path(strjoin(p(~sh),pathsep)); end
        d = '';
        return
    case {'lpisolve','poslpivar'}
    otherwise
        error('cuadmm_shadow:name','unknown shadow "%s"',name);
end
src = fullfile(fileparts(mfilename('fullpath')),'shadows',[name '.m.txt']);
if ~exist(src,'file'), error('cuadmm_shadow:missing','no shadow source %s',src); end
d = fullfile(tempdir,['cuadmm_shadow_' name]);
if ~exist(d,'dir'), mkdir(d); end
copyfile(src,fullfile(d,[name '.m']),'f');
addpath(d,'-begin');
rehash;
warning('cuadmm_shadow:active', ...
    '%s is SHADOWED by %s until cuadmm_shadow(''off'') or cuadmm_path is called.', ...
    name, d);
end
