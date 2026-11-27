function h = d3_orig(name)
% H = D3_ORIG(NAME), 11/27/2026: a handle to the function NAME that the
% path resolves to without this capture folder. A named function handle
% binds the file it resolves to when created (the phase-0 capture of
% PIESIM_time_integrate relies on the same), so the folder is taken off
% the path for the creation and put back first. Errors if the handle still
% points into this folder.
me = fileparts(mfilename('fullpath'));
warning('off','MATLAB:rmpath:DirNotFound');
rmpath(me);
back = onCleanup(@() addpath(me));
h = str2func(name);
f = functions(h);
if isempty(f.file) || startsWith(f.file,me)
    error('d3_orig:self','%s: no original outside the capture folder (%s).',name,f.file);
end
end
