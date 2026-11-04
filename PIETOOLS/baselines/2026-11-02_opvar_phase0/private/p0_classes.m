function p0_classes(sub)
% P0_CLASSES(SUB) lists the classes of the objects stored in the fixture
% files under SUB, with the number of files that hold each class; for the
% classes opvar, dopvar, opvar2d and dopvar2d also the first 5 file names.
%
% INPUT
% - sub:    char, a folder relative to the fixture root, e.g. 'F4/pass1'.
%
% NOTES
% - Uses the current path. An object whose class is not on the path can
%   load as a 0x0 double without a warning (measured 11/02 with MATLAB's
%   factory path) and is then not listed: put the classes on the path
%   first (p0_path, or addpath(genpath(<tree>))).
% - Walks structs (first element), cells (first 50 entries) and object
%   properties to depth 5.
root = fileparts(fileparts(mfilename('fullpath')));
d = fullfile(root,sub);
f = dir(fullfile(d,'*.mat'));
ws = warning('off','all');
M = containers.Map();
for k = 1:numel(f)
    S = load(fullfile(d,f(k).name));
    c = {};
    fn = fieldnames(S);
    for j = 1:numel(fn),    c = p0_scan(S.(fn{j}), c, 0);    end
    c = unique(c);
    for j = 1:numel(c)
        if M.isKey(c{j}),   M(c{j}) = [M(c{j}), {f(k).name}];
        else,               M(c{j}) = {f(k).name};
        end
    end
end
warning(ws);
ks = M.keys;
fprintf('%s: %d files\n', sub, numel(f));
for j = 1:numel(ks)
    L = M(ks{j});
    fprintf('  %s %d', ks{j}, numel(L));
    if any(strcmp(ks{j},{'opvar','dopvar','opvar2d','dopvar2d'}))
        fprintf(' [%s]', strjoin(L(1:min(5,end)),' '));
    end
    fprintf('\n');
end
end

function c = p0_scan(x, c, depth)
% Appends to cell C the class of every object found in X, to depth 5.
if depth > 5 || isempty(x),  return,    end
if isstruct(x)
    fn = fieldnames(x);
    for j = 1:numel(fn),    c = p0_scan(x(1).(fn{j}), c, depth+1);   end
elseif iscell(x)
    for e = 1:min(numel(x),50),     c = p0_scan(x{e}, c, depth+1);   end
elseif isobject(x) && ~isstring(x)
    c{end+1} = class(x);
    try                                 % properties of a class may be hidden
        p = properties(x);
        for j = 1:numel(p),     c = p0_scan(x(1).(p{j}), c, depth+1);    end
    catch
    end
end
end
