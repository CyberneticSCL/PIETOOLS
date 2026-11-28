function namecheck_sop(listfile)
% Resolution of each proposed name on the saved path, then with PIETOOLS added.
names = strtrim(splitlines(fileread(listfile))); names = names(~cellfun(@isempty,names));
report('saved path', names);
addpath(genpath('C:\Users\mpeet\OneDrive\Documents\GitHub\PIETOOLS\PIETOOLS'));
report('saved path + genpath(PIETOOLS)', names);
end
function report(tag, names)
fprintf('== %s\n', tag); n = 0;
for k = 1:numel(names)
    w = which('-all', names{k});
    if ~isempty(w), n = n+1; fprintf('%-32s %d: %s\n', names{k}, numel(w), strjoin(w, ' | ')); end
end
fprintf('(%d of %d names resolve; the others resolve nowhere)\n', n, numel(names));
end
