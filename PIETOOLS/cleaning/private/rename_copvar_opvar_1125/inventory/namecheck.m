% Resolution of each proposed new name: on the saved (startup) path, then with PIETOOLS added.
names = strtrim(splitlines(fileread(fullfile(pwd,'newnames.txt')))); names = names(~cellfun(@isempty,names));
report('saved path', names);
root = 'C:\Users\mpeet\OneDrive\Documents\GitHub\PIETOOLS\PIETOOLS';
addpath(genpath(root));                       % what pietools_path_update does
report('saved path + genpath(PIETOOLS)', names);
function report(tag, names)
  fprintf('== %s (%d path entries)\n', tag, numel(strsplit(path,pathsep)));
  for k = 1:numel(names)
    w = which('-all', names{k});
    if ~isempty(w), fprintf('%-22s %d: %s\n', names{k}, numel(w), strjoin(w, ' | ')); end
  end
  fprintf('(names not listed resolve nowhere)\n');
end
