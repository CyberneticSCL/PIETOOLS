function sop_mcheck(base, new, paths, changed)
% SOP_MCHECK(BASE,NEW,PATHS,CHANGED) - checks 3d and 3e of Phase 1 (../PLAN.md), 11/28/2026;
% adapted from ../../rename_copvar_opvar_1125/phase1/ren_mcheck.m.
% BASE, NEW: the PIETOOLS folder before and after; PATHS: "old<TAB>new" moves; CHANGED: the
% changed .m files, by their path in NEW.
%   3d. on the factory path plus genpath(NEW), for each moved file: on the path (outside
%       testing_codes/ and private folders) its new name resolves exactly once, at its new path,
%       and a changed old name resolves nowhere; off the path the file exists and its name
%       resolves nowhere but to itself (a private function). The same old names resolve on
%       genpath(BASE), each once.
%   3e. no mtree error node and no checkcode syntax error in the changed files.
L = strtrim(splitlines(fileread(paths)));  L = L(~cellfun(@isempty,L));
P = cellfun(@(s) split(s,sprintf('\t'))', L, 'UniformOutput', false);  P = vertcat(P{:});
onpath = @(p) ~startsWith(p,'testing_codes/') && ~contains(p,'/private/');
nbad = 0;
% BASE: each old name on the path resolves once, inside BASE
restoredefaultpath; addpath(genpath(base));
for k = 1:size(P,1)
    if ~onpath(P{k,1}), continue, end
    [~,o] = fileparts(P{k,1});  w = which(o,'-all');
    if numel(w) ~= 1, fprintf('  3d BASE %s: %d resolution(s)\n', o, numel(w)); nbad = nbad + 1; end
end
restoredefaultpath; addpath(genpath(new));
non = 0;  noff = 0;
for k = 1:size(P,1)
    [~,o] = fileparts(P{k,1});  [~,n] = fileparts(P{k,2});
    f = strrep(fullfile(new,P{k,2}),'\','/');
    if ~isfile(f), fprintf('  3d MISSING %s\n', P{k,2}); nbad = nbad + 1; continue, end
    w = strrep(which(n,'-all'),'\','/');
    if onpath(P{k,2})
        non = non + 1;
        if numel(w) ~= 1 || ~strcmp(w{1},f)
            fprintf('  3d NEW %s: %d resolution(s) %s\n', n, numel(w), strjoin(w,' | ')); nbad = nbad + 1;
        end
    else
        noff = noff + 1;
        w = w(~strcmp(w,f));                     % which -all also lists a private function itself
        if ~isempty(w), fprintf('  3d OFFPATH %s resolves: %s\n', n, strjoin(w,' | ')); nbad = nbad + 1; end
    end
    if ~strcmp(o,n)
        w = which(o,'-all');
        if ~isempty(w), fprintf('  3d OLD %s: %s\n', o, strjoin(w,' | ')); nbad = nbad + 1; end
    end
end
fprintf('RESULT 3d: %d moves (%d on the path, %d off it); %d defects\n', size(P,1), non, noff, nbad);
% 3e
f = strtrim(splitlines(fileread(changed)));  f = f(~cellfun(@isempty,f));
nerr = 0;
for k = 1:numel(f)
    p = fullfile(new, f{k});
    e = count(mtfind(mtree(p,'-file'),'Kind','ERR'));
    m = checkcode(p,'-id');  ids = {m.id};
    bad = ismember(ids,{'SYNER','NOPAR','EOFER','MCDUP','ENDCT','ENDCT2'});
    if e > 0 || any(bad)
        nerr = nerr + 1;  fprintf('  3e ERROR %s (mtree %d, checkcode %s)\n', f{k}, e, strjoin(ids(bad),','));
    end
end
fprintf('RESULT 3e: %d files, %d with errors\n', numel(f), nerr);
end
