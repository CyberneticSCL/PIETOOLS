function ren_ccmp(base, new, pairs)
% REN_CCMP(BASE,NEW,PAIRS) - part of check 3e of Phase 1, 11/26/2026: the checkcode message
% identifiers of each changed file, before (BASE, old path) and after (NEW, new path), with
% counts. PAIRS: "old<TAB>new" per changed file (ren_check.sh output).
L = strtrim(splitlines(fileread(pairs)));  L = L(~cellfun(@isempty,L));
ndiff = 0;
for k = 1:numel(L)
    p = split(L{k}, sprintf('\t'));
    a = checkcode(fullfile(base,p{1}),'-id');  b = checkcode(fullfile(new,p{2}),'-id');
    ia = sort({a.id});  ib = sort({b.id});
    if ~isequal(ia, ib)
        ndiff = ndiff + 1;
        fprintf('  3e CHECKCODE %s: before %s | after %s\n', p{2}, strjoin(setdiff(ia,ib),','), strjoin(setdiff(ib,ia),','));
        for m = 1:numel(b), if ~any(strcmp(b(m).id, ia)), fprintf('      line %d: %s\n', b(m).line, b(m).message); end, end
    end
end
fprintf('RESULT 3e checkcode: %d files, %d with other messages after\n', numel(L), ndiff);
end
