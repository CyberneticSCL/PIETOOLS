function syntax_check(P, listfile, outfile)
% SYNTAX_CHECK(P,LISTFILE,OUTFILE): runs checkcode on every .m file named in
% LISTFILE (one path per line, relative to the PIETOOLS folder P) and writes
% to OUTFILE one line per file with a parse problem, then a summary line.
% A parse problem: the parse tree of the file (mtree) has an error root, or a
% checkcode message reads like a syntax error (case-insensitive: 'parse
% error', 'invalid syntax', 'unterminated', 'might be missing',
% 'unbalanced'). Code is not run.
L = strtrim(strsplit(fileread(listfile), newline));
L = L(~cellfun(@isempty, L));
fid = fopen(outfile, 'w');
bad = 0;
pat = {'parse error','invalid syntax','unterminated','might be missing','unbalanced'};
for k = 1:numel(L)
    f = fullfile(P, L{k});
    if ~exist(f, 'file'), fprintf(fid, 'MISSING\t%s\n', L{k}); bad = bad + 1; continue; end
    hit = {};
    % The parse tree has an error root when the file does not parse
    t = mtree(fileread(f));
    if strcmp(t.root.kind, 'ERR'), hit{end+1} = ['parse tree: ' t.root.string]; end
    m = checkcode(f, '-id');
    for j = 1:numel(m)
        if any(cellfun(@(p) contains(lower(m(j).message), p), pat))
            hit{end+1} = sprintf('line %d: %s (%s)', m(j).line, m(j).message, m(j).id); %#ok<AGROW>
        end
    end
    if ~isempty(hit)
        bad = bad + 1;
        fprintf(fid, 'PARSE\t%s\t%s\n', L{k}, strjoin(hit, ' | '));
    end
end
fprintf(fid, 'SUMMARY\t%d files checked, %d with parse problems\n', numel(L), bad);
fclose(fid);
fprintf('SYNTAX_CHECK %d files, %d with parse problems\n', numel(L), bad);
end
