function cc_1b(root)
% checkcode messages and mtree errors of the files of Phase 1b (11/28/2026).
f = [{'executives/utility_functions/getController.m','executives/utility_functions/getObserver.m', ...
      'lpi_programming/lpi_ineq.m','testing_codes/private/executives/tests/test_retired2025.m'}, ...
     cellfun(@(s) ['executives/retired_2025/' s], {dir(fullfile(root,'executives','retired_2025','*.m')).name}, 'UniformOutput', false)];
for k = 1:numel(f)
    p = fullfile(root, f{k});
    e = count(mtfind(mtree(p,'-file'),'Kind','ERR'));
    m = checkcode(p,'-id');
    fprintf('%-62s mtree errors %d; checkcode %d: %s\n', f{k}, e, numel(m), strjoin(unique({m.id}),','));
end
end
