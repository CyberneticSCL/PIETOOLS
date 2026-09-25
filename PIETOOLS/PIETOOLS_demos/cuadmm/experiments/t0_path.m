% t0_path.m -- establish a clean, SINGLE-resolution PIETOOLS path.
% Nested .claude/worktrees copies and a saved pathdef both silently shadow the
% repo under test, so every function this evaluation depends on is asserted to
% resolve EXACTLY once before anything is measured.

cuadmm_path('reset');   % CLEARS the session path on purpose -- that is this script's job --
                        % then rebuilds it from THIS checkout. Was a hardcoded workstation
                        % REPO path plus a loose '.claude|worktrees' strip.

need = {'poslpivar','lpi_eq','lpiprogram','monomials','sossolve', ...
        'poslpivar_2d','lpi_eq_2d','lpivar_2d','opvar2d'};
nbad = 0;
for i = 1:numel(need)
    w = which(need{i},'-all');
    % a classdef and its methods legitimately report several hits; count only
    % distinct FILES that define the name at top level
    w = unique(w);
    fprintf('T0 resolve %-14s %d\n', need{i}, numel(w));
    if numel(w) ~= 1
        nbad = nbad + 1;
        for j = 1:min(numel(w),4), fprintf('T0    %s\n', w{j}); end
    end
end
fprintf('T0 multi_resolved %d\n', nbad);

fprintf('T0 matlabver %s\n', version);
fprintf('T0 sedumi %s\n', mat2str(~isempty(which('sedumi'))));
fprintf('T0 mosek %s\n',  mat2str(~isempty(which('mosekopt'))));
fprintf('T0DONE\n');
