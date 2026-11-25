function test_examples_FI_library()
% TEST_EXAMPLES_FI_LIBRARY runs the 15 examples of
% examples_FI_library_PIETOOLS and checks each bound against the constant
% the library reports for it; it errors, listing every example that fails.
%
% Each example returns an upper bound C on the constant of its inequality,
% the solution of an LPI. Where the sharp constant is known (info.sharp),
% C must lie in [sharp*(1-1e-6), sharp*(1+1e-4)]: below that interval the
% LPI would certify a false inequality beyond the solver tolerance; above
% it the bound has lost accuracy (the largest gap measured on 11/18 and
% 11/22/2026 is 6.3e-5 of the sharp constant, example 2). Where only the
% constant of Barraza Martinez et al. is known (info.paper; examples 5 and
% 6), C must not exceed it.
%
% INPUT
% - none.
% OUTPUT
% - none; one line per example, beginning '## FI'; an error if any fails.
%
% MMP, 11/22/2026.
fails = {};
for k = 1:15
    t = tic;
    [C,info] = examples_FI_library_PIETOOLS(k);
    if ~isnan(info.sharp)
        ok = C >= info.sharp*(1-1e-6) && C <= info.sharp*(1+1e-4);
        ref = sprintf('sharp %.8g',info.sharp);
    elseif ~isnan(info.paper)
        ok = C <= info.paper;
        ref = sprintf('paper %.8g',info.paper);
    else                                % no reference: a finite positive bound
        ok = isfinite(C) && C > 0;
        ref = 'no reference';
    end
    if ok,  verdict = 'pass';  else,  verdict = 'FAIL';  fails{end+1} = sprintf('%d (%s)',k,info.name);  end %#ok<AGROW>
    fprintf('## FI %2d  %s  C %.8g  %s  (%.1f s)  %s\n',k,verdict,C,ref,toc(t),info.name);
end
if ~isempty(fails)
    error('test_examples_FI_library:fail','Failed: %s',strjoin(fails,'; '));
end
fprintf('## FI all 15 examples pass\n');
end
