function pc_sanity()
% PC_SANITY() checks the worktree run of 11/20/2026 before the queues start:
% the path resolves to this checkout, the statement split of ex_search keeps
% column vectors, the stock 2-D storage levels exist, and one job runs.
% MMP, 11/20/2026.
here = fileparts(mfilename('fullpath'));
lb_setup();   cd(here);
root = fileparts(fileparts(fileparts(fileparts(fileparts(here)))));
fns = {'lpisettings_sop','ex_search','diag_amt','pc_run','lb_exec','sossolve', ...
       'PIETOOLS_stability_sop','PIETOOLS_H2_norm_sop','settings_2d_stock','get_lift_degs', ...
       'lpi_ineq_sop','poscopvar','lpisolve'};
bad = 0;
for k = 1:numel(fns)
    w = which(fns{k},'-all');   w = w(~contains(w,'Shadowed'));
    ok = ~isempty(w) && startsWith(lower(strrep(w{1},'\','/')),lower(strrep(root,'\','/')));
    bad = bad + ~ok;
    fprintf('## which %-24s %d %s %s\n',fns{k},numel(w),repmat('OK',1,ok),w{1});
end
fprintf('## mosekopt %s\n## sedumi %s\n',which('mosekopt'),which('sedumi'));
mods = 'st.lpi.pos.len_sos{1} = [2;1;2]; st.lpi.ineq.dD = 0; st.x = {''a'';''b''}';
parts = strtrim(regexp(mods,';(?![^\[\{\(]*[\]\}\)])','split'));
fprintf('## split: %d parts: %s\n',numel(parts),strjoin(parts,' | '));
if numel(parts) ~= 3,  bad = bad + 1;  end
for lv = {'light','heavy'}
    try,  s2 = settings_2d_stock(lv{1});  fprintf('## settings_2d_stock(%s) ok (%s)\n',lv{1},class(s2));
    catch ME,  fprintf('## settings_2d_stock(%s) ERROR %s\n',lv{1},ME.message);  bad = bad + 1;
    end
end
fprintf('## SANITY checks: %d failed\n',bad);
if bad,  error('pc_sanity:failed','%d checks failed',bad);  end
pc_run('SANITY');
end
