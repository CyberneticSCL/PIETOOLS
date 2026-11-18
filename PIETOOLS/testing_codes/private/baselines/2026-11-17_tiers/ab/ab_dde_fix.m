function ab_dde_fix(tree)
% AB_DDE_FIX(TREE) the ten 1-D bisection rows of the library (DDE01-08,
% NDS05, NDS08) under the light tier with lpi.var.deg{1} = [4 2 2] (the
% heavy entry) in place of [2 2 2], the harness bisection (AB_BISECT =
% lb_bisect of lb_run.m, reltol 1e-3, 14 solves): the thresholds and the
% last program's size, against the pass rows. TREE is a label only.
here = pwd;   abdir = fileparts(mfilename('fullpath'));
info = lb_setup();  cd(here);  addpath(fullfile(here,'gate'));  addpath(abdir);
setappdata(0,'LB_NDV_CAP',2e6);   setappdata(0,'LB_M_CAP',45000);
S = load(fullfile(here,'catalog.mat'));  C = S.C;
keys = {'DDE:DDE01','DDE:DDE02','DDE:DDE03','DDE:DDE04','DDE:DDE05','DDE:DDE06','DDE:DDE07','DDE:DDE08','NDSDDF:NDS05','NDSDDF:NDS08'};
opts = struct('reltol',1e-3,'maxsolves',14);
fprintf('@@AB git %s\n', info.git);
for k = 1:numel(keys)
    c = C(strcmp({C.key},keys{k}));   c.exec = ['sop:',c.exec];
    S0 = lpisettings_sop('light',0,'','mosek');
    S0.lpi.var = struct('deg',{{[4 2 2],1,1}});
    t0 = tic;
    try
        out = ab_bisect(c, S0, opts);
        fprintf('@@AB %-14s %-10s %-10s result %.8g  n=%d ndv=%d  %.1fs | %s\n', ...
            tree, keys{k}, out.status, out.result, out.n, out.ndec, toc(t0), out.note);
    catch e
        fprintf('@@AB %-14s %-10s ERROR %s: %s\n', tree, keys{k}, e.identifier, strtok(e.message,newline));
    end
end
fprintf('@@AB done %s\n', char(datetime('now')));
end
