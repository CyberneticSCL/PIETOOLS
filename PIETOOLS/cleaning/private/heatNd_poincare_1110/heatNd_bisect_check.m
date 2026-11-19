function heatNd_bisect_check(P, tag, logfile)
% HEATND_BISECT_CHECK(P,TAG,LOGFILE) runs test_heatNd_poincare in the export P (its
% PIETOOLS folder) on a clean path; appends the console output and a
% one-line verdict to LOGFILE. Same steps as tp_check (oc_census), plus a
% guard that the family-build routines resolve inside P, and a diary so
% the check table survives an error.
P = strrep(P,'/',filesep);
restoredefaultpath;
old = cd(P);  evalc('pietools_path_update');  cd(tempdir);
warning('off','MATLAB:mpath:privateDirectoriesNotAllowedOnPath');
TP = fullfile(P,'testing_codes','private');
if exist(TP,'dir'),  addpath(genpath(TP));  end
diary(logfile);
fprintf('==== [%s] %s  %s\n', tag, P, char(datetime('now')));
fns = {'test_heatNd_poincare','heatNd_poincare','heatNd_poincare_ops','heatNd_sdp', ...
       'heatNd_path','heatNd_solve','lpi_eq_cdopvar','poscopvar','copquadvar', ...
       'lpiprogram','sossolve'};
for k = 1:numel(fns)
    w = which(fns{k},'-all');   w = w(endsWith(w,'.m'));
    if isempty(w) || ~startsWith(w{1},P)
        fprintf('PATH GUARD FAILED: %s -> %s\n', fns{k}, strjoin(w,' | '));
        diary off;  cd(old);  return
    end
end
t0 = tic;
try
    R = test_heatNd_poincare();
    s = sprintf('PASS, %d checks', numel(R));
catch e
    s = sprintf('ERROR [%s] %s (at %s %d)', e.identifier, e.message, e.stack(1).name, e.stack(1).line);
end
fprintf('[%s] test_heatNd_poincare: %s, %.1f s\n', tag, s, toc(t0));
diary off;
cd(old);
end
