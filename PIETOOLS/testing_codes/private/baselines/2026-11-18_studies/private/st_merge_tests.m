function st_merge_tests(what)
% ST_MERGE_TESTS(WHAT) runs the checks of the fast-forward to origin/ndopvar
% 7764389f (11/18/2026) that the harness rows (st_run 'M1L', 'M2L') do not
% cover: the container unit tests that 7d0acf2c (opvar size) updated, plus
% the QP and getsol tests; the 17 examples of examples_FI_library_PIETOOLS
% (12fd1c3e); and DEMO11_poincare_inequality_2D (7764389f). Results are
% printed to the log only, on lines that begin with '## '.
% WHAT: 'all' (default) or 'tests' (the unit tests only).
% MMP, 11/18/2026.
% MMP, 11/23/2026: the root R is the checkout that holds this file (was the
%   live repo: a run from an export or a worktree put the live tree on the
%   path; stress test of 11/22/2026).
if nargin<1,  what = 'all';  end
% R = 'C:/Users/mpeet/OneDrive/Documents/GitHub/PIETOOLS/PIETOOLS';         % MMP, 11/23/2026 (was)
% R: the checkout holding this file, its first ancestor with                % MMP, 11/23/2026
% pietools_path_update.m (the loop ends at the drive root if none);         % MMP, 11/23/2026
% '/' separators as in the literal, for the startsWith check below          % MMP, 11/23/2026
R = fileparts(mfilename('fullpath'));                                       % MMP, 11/23/2026
while ~isfile([R '/pietools_path_update.m']) && ~strcmp(R,fileparts(R))     % MMP, 11/23/2026
    R = fileparts(R);                                                       % MMP, 11/23/2026
end                                                                         % MMP, 11/23/2026
R = strrep(R,'\','/');                                                      % MMP, 11/23/2026
cd(R);  pietools_path_update;
% The helpers of the tests are off the path with the rest of
% testing_codes/private: rand_opvar, rand_dopvar (Testfolder); rand_sopvar,
% rand_sdopvar (Testfolder/sdopvar); pi_poly_grid (sdopvar/claude_tests);
% cx_plant (claude_tests/cx_exec); subs_dvar (lpi_programming/misc)
TF = fullfile(R,'testing_codes','private','sopvar','Testfolder');
addpath(TF,fullfile(TF,'sdopvar'),fullfile(TF,'sdopvar','claude_tests'), ...
        fullfile(TF,'sdopvar','claude_tests','cx_exec'), ...
        fullfile(R,'testing_codes','private','lpi_programming','misc'));
% One resolution of the container class, inside the main checkout
w = strrep(which('opvar'),'\','/');    fprintf('## opvar resolves to %s\n',w);
if ~startsWith(w,R),  error('st_merge_tests:path','opvar outside the main checkout');  end
tests = {'sopvar/Testfolder','Test_opvar';
         'sopvar/Testfolder','Test_dopvar';
         'sopvar/Testfolder','Test_int2sopvar';
         'sopvar/Testfolder','Test_qp_opvar_nd';
         'lpi_programming/tests','test_dpvar_ops';
         'lpi_programming/tests','test_getsol'};
for k = 1:size(tests,1)
    rng(0);     % the Test_* scripts draw random operators without a seed
    run_item(fullfile(R,'testing_codes','private',tests{k,1},[tests{k,2},'.m']),['TEST ',tests{k,2}]);
end
if strcmp(what,'tests'),  fprintf('## ST_MERGE_TESTS DONE %s\n',char(datetime('now')));  return,  end
run_item(fullfile(R,'PIETOOLS_demos','DEMO11_poincare_inequality_2D.m'),'DEMO11');
for k = 1:17
    t0 = tic;
    try
        [C,info] = examples_FI_library_PIETOOLS(k,'plot=0');
        fprintf('## FI %2d  C %.6g  sharp %.6g  paper %.6g  %.1f s  %s\n', ...
                k,C,info.sharp,info.paper,toc(t0),info.name);
    catch ME
        fprintf('## FI %2d  ERROR %s (%.1f s)\n',k,strrep(ME.message,newline,' '),toc(t0));
    end
    close all
end
fprintf('## ST_MERGE_TESTS DONE %s\n',char(datetime('now')));
end

function run_item(f,label)
% runs the script or function file F in this workspace; the label and the
% timer sit in appdata, since a script may begin with 'clear'
setappdata(0,'MT_LABEL',label);     setappdata(0,'MT_T0',tic);
fprintf('## %s begin\n',label);
try
    run(f);
    fprintf('## %s end %.1f s\n',getappdata(0,'MT_LABEL'),toc(getappdata(0,'MT_T0')));
catch ME
    fprintf('## %s ERROR %s (%.1f s)\n',getappdata(0,'MT_LABEL'), ...
            strrep(ME.message,newline,' '),toc(getappdata(0,'MT_T0')));
end
close all
end
