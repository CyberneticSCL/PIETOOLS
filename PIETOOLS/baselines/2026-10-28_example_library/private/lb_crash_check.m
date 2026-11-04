function lb_crash_check()
% LB_CRASH_CHECK, 10/28/2026: the crash-recovery branch of lb_run on a
% copy of results.csv (a key started twice is recorded as 'crash', not
% run), then the dry run of the reduced plan. runs.log is restored.
here = fileparts(mfilename('fullpath'));
T = fullfile(here,'tmp_crash');
if ~exist(T,'dir'), mkdir(T); end
copyfile(fullfile(here,'results.csv'),T);
fid = fopen(fullfile(T,'started.txt'),'w');
fprintf(fid,'PDE|50|extreme\nPDE|50|extreme\n');  fclose(fid);
lg = fileread(fullfile(here,'runs.log'));
lb_run({'PDE:50'},{'extreme'},struct('results',fullfile(T,'results.csv')));
fid = fopen(fullfile(here,'runs.log'),'w');  fprintf(fid,'%s',lg);  fclose(fid);
L = splitlines(strtrim(fileread(fullfile(T,'results.csv'))));
fprintf('CHECK last line: %s\n',L{end});
rmdir(T,'s');
lb_full('reduced','',true);
end
