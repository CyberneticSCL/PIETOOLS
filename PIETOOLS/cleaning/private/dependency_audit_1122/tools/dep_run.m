function dep_run(listfile, outfile)
% DEP_RUN(LISTFILE,OUTFILE): direct dependencies ('toponly') of each toolbox file in
% LISTFILE (paths relative to PIETOOLS/), with the toolbox first on the path and
% cleaning/, orphaned_file/, testing_codes/ (subfolders of their private/) after it, so
% that a name resolves into those folders only when the toolbox does not define it.
R = 'C:\Users\mpeet\OneDrive\Documents\GitHub\PIETOOLS\PIETOOLS';
restoredefaultpath;
warning('off','all');
p = strsplit(genpath(R),pathsep);
isz = contains(p,{[filesep 'cleaning'],[filesep 'orphaned_file'],[filesep 'testing_codes'],[filesep '.claude']});
addpath(strjoin(p(~isz & ~cellfun(@isempty,p)),pathsep));
for z = {'cleaning','orphaned_file','testing_codes'}
    addpath(genpath(fullfile(R,z{1})),'-end');
    q = fullfile(R,z{1},'private');
    if exist(q,'dir'), addpath(genpath(q),'-end'); end
end
addpath('C:\Program Files\Mosek\11.0\toolbox\r2019b','-end');
addpath(genpath('C:\YALMIP\yalmip_install\tbxmanager\toolboxes\sedumi\1.3\pcwin64'),'-end');
L = strtrim(splitlines(fileread(listfile)));  L = L(~cellfun(@isempty,L));
fid = fopen(outfile,'w');
t0 = tic;
for i = 1:numel(L)
    f = fullfile(R,strrep(L{i},'/',filesep));
    try
        fl = matlab.codetools.requiredFilesAndProducts(f,'toponly');
        z = fl(contains(fl,{[filesep 'cleaning' filesep],[filesep 'orphaned_file' filesep],[filesep 'testing_codes' filesep]}));
        for k = 1:numel(z), fprintf(fid,'DEP\t%s\t%s\n',L{i},strrep(z{k},[R filesep],'')); end
    catch ME
        fprintf(fid,'ERR\t%s\t%s\n',L{i},regexprep(ME.message,'\s+',' '));
    end
end
fprintf(fid,'DONE\t%d files\t%.1f s\n',numel(L),toc(t0));
fclose(fid);
end
