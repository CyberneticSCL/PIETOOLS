function d6_run(root,testdir,record,outfile,cmpfile,ccroot)
% D6_RUN(ROOT,TESTDIR,RECORD,OUTFILE,CMPFILE,CCROOT) - scratch harness, 11/27/2026.
% Sets the path to the tree ROOT only (factory path + genpath(ROOT), less
% cleaning/, orphaned_file/ and two research harnesses), checks that the
% names used resolve once, inside ROOT, runs test_call2025_sop(RECORD) from
% TESTDIR and saves R to OUTFILE. CMPFILE (optional): an R of another tree;
% the outputs of the current forms are compared with isequaln. CCROOT
% (optional): another tree; checkcode messages of the three edited files
% are printed for both trees.
root = strrep(root,'/','\');
restoredefaultpath;
addpath(genpath(root));
warning('off','MATLAB:rmpath:DirNotFound');
for f = {'PIETOOLS_demos\cuadmm','PIETOOLS_demos\lowrank_2d_stability','cleaning','orphaned_file'}
    rmpath(genpath(fullfile(root,f{1})));
end
addpath(testdir);
addpath('C:\Program Files\Mosek\11.0\toolbox\r2019b');     % sosprogram needs a solver
names = {'poslpivar','lpivar','mat2opvar','lpiprogram','lpidecvar','parse_opvar_spaces',...
         'poslpivar_direct','test_call2025_sop'};
for k = 1:numel(names)
    w = which(names{k},'-all');     w = w(~contains(w,'Shadowed'));
    ok = numel(w)==1 && (k==numel(names) || startsWith(w{1},root,'IgnoreCase',true));
    fprintf('%-20s %s\n',names{k},strjoin(w,' | '));
    if ~ok,     error('d6_run:resolve','%s resolves to %d files',names{k},numel(w));    end
end
fprintf('opvar: %s\ndopvar: %s\n',which('opvar'),which('dopvar'));
t = tic;
R = test_call2025_sop(record);
fprintf('test time %.1f s\n',toc(t));
save(outfile,'R');
for k = 1:numel(R.old)
    fprintf('%2d  %-46s %-22s %s\n',k,R.old{k},R.errs(k).identifier,shorten(R.errs(k).message));
end
if nargin>=5 && ~isempty(cmpfile)
    Q = load(cmpfile);  Q = Q.R;    nsame = 0;
    for k = 1:numel(R.out)
        same = isequaln(R.out{k},Q.out{k});     nsame = nsame + same;
        if ~same,   fprintf('DIFFERS: current case %d, %s\n',k,R.cur{k});   end
    end
    fprintf('current forms: %d of %d outputs equal to those of %s\n',nsame,numel(R.out),cmpfile);
end
if nargin>=6 && ~isempty(ccroot)
    for f = {'lpi_programming_sopvar\poslpivar.m','lpi_programming_sopvar\lpivar.m',...
             'sopvar\constructors\mat2opvar.m'}
        for r = {strrep(ccroot,'/','\'),root}
            m = checkcode(fullfile(r{1},f{1}));
            fprintf('checkcode %-40s %-50s %d: %s\n',f{1},r{1},numel(m),strjoin({m.message},' | '));
        end
    end
end
end

function s = shorten(m)
s = regexprep(m,'\s+',' ');     if numel(s)>120,    s = [s(1:117) '...'];   end
end
