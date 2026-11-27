function ren_mcheck(base, new, changed)
% REN_MCHECK(BASE,NEW,CHANGED) - checks 3d and 3e of Phase 1 (PLAN.md), 11/26/2026.
% BASE, NEW: the PIETOOLS folder before and after the rename; CHANGED: list of the changed .m
% files, by their path in NEW.
%   3d. on the factory path plus genpath(NEW): opvar, dopvar and each new name resolve exactly
%       once, inside NEW; each old name resolves nowhere; the test utilities, which lie in
%       private folders, exist as files and resolve nowhere. Both classes load, and their
%       method lists equal those of copvar and cdopvar on genpath(BASE), up to the
%       constructor name.
%   3e. no mtree error node and no checkcode syntax error in the changed files.
old = {'copvar','cdopvar','poscopvar','lpivar_cdopvar','mat2copvar_sop','eye_copvar_sop', ...
    'zeros_copvar_sop','pi4_copvar_sop','int2copvar','qp_copvar','lpi_eq_cdopvar', ...
    'poscopvar_direct','poscopvar_lift','posmult_cdopvar','copquadvar','lift_copvar', ...
    'cat_copvar_grid','derive_copvar_meta','dpvar_op_copvar','eq_copvar','mat2copvar_grid', ...
    'merge_copvar_registry','split_copvar','verify_copvar_meta','copvar_space_list', ...
    'parse_copvar_spaces','copvar2tensopmat','tensopmat2copvar','poslpivar_sop'};
nw  = {'opvar','dopvar','poslpivar','lpivar','mat2opvar','eye_opvar', ...
    'zeros_opvar','pi4_opvar','int2opvar','qp_opvar','lpi_eq_dopvar', ...
    'poslpivar_direct','poslpivar_lift','posmult_dopvar','opquadvar','lift_opvar', ...
    'cat_opvar_grid','derive_opvar_meta','dpvar_op_opvar','eq_opvar','mat2opvar_grid', ...
    'merge_opvar_registry','split_opvar','verify_opvar_meta','opvar_space_list', ...
    'parse_opvar_spaces','opvar2tensopmat','tensopmat2opvar','poslpivar_legacy'};
priv = {'parser/@pde_struct/private/mtimes_opvar.m', ...
    'testing_codes/private/sopvar/Testfolder/rand_opvar.m', ...
    'testing_codes/private/sopvar/Testfolder/rand_dopvar.m', ...
    'testing_codes/private/sopvar/Testfolder/Test_opvar.m', ...
    'testing_codes/private/sopvar/Testfolder/Test_dopvar.m', ...
    'testing_codes/private/sopvar/Testfolder/Test_qp_opvar_nd.m', ...
    'testing_codes/private/PIETOOLS_demos/sopvar_demos/test_opquadvar_faces.m'};
% method lists on BASE
restoredefaultpath; addpath(genpath(base));
mb = {methods('copvar'), methods('cdopvar')};
nbad = 0;
if ~isempty(which('opvar')) || ~isempty(which('dopvar'))
    fprintf('  3d BASE: opvar or dopvar resolves on BASE\n'); nbad = nbad + 1;
end
setappdata(0,'ren_mcheck',{base,new,changed,mb,nbad,old,nw,priv});
clear classes; %#ok<CLCLS>   % clears this workspace too (sopvar, pde_struct, ... exist in both trees)
c = getappdata(0,'ren_mcheck'); [base,new,changed,mb,nbad,old,nw,priv] = c{:}; %#ok<ASGLU>
restoredefaultpath; addpath(genpath(new));
for k = 1:numel(nw)
    w = which(nw{k},'-all');
    if numel(w) ~= 1 || ~startsWith(strrep(w{1},'\','/'), strrep(new,'\','/'))
        fprintf('  3d NEW %s: %d resolution(s) %s\n', nw{k}, numel(w), strjoin(w,' | ')); nbad = nbad + 1;
    end
    w = which(old{k},'-all');
    if ~isempty(w), fprintf('  3d OLD %s: %s\n', old{k}, strjoin(w,' | ')); nbad = nbad + 1; end
end
for k = 1:numel(priv)
    [~,nm] = fileparts(priv{k});
    if ~isfile(fullfile(new,priv{k})), fprintf('  3d MISSING %s\n', priv{k}); nbad = nbad + 1; end
    w = which(nm,'-all');                    % which -all also lists a class private function
    w = w(~strcmp(strrep(w,'\','/'), strrep(fullfile(new,priv{k}),'\','/')));
    if ~isempty(w), fprintf('  3d PRIVATE %s resolves: %s\n', nm, strjoin(w,' | ')); nbad = nbad + 1; end
end
P = opvar(); D = dopvar(); %#ok<NASGU>
mc = meta.class.fromName('pde_struct'); %#ok<NASGU>   % its InferiorClasses names ?opvar (H4)
mn = {methods('opvar'), methods('dopvar')}; cls = {'copvar','cdopvar'; 'opvar','dopvar'};
for k = 1:2
    a = sort(strrep(mb{k}, cls{1,k}, cls{2,k}));  b = sort(mn{k});
    if ~isequal(a,b)
        fprintf('  3d METHODS %s: %s\n', cls{2,k}, strjoin(setxor(a,b),' ')); nbad = nbad + 1;
    end
end
fprintf('RESULT 3d: %d names, %d private files; class(opvar()) = %s, class(dopvar()) = %s; methods %d, %d; %d defects\n', ...
    numel(nw), numel(priv), class(P), class(D), numel(mn{1}), numel(mn{2}), nbad);
% 3e
f = strtrim(splitlines(fileread(changed)));  f = f(~cellfun(@isempty,f));
nerr = 0;
for k = 1:numel(f)
    p = fullfile(new, f{k});
    e = count(mtfind(mtree(p,'-file'),'Kind','ERR'));
    m = checkcode(p,'-id');  ids = {m.id};
    bad = ismember(ids,{'SYNER','NOPAR','EOFER','MCDUP','ENDCT','ENDCT2'});
    if e > 0 || any(bad)
        nerr = nerr + 1;  fprintf('  3e ERROR %s (mtree %d, checkcode %s)\n', f{k}, e, strjoin(ids(bad),','));
    end
end
fprintf('RESULT 3e: %d files, %d with errors\n', numel(f), nerr);
end
