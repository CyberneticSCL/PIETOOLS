function ch3_record(root, file, outfile)
% CH3_RECORD(ROOT,FILE,OUTFILE), 11/27/2026: runs the cells (%%) of FILE
% (the Chapter 3 snippet of PIETOOLS 2025, v5.0) one after another in one
% workspace, with the path set to the tree ROOT only, and records for each
% cell the error it raises, or the class and size of the variables it
% creates; then the display of each operator created. Phase 3 record of
% what the 2025 1-D and 2-D field assignments do after the rename.
root = strrep(root,'/','\');
restoredefaultpath;     addpath(genpath(root));
warning('off','MATLAB:rmpath:DirNotFound');
for f = {'PIETOOLS_demos\cuadmm','PIETOOLS_demos\lowrank_2d_stability','cleaning','orphaned_file'}
    rmpath(genpath(fullfile(root,f{1})));
end
addpath('C:\Program Files\Mosek\11.0\toolbox\r2019b');
txt = fileread(file);
cells = regexp(txt,'(?m)^%%','split');
fid = fopen(outfile,'w');
fprintf(fid,'tree %s; opvar = %s\n',root,which('opvar'));
ws = struct();
for k = 1:numel(cells)
    code = ['%%',cells{k}];     % the split drops the marker of the title line
    if isempty(strtrim(regexprep(code,'(?m)^[ \t]*%[^\n]*$','')))
        continue
    end
    first = strtrim(regexp(code,'^[^\n]*','match','once'));
    before = fieldnames(ws);
    [ws,err] = runcell(code,ws);
    fprintf(fid,'\n== cell %d: %s\n',k,first);
    if ~isempty(err)
        fprintf(fid,'   ERROR: %s\n',regexprep(err,'\s+',' '));
    end
    new = setdiff(fieldnames(ws),before,'stable');
    for j = 1:numel(new)
        v = ws.(new{j});
        if isobject(v)
            fprintf(fid,'   %s: %s, size %s\n',new{j},class(v),mat2str(size(v)));
        end
    end
end
fclose(fid);
type(outfile);
end

function [ws,err] = runcell(code,ws)
% [WS,ERR] = RUNCELL(CODE,WS): evaluates CODE with the variables of WS and
% returns the variables after it, and the error message ('' for none).
% The snippet calls clear: every local of this function may vanish, so the
% message is assigned after evalc returns and ws is rebuilt from who
f = fieldnames(ws);
for k = 1:numel(f),     eval([f{k},' = ws.(f{k});']);   end
clear k f
try
    evalc(code);
    e_rr__ = '';
catch ME_rr__
    e_rr__ = ME_rr__.message;
end
w_rr__ = who;
ws = struct();
for k_rr__ = 1:numel(w_rr__)
    if ~any(strcmp(w_rr__{k_rr__},{'ws','code','e_rr__','ME_rr__','w_rr__','k_rr__'}))
        ws.(w_rr__{k_rr__}) = eval(w_rr__{k_rr__});
    end
end
err = e_rr__;
end
