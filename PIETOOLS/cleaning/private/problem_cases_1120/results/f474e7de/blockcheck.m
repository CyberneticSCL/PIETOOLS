function blockcheck(which_lib, ids, outfile)
% BLOCKCHECK(WHICH_LIB, IDS, OUTFILE) runs blocks of the DDE or NDS/DDF
% library through PIETOOLS_DDE.m or PIETOOLS_DDF.m as a user who uncomments
% a block does: a copy of the library with the block (its code, the example
% settings st and the forms of 11/19/2026) uncommented one level is put first
% on the path, the script runs in the base workspace, and the settings, the
% form and the verdict of the test are read back. One run per flagged test:
% the other flags are set to 0 after the block, since the script keeps only
% the prog of its last executive.
%
% INPUT
%   which_lib: 'DDE' or 'NDSDDF'
%   ids:       cellstr of block ids (harness copies blocks/<id>.m)
%   outfile:   CSV, appended: lib, id, test, form, settings_match, verdict,
%              feasratio, pinf, dinf, numerr, gamma, ndv, seconds, note
%
% NOTES
%   - PIESIM is shadowed by a stub: PIETOOLS_DDE.m simulates before the
%     LPIs, and the simulation is not under test.
%   - Code lines outside the block are commented out in the copy: each
%     library ships with one block active (DDE: an H2 example at its end;
%     NDS/DDF: DDF03), which would otherwise run after the block.
%   - Verdict as lb_exec: certified when |feasratio-1| <= 0.3 and no pinf,
%     dinf or numerr.
%   - settings_match: the settings the executive received equal the st of
%     the block (isequaln); 'no st' when the block sets none.
% MMP, 11/19/2026.
here = fileparts(mfilename('fullpath'));
switch which_lib
    case 'DDE',     libname = 'examples_DDE_library_PIETOOLS';     script = 'PIETOOLS_DDE';
    case 'NDSDDF',  libname = 'examples_NDSDDF_library_PIETOOLS';  script = 'PIETOOLS_DDF';
end
libpath = which(libname);
L = splitlines(fileread(libpath));
S = regexprep(L,'^(\s*)%','$1','once');             % one comment level removed
if ~exist(outfile,'file')
    fid = fopen(outfile,'w');
    fprintf(fid,'lib,id,test,form,settings_match,verdict,feasratio,pinf,dinf,numerr,gamma,ndv,seconds,note\n');
    fclose(fid);
end
flags = {'stability','stability_dual','Hinf_gain','Hinf_gain_dual','Hinf_control','Hinf_estimator','H2_norm','H2_norm_dual'};
for b = 1:numel(ids)
    id = ids{b};
    B = splitlines(fileread(fullfile(here,'blocks',[id,'.m'])));  B = B(2:end);
    code = strtrim(B(~startsWith(strtrim(B),'%') & ~cellfun(@isempty,strtrim(B))));
    % the block in the library: its code lines in order, prose between them;
    % first one comment level removed (a commented block), then the lines as
    % they are (the active block, whose one-level comments, e.g. commented
    % alternatives, would read as code with a level removed)
    hit = [];
    for T = {S, L}
        T = T{1}; %#ok<FXSET>
        for i = find(strcmp(strtrim(T),code{1}))'
            h = i;  j = i + 1;  ok = true;
            for c = 2:numel(code)
                while j <= numel(L) && ~strcmp(strtrim(T{j}),code{c}) && (startsWith(strtrim(T{j}),'%') || isempty(strtrim(T{j})))
                    j = j + 1;
                end
                if j > numel(L) || ~strcmp(strtrim(T{j}),code{c}),  ok = false;  break,  end
                h(end+1) = j;  j = j + 1; %#ok<AGROW>
            end
            if ok,  hit = h;  break,  end
        end
        if ~isempty(hit),  break,  end
    end
    if isempty(hit),  wline(outfile,which_lib,id,'','','','not found',NaN,NaN,NaN,NaN,NaN,NaN,0,'');  continue,  end
    k = hit(end);
    % while k < numel(L) && contains(L{k+1},'MMP, 11/19/2026'),  k = k + 1;  end  % MMP, 11/20/2026 (was)
    while k < numel(L) && (contains(L{k+1},'MMP, 11/19/2026') || contains(L{k+1},'MMP, 11/20/2026')),  k = k + 1;  end  % MMP, 11/20/2026: worktree test lines of 11/20 join the block
    blk = S(hit(1):k);
    if ~startsWith(strtrim(L{hit(end)}),'%'),  blk = L(hit(1):k);  end   % the active block: as it is
    % the st and forms the block sets
    stl = blk(~startsWith(strtrim(blk),'%') & (startsWith(strtrim(blk),'st ') | startsWith(strtrim(blk),'st.') | startsWith(strtrim(blk),'st=')));
    st_exp = [];
    if ~isempty(stl),  st = [];  eval(strjoin(stl,newline));  st_exp = st;  end
    on = flags(cellfun(@(f) any(~cellfun(@isempty,regexp(strtrim(blk),['^',f,'\s*=\s*1\s*;?\s*$']))), flags));
    for t = 1:numel(on)
        tmp = fullfile(tempdir,['blockcheck_',id]);
        if exist(tmp,'dir'),  rmdir(tmp,'s');  end
        mkdir(tmp);
        off = setdiff(flags,on(t));
        % the block active in the shipped library (DDE: an H2 example at the
        % end; NDS/DDF: DDF03) commented out, as a user who switches
        % examples does
        O = [L(1:hit(1)-1); {''}; L(k+1:end)];
        act = ~startsWith(strtrim(O),'%') & ~cellfun(@isempty,strtrim(O));
        O(act) = strcat({'% '},O(act));
        M = [O(1:hit(1)-1); blk; {strjoin(cellfun(@(f) [f,'=0;'],off,'UniformOutput',false),' ')}; O(hit(1)+1:end)];
        fid = fopen(fullfile(tmp,[libname,'.m']),'w');  fprintf(fid,'%s\n',M{:});  fclose(fid);
        fid = fopen(fullfile(tmp,'PIESIM.m'),'w');
        fprintf(fid,'function solution = PIESIM(varargin)\n%% stub (blockcheck): no simulation\nsolution = [];\nend\n');
        fclose(fid);
        addpath(tmp,'-begin');
        assert(startsWith(which(libname),tmp),'blockcheck: the library copy does not shadow the library');
        t0 = tic;  note = '';
        try
            evalin('base','clear');
            evalin('base',[script,';']);
        catch ME
            note = strrep(ME.message,newline,' ');
        end
        sec = toc(t0);
        rmpath(tmp);  rmdir(tmp,'s');
        % read back
        form = '';  fv = {'stability_form','stability_form','Hinf_gain_form','Hinf_gain_form','','','',''};  % H2: the script fixes 'o' and 'c'
        fi = find(strcmp(flags,on{t}));
        if ~isempty(fv{fi}) && evalin('base',sprintf('exist(''%s'',''var'')',fv{fi}))
            form = evalin('base',fv{fi});  if endsWith(on{t},'_dual'),  form = [form,'d'];  end
        end
        match = 'no st';
        if ~isempty(st_exp) && evalin('base','exist(''settings'',''var'')')
            match = mat2str(isequaln(evalin('base','settings'),st_exp));
        end
        [verdict,fr,pinf,dinf,nerr,gam,ndv] = deal('',NaN,NaN,NaN,NaN,NaN,NaN);
        if evalin('base','exist(''prog'',''var'')')
            prog = evalin('base','prog');
            if isfield(prog,'decvartable'),  ndv = numel(prog.decvartable);  end
            if isfield(prog,'solinfo') && isfield(prog.solinfo,'info')
                I = prog.solinfo.info;
                fr = getf(I,'feasratio');  pinf = getf(I,'pinf');  dinf = getf(I,'dinf');  nerr = getf(I,'numerr');
                if abs(fr-1) <= 0.3 && ~pinf && ~dinf && ~nerr,  verdict = 'certified';
                elseif pinf || dinf || abs(fr+1) <= 0.1,      verdict = 'infeasible';
                else,                                          verdict = 'unknown';  end
            end
        end
        if evalin('base','exist(''gamma'',''var'')'),  g = evalin('base','gamma');  if isnumeric(g) && isscalar(g),  gam = double(g);  end,  end
        wline(outfile,which_lib,id,on{t},form,match,verdict,fr,pinf,dinf,nerr,gam,ndv,sec,note);
    end
end
fprintf('## BLOCKCHECK %s DONE %s\n',which_lib,char(datetime('now')));
end

function v = getf(I,f)
v = NaN;  if isfield(I,f) && ~isempty(I.(f)),  v = double(I.(f));  end
end

function wline(outfile,lib,id,test,form,match,verdict,fr,pinf,dinf,nerr,gam,ndv,sec,note)
fid = fopen(outfile,'a');
fprintf(fid,'%s,%s,%s,%s,%s,%s,%.4g,%g,%g,%g,%.8g,%g,%.1f,"%s"\n',lib,id,test,form,match,verdict,fr,pinf,dinf,nerr,gam,ndv,sec,strrep(note,'"',''''));
fclose(fid);
fprintf('##### %s %s %s: %s (form %s, settings %s)\n',lib,id,test,verdict,form,match);
end
