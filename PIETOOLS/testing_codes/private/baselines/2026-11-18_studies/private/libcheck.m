function libcheck(casefile, outfile)
% LIBCHECK(CASEFILE, OUTFILE) runs each example of the PDE library the way
% examples_PDE_library_PIETOOLS runs it after the answer 'y', without its
% prompts: the example file called as the library calls it (the output st
% collected), the PDE initialized and converted, the example settings put in
% the base workspace as 'settings', and PIETOOLS_auto_execute run on the
% executives the example file flags (with the forms it names,
% stability_form and Hinf_gain_form). It checks the wiring of the
% example-specific settings of 11/19/2026 from the example file to the
% executive.
%
% INPUT
%   casefile: CSV with columns id, file, index (index: the first input of
%             an example function that takes one, else empty)
%   outfile:  CSV written here: case, file, executive, status, value,
%             decision variables, rows m, seconds, settings source, note
%
% NOTES
%   - Flags and forms are reset in the base workspace before each example,
%     as process_inputs of the library resets them.
%   - An example file without the output st leaves 'settings' unset, so
%     the auto-executor would ask for a tier; such cases are run with
%     lpisettings('light') and reported as 'no st'.
% MMP, 11/19/2026.
here = fileparts(mfilename('fullpath'));
T = readtable(casefile,'Delimiter',',','TextType','char');
done = {};
if exist(outfile,'file')
    D = readtable(outfile,'Delimiter',',','TextType','char');  done = D.id;  if isnumeric(done),  done = arrayfun(@num2str,done,'UniformOutput',false);  end
else
    fid = fopen(outfile,'w');
    fprintf(fid,'id,file,executive,status,value,ndv,m,seconds,source,note\n');  fclose(fid);
end
% ids listed in libcheck_redo.txt (one per line) run again after a fix; the
% later lines of the outfile supersede the earlier ones
rf = fullfile(here,'libcheck_redo.txt');
if exist(rf,'file'),  done = setdiff(done,strtrim(splitlines(fileread(rf))));  end
flags = {'stability','stability_dual','Hinf_gain','Hinf_gain_dual','Hinf_estimator', ...
         'Hinf_control','H2_norm','H2_norm_dual','H2_estimator','H2_control'};
outs = struct('stability','info_stability','stability_dual','info_stability_d', ...
    'Hinf_gain','info_Hinf_gain','Hinf_gain_dual','info_Hinf_gain_d', ...
    'Hinf_control','info_control','Hinf_estimator','info_estimator', ...
    'H2_norm','info_H2_norm','H2_norm_dual','info_H2_norm_d', ...
    'H2_control','info_control','H2_estimator','info_estimator');
vals = struct('Hinf_gain','Hinf_gain','Hinf_gain_dual','Hinf_gain_dual', ...
    'Hinf_control','Hinf_gain_control','Hinf_estimator','Hinf_gain_estimator', ...
    'H2_norm','H2_norm','H2_norm_dual','H2_norm_dual', ...
    'H2_control','H2_norm_control','H2_estimator','H2_norm_estimator');
progs = struct('stability','prog_stability','stability_dual','prog_stability_d', ...
    'Hinf_gain','prog_Hinf_gain','Hinf_gain_dual','prog_Hinf_gain_d', ...
    'Hinf_control','prog_control','Hinf_estimator','prog_estimator', ...
    'H2_norm','prog_H2_norm','H2_norm_dual','prog_H2_norm_d', ...
    'H2_control','prog_control','H2_estimator','prog_estimator');
for i = 1:height(T)
    if iscell(T.id),  cs = T.id{i};  else,  cs = num2str(T.id(i));  end
    if any(strcmp(done,cs)),  continue,  end
    fn = T.file{i};  ix = T.index(i);  if iscell(T.index),  ix = str2double(T.index{i});  end
    fprintf('\n##### case %s  %s  %s\n',cs,fn,char(datetime('now')));
    evalin('base','clear');
    evalin('base','silent_initialize_pde = true;');
    for f = flags,  evalin('base',[f{1},' = 0;']);  end
    evalin('base','stability_form = '''';');  evalin('base','Hinf_gain_form = '''';');
    evalin('base','H2_norm_form = '''';');
    t0 = tic;  src = 'st';  note = '';
    try
        nout = nargout(fn);  args = {0,{}};
        if isfinite(ix),  args = [{ix},args];  end %#ok<AGROW>
        o = cell(1,max(nout,1));  [o{:}] = feval(fn,args{:});
        PDE_t = o{1};  st = [];
        if nout >= 2 && isstruct(o{end}) && isfield(o{end},'lpi'),  st = o{end};  end
        if isempty(st),  st = lpisettings('light');  src = 'no st';  end
        PDE = initialize_PIETOOLS_PDE(PDE_t,true);
        PIE = convert_PIETOOLS_PDE(PDE);
        assignin('base','PIE',PIE);  assignin('base','settings',st);
        % the flags the example set, read before the run: the executive
        % stores the gain and the H2 norm under the names of their flags
        on = flags(cellfun(@(f) isequal(evalin('base',f),1), flags));
        evalin('base','PIETOOLS_auto_execute');
        for f = on
            [status, val, ndv, m] = readout(f{1}, outs, vals, progs);
            line(outfile, cs, fn, f{1}, status, val, ndv, m, toc(t0), src, note);
        end
    catch ME
        line(outfile, cs, fn, '', 'error', NaN, NaN, NaN, toc(t0), src, strrep(ME.message,newline,' '));
    end
end
fprintf('## LIBCHECK DONE %s\n',char(datetime('now')));
end

function [status, val, ndv, m] = readout(f, outs, vals, progs)
% the verdict, value and size the auto-executor left in the base workspace
status = '';  val = NaN;  ndv = NaN;  m = NaN;
if isfield(outs,f) && evalin('base',sprintf('exist(''%s'',''var'')',outs.(f)))
    info = evalin('base',outs.(f));
    if isstruct(info) && isfield(info,'status'),  status = info.status;  end
end
if isfield(vals,f) && evalin('base',sprintf('exist(''%s'',''var'')',vals.(f)))
    v = evalin('base',vals.(f));  if isnumeric(v) && isscalar(v),  val = double(v);  end
end
if isfield(progs,f) && evalin('base',sprintf('exist(''%s'',''var'')',progs.(f)))
    prog = evalin('base',progs.(f));
    if isstruct(prog) && isfield(prog,'decvartable'),  ndv = numel(prog.decvartable);  end
    if isstruct(prog) && isfield(prog,'solinfo') && isfield(prog.solinfo,'y'),  m = numel(prog.solinfo.y);  end
end
end

function line(outfile, cs, fn, ex, status, val, ndv, m, sec, src, note)
fid = fopen(outfile,'a');
fprintf(fid,'%s,%s,%s,%s,%.8g,%.0f,%.0f,%.1f,%s,"%s"\n',cs,fn,ex,status,val,ndv,m,sec,src,strrep(note,'"',''''));
fclose(fid);
fprintf('##### case %s  %s: %s %g (ndv %g, %s)\n',cs,ex,status,val,ndv,src);
end
