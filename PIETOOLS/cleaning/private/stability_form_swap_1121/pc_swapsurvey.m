function pc_swapsurvey()
% PC_SWAPSURVEY() lists the tests of the example libraries whose form the
% PIE does not admit, which the executives answer with an error since
% 11/21/2026 (before: a substituted form):
%   stability_dual, Hinf_gain_dual ('Pd', 'Qd') on a PIE with K;
%   Hinf_gain 'Q' on a PIE with Tw nonzero;
%   Hinf_gain_dual 'Qd' on a PIE in two or more variables;
%   H2_norm_dual ('c', 'cco') on a PIE with K.
% The form of a test is the one the example names (stability_form,
% Hinf_gain_form, H2_norm_form in its source, uncommented), else the default
% of PIETOOLS_auto_execute_sop ('Q', 'Q', 'o'; a dual appends 'd', or for H2
% turns 'o' into 'c'). Each catalog row builds its PIE (conversion only).
% results_SWAPSURVEY.csv in this folder.
% MMP, 11/21/2026.
here = fileparts(mfilename('fullpath'));
lb_setup();   cd(here);
Cat = load(fullfile(here,'catalog.mat'));  Cat = Cat.C;
of = fullfile(here,'results_SWAPSURVEY.csv');
fid = fopen(of,'w');  fprintf(fid,'key,lib,fn,exec,form,named,nv,Tw0,hasK,new_error,note\n');  fclose(fid);
cache = struct('k',{},'nv',{},'Tw0',{},'hasK',{},'note',{});
for i = 1:numel(Cat)
    c = Cat(i);
    [fs,fg,fh,named] = example_forms(c,here);
    switch c.exec
        case 'PIETOOLS_PIE2PDEstability',       form = fs;
        case 'PIETOOLS_PIE2PDEstability_dual',  form = [fs,'d'];
        case 'PIETOOLS_Hinf_gain',              form = fg;
        case 'PIETOOLS_Hinf_gain_dual',         form = [fg,'d'];
        case 'PIETOOLS_H2_norm_o',              form = fh;
        case 'PIETOOLS_H2_norm_c',              form = ['c',fh(2:end)];
        otherwise,                              continue
    end
    % The PIE of the row, once per example and index
    k = [c.lib,':',c.fn,':',num2str(c.id)];
    j = find(strcmp({cache.k},k),1);
    if isempty(j)
        e = struct('k',k,'nv',NaN,'Tw0',NaN,'hasK',NaN,'note','');
        try
            PIE = exec_pie_sop(lb_build_system(c, [], true));
            e.nv = size(PIE.vars,1);
            e.Tw0 = double(all(PIE.Tw==0));
            e.hasK = double(isa(PIE.K,'copvar') && size(PIE.K,1)>0);
        catch ME
            e.note = strrep(strrep(ME.message,',',';'),newline,' ');
        end
        cache(end+1) = e;  %#ok<AGROW>
        j = numel(cache);
    end
    e = cache(j);
    err = '';
    if ~isnan(e.nv)
        isg = contains(c.exec,'Hinf');
        if any(strcmp(form,{'Pd','Qd','c','cco'})) && e.hasK,   err = [err,' K'];   end
        if isg && strcmp(form,'Q') && ~e.Tw0,                   err = [err,' Tw'];  end
        if isg && strcmp(form,'Qd') && e.nv>=2,                 err = [err,' Qd-ND']; end
    end
    fid = fopen(of,'a');
    fprintf(fid,'%s,%s,%s,%s,%s,%d,%g,%g,%g,%s,"%s"\n',c.key,c.lib,c.fn,c.exec,form,named,e.nv,e.Tw0,e.hasK,strtrim(err),e.note);
    fclose(fid);
end
fprintf('## PC_SWAPSURVEY DONE %s: %d rows\n',char(datetime('now')),numel(Cat));
end

function [fs,fg,fh,named] = example_forms(c,here)
% The forms an example names in its source (uncommented evalin lines), else
% the defaults of PIETOOLS_auto_execute_sop; named counts the forms found.
fs = 'Q';  fg = 'Q';  fh = 'o';  named = 0;
try
    if strcmp(c.lib,'PDE'),  src = fileread(which(c.fn));
    else,                    src = fileread(fullfile(here,'blocks',[c.id,'.m']));
    end
catch
    return
end
L = splitlines(src);
for k = 1:numel(L)
    t = regexp(L{k},'^\s*evalin\(\s*''base''\s*,\s*''(stability_form|Hinf_gain_form|H2_norm_form)\s*=\s*''''(\w+)''''','tokens','once');
    if isempty(t),  continue,  end
    named = named + 1;
    switch t{1}
        case 'stability_form',  fs = t{2};
        case 'Hinf_gain_form',  fg = t{2};
        case 'H2_norm_form',    fh = t{2};
    end
end
end
