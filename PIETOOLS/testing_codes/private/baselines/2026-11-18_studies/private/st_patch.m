function st_patch(rows)
% ST_PATCH(ROWS) prints the final settings of rows of the settings rerun
% (cleaning/private/settings_rerun_1121/PLAN.md) as statements relative to
% the tier of each example file, the form a patch of its st block takes.
%
% INPUT
% - rows: n x 4 'cell' {tag, sfx, key, tier}: the ex_search tag and
%         opts.statedir of the run, the catalog key, the tier the example
%         file starts from (lpisettings_sop(tier)).
% OUTPUT
% Prints KEY | TIER | statements from the tier | statements from the start.
% The final settings of a row are those of exs_<tag>_final.csv: the end of
% the descent (st.S) or its smallest direct certificate (st.Sdir), told
% apart by the number of variables; the solver field is left out (the tiers
% set none). The statements from the start (st.S0, the st of the example
% file for an 'examplest' start) are the patch of the file: those from the
% tier omit a field the search returned to its tier value while the file
% sets it (case 22: dD 0 -> 1, the stripped value).
% MMP, 11/21/2026.
here = fileparts(mfilename('fullpath'));
for i = 1:size(rows,1)
    [tag,sfx,key,tier] = rows{i,:};
    F = readtable(fullfile(here,['exs_',tag,'_final.csv']),'Delimiter',',','TextType','char');
    j = find(strcmp(F.key,key),1);
    L0 = load(fullfile(here,['exs_state',sfx],[strrep(strrep(key,':','_'),'/','_'),'.mat']));
    st = L0.st;     S = st.S;
    if isfield(st,'Sdir') && ~isempty(st.Sdir) && st.costdir(1)==F.ndv(j),  S = st.Sdir;  end
    L = st_diff(lpisettings_sop(tier,0,'','mosek'),S,'st');
    L = L(~startsWith(L,'st.sos_opts.solver'));
    L0 = st_diff(st.S0,S,'st');     L0 = L0(~startsWith(L0,'st.sos_opts.solver'));
    fprintf('%s | %s | %d variables | from the tier: %s | from the start: %s\n',key,tier,F.ndv(j), ...
        strjoin(L,' '),strjoin(L0,' '));
end
end

function L = st_diff(A,B,pre)
% MATLAB statements that turn the settings A into B (fields of B only); the
% routine of ex_search.
L = {};
if isstruct(B) && isstruct(A) && isscalar(B) && isscalar(A)
    for f = fieldnames(B)'
        if isfield(A,f{1})
            L = [L, st_diff(A.(f{1}),B.(f{1}),[pre,'.',f{1}])]; %#ok<AGROW>
        else
            L{end+1} = sprintf('%s.%s = %s;',pre,f{1},lit(B.(f{1}))); %#ok<AGROW>
        end
    end
elseif iscell(B) && iscell(A) && numel(A)==numel(B)
    for i = 1:numel(B)
        if ~isequaln(A{i},B{i}),  L{end+1} = sprintf('%s{%d} = %s;',pre,i,lit(B{i}));  end %#ok<AGROW>
    end
elseif ~isequaln(A,B)
    L{end+1} = sprintf('%s = %s;',pre,lit(B));
end
end

function s = lit(v)
if ischar(v),  s = ['''',v,''''];
elseif isnumeric(v) || islogical(v),  s = mat2str(v);
elseif iscell(v),  s = ['{',strjoin(cellfun(@lit,v,'UniformOutput',false),','),'}'];
else,  s = '<struct>';
end
end
