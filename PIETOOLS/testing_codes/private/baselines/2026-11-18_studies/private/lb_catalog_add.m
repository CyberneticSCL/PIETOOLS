function C = lb_catalog_add(ids)
% C = LB_CATALOG_ADD(IDS) adds the PDE rows with the given lb_meta ids
% (cellstr, e.g. {'53',...,'72'}) to catalog.mat of this folder, built as
% LB_BUILD_CATALOG builds a row: the system built (2-D and 3-D initialized,
% not converted), the executive toggles read from the example file, one row
% per toggle ('LIB:ID' for the first, 'LIB:ID/toggle' after). Rows already
% present under those ids are replaced; every other row is kept as recorded
% 10/28/2026. catalog.csv is not rewritten (its writer is local to
% LB_BUILD_CATALOG). MMP, 11/17/2026, for the tier pass.
info = lb_setup();
M = lb_meta();
M = M(strcmp({M.lib},'PDE') & ismember({M.id},ids));
assert(numel(M)==numel(ids),'lb_catalog_add: %d of %d ids found in lb_meta.',numel(M),numel(ids));
S = load(fullfile(info.here,'catalog.mat'));   C = S.C;
keep = true(1,numel(C));
for k = 1:numel(C)
    id = regexprep(C(k).key(numel(C(k).lib)+2:end),'/.*$','');
    if strcmp(C(k).lib,'PDE') && ismember(id,ids),    keep(k) = false;    end
end
C = C(keep);
[names, execs] = lb_toggles();
fam = @(t) erase(t,'_dual');
for k = 1:numel(M)
    m = M(k);
    base = m;
    base.key = [m.lib,':',m.id];  base.type = '';  base.exec = '';
    base.dim = NaN;  base.size = '';  base.primary = true;
    base.setter_effective = NaN;
    fprintf('\n=== catalog: %s (%s)\n', base.key, m.fn);
    try
        [~, tg, si] = lb_build_system(m, [], 'auto');
        base.dim = si.dim;  base.size = si.size;
    catch ME
        base.exec = 'ERROR';
        base.note = strtrim([base.note,' | catalog build error: ',strrep(ME.message,newline,' ')]);
        C(end+1) = orderfields(base,C(1));   continue                       %#ok<AGROW>
    end
    on = [];
    for j = 1:numel(names)
        if tg.(names{j})==1,  on(end+1) = j;  end                           %#ok<AGROW>
    end
    if isempty(on)
        base.note = strtrim([base.note,' | no executive toggle set']);
        base.run = false;
        C(end+1) = orderfields(base,C(1));   continue                       %#ok<AGROW>
    end
    for j = 1:numel(on)
        row = base;
        row.type = names{on(j)};  row.exec = execs{on(j)};
        if j>1
            row.key = [base.key,'/',names{on(j)}];
            row.primary = false;   row.bisect = false;
            if contains(row.type,'_dual') && ~isempty(m.docdual)
                row.doc = m.docdual;
            elseif ~strcmp(fam(row.type),fam(names{on(1)}))
                row.doc = '';
            end
        end
        C(end+1) = orderfields(row,C(1));                                   %#ok<AGROW>
        fprintf('    row %s: %s, dim %d\n', row.key, row.exec, row.dim);
    end
end
save(fullfile(info.here,'catalog.mat'),'C');
fprintf('catalog.mat: %d rows\n', numel(C));
end
