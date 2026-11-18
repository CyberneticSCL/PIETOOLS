function C = pf_cfg(names,nv,extra)
% C = PF_CFG(NAMES,NV,EXTRA) a struct array of term-set configurations
% (name, rows, offs, pd) from the PF_TERMS names in NAMES, followed by the
% explicit configurations EXTRA (struct array with the same fields), for
% PF_HEAT and PF_EXEC. SCRATCH, MMP, 11/15/2026.
C = struct('name',{},'rows',{},'offs',{},'pd',{});
for i = 1:numel(names)
    [rows,offs,pd] = pf_terms(names{i},nv,[]);
    C(end+1) = struct('name',names{i},'rows',rows,'offs',offs,'pd',pd);     %#ok<AGROW>
end
if nargin>2 && ~isempty(extra)
    for i = 1:numel(extra),     C(end+1) = extra(i);    end                 %#ok<AGROW>
end
end
