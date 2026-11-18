% PF_RUN_P2_LIB  Phase 2, single-domain baselines from the example library
% under the row storages: 1-D rows (stability 'P', light) and the 2-D rows
% 33, 34, 35, 37, 39, 40 (stability 'P'), 45, 46, 47 (gain 'P') and 48
% (estimator), light; storages STO (set by the caller; default stock, Rdef,
% Rlight). The catalog of the routing session (psatz_mix/route4/lb).
% SCRATCH, MMP, 11/16/2026.
if ~exist('STO','var'),     STO = {'stock','Rdef','Rlight'};   end
if ~exist('IDS1','var'),    IDS1 = {'1','2','3','5','8','13','17','21'};  end
if ~exist('IDS2','var'),    IDS2 = {'33','34','35','37','39','40','45','46','47','48'};  end
S = pf_setup(false);
H = 'C:/Users/mpeet/AppData/Local/Temp/claude/363c499c6332a503/af90e0a1-bcff-4e8d-9a24-9aa4c52491cc/scratchpad';
LB = [H '/psatz_mix/route4/lb'];    addpath(LB);
Cat = load([LB '/catalog.mat']);  C = Cat.C;
ids = [IDS1, IDS2];
rows = C(strcmp({C.lib},'PDE') & ismember({C.id},ids));
[~,iu] = unique({rows.id},'stable');    rows = rows(iu);
fprintf('@@ P2LIB %d rows: %s\n',numel(rows),strjoin({rows.id},' '));
out = sprintf('%s/results/p2_lib.mat',S);   Out = [];
for i = 1:numel(rows)
    m = rows(i);
    try
        params = {};
        if m.idxarg,   PDE_t = feval(m.fn,str2double(m.id),0,params);
        else,          PDE_t = feval(m.fn,0,params);
        end
        PDE = initialize_PIETOOLS_PDE(PDE_t,true);
        evalc('PIE = initialize(convert(PDE,{''pie'',''silent''}));');
    catch ME
        fprintf('@@ EX PDE%s build ERROR %s\n',m.id,strrep(ME.message,newline,' '));  continue
    end
    switch m.id
        case {'45','46','47'},  tasks = {'hinfP'};
        case '48',              tasks = {'est'};
        otherwise,              tasks = {'stabP'};
    end
    for tk = tasks
        for sto = STO
            r = pf_exec_row(['PDE' m.id],PIE,tk{1},sto{1},'light');
            if isempty(Out),    Out = r;  else,   Out(end+1) = r;   end     %#ok<SAGROW>
            save(out,'Out');
        end
    end
end
disp('## P2LIB DONE');
