function r = lb_exec(exec, PIE, S)
% R = LB_EXEC(EXEC, PIE, S) runs one stock executive on a PIE and reads
% the solver exit with the criterion the stock 2-D executives report.
%
% INPUT
%   exec:   name of the executive, e.g. 'PIETOOLS_PIE2PDEstability'; with
%           the prefix 'sop:', its _sop equivalent (table below)
%   PIE:    the PIE
%   S:      settings struct from lpisettings (lpisettings_sop for 'sop:')
%
% OUTPUT
%   r:  struct with fields
%       status  - 'certified', 'inaccurate', 'infeasible', 'unknown', or
%                 'nosolinfo' (no solver exit in prog.solinfo.info)
%       val     - stock solve_val: 1, 0.5, 0, 0 (0 for nosolinfo)
%       gam     - third output of a gain/norm executive (NaN otherwise)
%       ndec    - numel(prog.decvartable)
%       m, nx   - numel(solinfo.y), numel(solinfo.x): SDP rows, cone size
%       fr, pinf, dinf, numerr, cpusec, iter - solver exit (solinfo.info)
%       sec     - wall time of the executive call [s]
%
% NOTES
%   - Criterion of executives/2D/*.m (Conclusion block), feastol 0.3:
%     |feasratio-1|<=0.3, ~numerr, ~pinf, ~dinf          -> certified (1)
%     |feasratio-1|<=0.3, numerr, ~pinf, ~dinf           -> inaccurate (0.5)
%     pinf | dinf | |feasratio+1|<=0.1                   -> infeasible (0)
%     otherwise                                          -> unknown (0)
%     The 1-D executives print no verdict; the same criterion is applied.
%   - For MOSEK, sossolve sets feasratio = MSK_DINF_INTPNT_OPT_STATUS,
%     pinf/dinf from prosta, numerr 2 (prosta UNKNOWN) or 1 (solsta UNKNOWN).
%   - MMP, 11/21/2026: with appdata LB_KEEP_LAST true (sop path), appdata
%     LB_LAST holds the last solve - fn, form, prog, info, PIE, S - for
%     post_verify, which reads the operators of the solve with getsol.

t0 = tic;
gam = NaN;
r.called = exec;
if startsWith(exec,'sop:')
    % sopvar path: the _sop executive and form PIETOOLS_auto_execute_sop
    % uses for this stock executive; S from lpisettings_sop
    M = {'PIETOOLS_PIE2PDEstability',      'PIETOOLS_stability_sop',      {'Q'};
         'PIETOOLS_PIE2PDEstability_dual', 'PIETOOLS_stability_sop',      {'Qd'};
         'PIETOOLS_Hinf_gain',             'PIETOOLS_Hinf_gain_sop',      {'Q'};
         'PIETOOLS_Hinf_gain_dual',        'PIETOOLS_Hinf_gain_sop',      {'Qd'};
         'PIETOOLS_H2_norm_o',             'PIETOOLS_H2_norm_sop',        {'o'};
         'PIETOOLS_H2_norm_c',             'PIETOOLS_H2_norm_sop',        {'c'};
         'PIETOOLS_Hinf_estimator',        'PIETOOLS_Hinf_estimator_sop', {};
         'PIETOOLS_Hinf_control',          'PIETOOLS_Hinf_control_sop',   {};
         'PIETOOLS_H2_estimator',          'PIETOOLS_H2_estimator_sop',   {};
         'PIETOOLS_H2_control',            'PIETOOLS_H2_control_sop',     {}};
    j = find(strcmp(M(:,1),exec(5:end)),1);
    if isempty(j), error('lb_exec:sop','No _sop executive for %s.',exec(5:end)); end
    fn = M{j,2};  args = M{j,3};
    % BEGIN MMP, 11/18/2026: study hooks (2026-11-18_studies, PLAN_1118_studies.md)
    % LB_ST_HOOK: @(S,PIE,exec) -> S, a candidate tier without a tier file.
    if isappdata(0,'LB_ST_HOOK'),  S = feval(getappdata(0,'LB_ST_HOOK'),S,PIE,exec);  end
    % LB_FORM_MAP: struct mapping a form to another, e.g.
    % struct('Q','P','Qd','Pd'); stability rows only, or the executive in
    % the field 'only' (e.g. 'PIETOOLS_Hinf_gain_sop').
    % 'only' may be a cellstr: one map for several executives (Phase B of   % MMP, 11/19/2026
    % the per-example settings, cleaning/private/example_settings_1119)     % MMP, 11/19/2026
    if isappdata(0,'LB_FORM_MAP') && ~isempty(args)
        fm = getappdata(0,'LB_FORM_MAP');
        tgt = 'PIETOOLS_stability_sop';   if isfield(fm,'only'),  tgt = fm.only;  end
%       if strcmp(fn,tgt) && isfield(fm,args{1}),  args = {fm.(args{1})};  end  % MMP, 11/19/2026 (was)
        if any(strcmp(fn,tgt)) && isfield(fm,args{1}),  args = {fm.(args{1})};  end  % MMP, 11/19/2026
    end
    r.called = strjoin([{fn},args],'/');
    % Every output, so the executive verdict (info.status, eta of the last
    % pass) is recorded beside the harness gate; info is the last output.
    outs = cell(1,nargout(fn));
    % LB_KEEP_LAST true: keep the solve for an a posteriori check (ex_search, % MMP, 11/21/2026
    % post_verify, which reads its operators by getsol): the program, info, % MMP, 11/21/2026
    % the form, the PIE and the settings, in appdata LB_LAST                % MMP, 11/21/2026
    keep = isappdata(0,'LB_KEEP_LAST') && isequal(getappdata(0,'LB_KEEP_LAST'),true); % MMP, 11/21/2026
    if isappdata(0,'LB_LAST'),  rmappdata(0,'LB_LAST');  end                % MMP, 11/21/2026
    [outs{:}] = feval(fn, PIE, S, args{:});
    prog = outs{1};     info_ex = outs{end};
    if keep                                                                 % MMP, 11/21/2026
        setappdata(0,'LB_LAST',struct('fn',fn,'form',{args},'prog',prog,'info',info_ex, ...
            'PIE',PIE,'S',S));                                              % MMP, 11/21/2026
    end                                                                     % MMP, 11/21/2026
    if ~strcmp(fn,'PIETOOLS_stability_sop'),  gam = outs{3};  end
    r.exst = '';    r.eta = NaN;    r.passes = NaN;
    if isstruct(info_ex) && isfield(info_ex,'status'),  r.exst = info_ex.status;  end
    if isstruct(info_ex) && isfield(info_ex,'hist') && ~isempty(info_ex.hist)
        r.passes = numel(info_ex.hist);
        if isfield(info_ex.hist,'eta') && ~isempty(info_ex.hist(end).eta),  r.eta = double(info_ex.hist(end).eta);  end
    end
    % END MMP, 11/18/2026
%   if strcmp(fn,'PIETOOLS_stability_sop')                                  % MMP, 11/18/2026 (was)
%       prog = feval(fn, PIE, S, args{:});                                  % MMP, 11/18/2026 (was)
%   else                                                                    % MMP, 11/18/2026 (was)
%       [prog,~,gam] = feval(fn, PIE, S, args{:});                          % MMP, 11/18/2026 (was)
%   end                                                                     % MMP, 11/18/2026 (was)
else
    switch exec
        case {'PIETOOLS_PIE2PDEstability','PIETOOLS_PIE2PDEstability_dual'}
            prog = feval(exec, PIE, S);
        case 'PIETOOLS_well_posedness'
            [prog,~,~,gam] = feval(exec, PIE, S);
        otherwise
            [prog,~,gam] = feval(exec, PIE, S);
    end
end
r.sec = toc(t0);
r.gam = double(gam);
if numel(r.gam)~=1, r.gam = NaN; end
r.ndec = NaN;  r.m = NaN;  r.nx = NaN;
if isfield(prog,'decvartable'), r.ndec = numel(prog.decvartable); end
info = struct('feasratio',NaN,'pinf',NaN,'dinf',NaN,'numerr',NaN,'cpusec',NaN,'iter',NaN);
if isfield(prog,'solinfo') && isfield(prog.solinfo,'info')
    fi = fieldnames(info);
    for k = 1:numel(fi)
        if isfield(prog.solinfo.info,fi{k}), info.(fi{k}) = double(prog.solinfo.info.(fi{k})); end
    end
    if isfield(prog.solinfo,'y'), r.m = numel(prog.solinfo.y); end
    if isfield(prog.solinfo,'x'), r.nx = numel(prog.solinfo.x); end
end
r.fr = info.feasratio;  r.pinf = info.pinf;  r.dinf = info.dinf;
r.numerr = info.numerr;  r.cpusec = info.cpusec;  r.iter = info.iter;
if any(isnan([r.fr r.pinf r.dinf r.numerr]))
    r.status = 'nosolinfo';   r.val = 0;     % solver exit not recorded
elseif abs(r.fr-1)<=0.3 && ~r.numerr && ~r.pinf && ~r.dinf
    r.status = 'certified';   r.val = 1;
elseif abs(r.fr-1)<=0.3 && r.numerr && ~r.pinf && ~r.dinf
    r.status = 'inaccurate';  r.val = 0.5;
elseif r.pinf || r.dinf || abs(r.fr+1)<=0.1
    r.status = 'infeasible';  r.val = 0;
else
    r.status = 'unknown';     r.val = 0;
end
end
