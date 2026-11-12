function R = bench_lpi_degree(nv,dlist,series,problem)
% R = BENCH_LPI_DEGREE(NV,DLIST,SERIES[,PROBLEM]) times the parsing and the
% solving of an LPI as a function of the degree d of its Lyapunov (storage)
% variables, with the legacy classes and with the container classes.
% Parsing time = wall time of the executive minus the SDP solver time
% reported by MOSEK; with the stub MOSEKOPT of stub_mosek/ first on the path
% the solver returns at once and the run measures parsing only.
% PROBLEM 'stability' (default): the direct stability LPI, P > 0 and
% T'PA + A'PT < 0, for the heat equation
%     x_t = sum_{i<=NV} x_{s_i s_i} + lam x   on [0,1]^NV, x = 0 on the boundary
% (lam = pi^2/2 for NV = 1 and pi^2 for NV = 2, both stable).
% PROBLEM 'of': the H-infinity dynamic output feedback LPI of
% PIETOOLS_HINF_OFCONTROL_SOP for
%     x_t = sum_{i<=NV} x_{s_i s_i} + lam x + w1 + u,  x = 0 on the boundary,
%     z = [int x; u],  y = int x + w2
% (lam = 12 for NV = 1, as in DEMO12, and 24 for NV = 2; both unstable).
% INPUT
% - nv:      1 or 2, the number of spatial variables;
% - dlist:   vector of positive integers, the Lyapunov degrees d;
% - series:  cellstr, any of
%             'legacy'    LEGACY_PDESTABILITY / LEGACY_PDESTABILITY_2D (the
%                         stock executive) on PIE2OPVAR(PIE): opvar, dopvar
%                         (opvar2d, dopvar2d) and lpi_programming; stability
%                         only;
%             'sop_stock' the container executive with settings.sop.slack =
%                         'stock' (for stability in 1-D, the stock SDP);
%             'sop'       the container executive with its default slack;
% - problem: (optional) 'stability' (default) or 'of'.
% OUTPUT
% - R: struct array, one row per (d, series), fields problem, nv, d, series,
%      parse_only (true when MOSEKOPT is the stub), t_total, t_solve,
%      t_parse (median seconds over the repetitions), t_build (the LPI
%      construction) and t_lpisolve (the wall time of LPISOLVE, which
%      contains the MOSEK call), both from info.hist of the container
%      executives (NaN for legacy), nrep, parse_reps, build_reps, status,
%      ndv (decision variables), m (equality rows), Ks (PSD block sizes, as
%      text), sumKs2 (sum of squared PSD sizes), err.
% NOTES
% Settings, one family per dimension, d the Lyapunov degree:
%   1-D: sosineq_on = 0, dd1 = dd12 = {d,[d d 2d],[d d 2d]},
%        dd2 = {d+1,[d+2 d+2 2d+2],..}, dd3 = {d,[d+1 d+1 2d+1],..},
%        ddZ = [2d 2d 2d]; equal to lpisettings('veryheavy') at d = 3.
%   2-D: lpisettings('light') with every degree in settings_2d.LF_deg
%        multiplied by d (light at d = 1); the negativity by lpi_ineq.
% Container options: one pass (max_raise = max_raise_feas = 0), no witness,
% no dual read-back, not verbose. A run shorter than 10 s is repeated 3
% times; the medians are reported.

% Initial coding MMP, 11/06/2026: parsing-time baseline of the stability
% LPI against the Lyapunov degree (README.md of baselines/2026-11-05_class_ops).
% MMP, 11/06/2026: PROBLEM 'of', the output feedback LPI (Fig. 5 of the
%                README); PARSE_ONLY flags runs on the stub MOSEKOPT, whose
%                solve is too long to run; t_build and t_lpisolve from
%                info.hist of the container executives, since the paper's
%                parsing time (wall minus solver) also counts the work after
%                the solver call, which differs between stub and full runs.
% MMP, 11/12/2026: The container series get the 1-D settings in the native
%                form (local NATIVE_FROM_STOCK: lpi.pos, lpi.var, sop.stock),
%                since the container executives refuse the stock-format
%                fields (phase 2 of the knob reduction); the legacy series
%                keeps them. Same degrees, so the same programs.

% The PDE and its PIE (copvar operators), and the PIE with opvar operators
if nargin<4 || isempty(problem),    problem = 'stability';  end             % MMP, 11/06/2026
pvar s1 s2
% if nv==1                                                                  % MMP, 11/06/2026 (was)
if strcmp(problem,'stability') && nv==1                                     % MMP, 11/06/2026
    x = pde_var(1,s1,[0,1]);
    PDE = [diff(x,'t')==diff(x,s1,2)+pi^2/2*x; subs(x,s1,0)==0; subs(x,s1,1)==0];
% else                                                                      % MMP, 11/06/2026 (was)
elseif strcmp(problem,'stability')                                          % MMP, 11/06/2026
    x = pde_var(1,[s1;s2],[0,1;0,1]);
    PDE = [diff(x,'t')==diff(x,s1,2)+diff(x,s2,2)+pi^2*x;
           subs(x,s1,0)==0; subs(x,s1,1)==0; subs(x,s2,0)==0; subs(x,s2,1)==0];
% MMP, 11/06/2026: start of the output feedback plants (new)
elseif nv==1    % 'of': the plant of DEMO12
    x = pde_var(1,s1,[0,1]);    w1 = pde_var('in');     w2 = pde_var('in');
    u = pde_var('control');     z = pde_var('out',2);   y = pde_var('sense');
    PDE = [diff(x,'t')==12*x+diff(x,s1,2)+w1+u; subs(x,s1,0)==0; subs(x,s1,1)==0;
           z==[int(x,s1,[0,1]); u];     y==int(x,s1,[0,1])+w2];
else            % 'of' in 2-D
    x = pde_var(1,[s1;s2],[0,1;0,1]);   w1 = pde_var('in');     w2 = pde_var('in');
    u = pde_var('control');     z = pde_var('out',2);   y = pde_var('sense');
    PDE = [diff(x,'t')==diff(x,s1,2)+diff(x,s2,2)+24*x+w1+u;
           subs(x,s1,0)==0; subs(x,s1,1)==0; subs(x,s2,0)==0; subs(x,s2,1)==0;
           z==[int(x,[s1;s2],[0,1;0,1]); u];   y==int(x,[s1;s2],[0,1;0,1])+w2];
% MMP, 11/06/2026: end of the output feedback plants
end
PDE = initialize(PDE);                                                      % MMP, 11/06/2026
PIE = convert(PDE,'pie');
% PIEo = pie2opvar(PIE);                                                    % MMP, 11/06/2026 (was)
PIEo = [];  if any(strcmp(series,'legacy')),    PIEo = pie2opvar(PIE);  end % MMP, 11/06/2026: only the legacy series reads it
parse_only = contains(which('mosekopt'),'stub_mosek');                      % MMP, 11/06/2026: the stub returns at once

R = [];
for d = dlist
    st = settings_for(nv,d);
    if strcmp(problem,'of')     % margins of DEMO12: the free Q, U rescale R, V % MMP, 11/06/2026
        st.eppos = 1e-7;    st.eppos2 = 1e-7;                               % MMP, 11/06/2026
    end                                                                     % MMP, 11/06/2026
    for k = 1:numel(series)
        ser = series{k};
        stk = st;
        stk.sop = struct('max_raise',0,'max_raise_feas',0,'witness',false,'keepdual',false,'verbose',false);
        if strcmp(ser,'sop_stock'),     stk.sop.slack = 'stock';    end
        % The container executives read the native settings since           % MMP, 11/12/2026
        % 11/12/2026: in 1-D the same degrees as lpi.pos, lpi.var, sop.stock % MMP, 11/12/2026
        if nv==1 && ~strcmp(ser,'legacy'),  stk = native_from_stock(stk);  end % MMP, 11/12/2026
        if strcmp(problem,'of') && strcmp(ser,'legacy')                     % MMP, 11/06/2026
            error('bench_lpi_degree:noLegacyOF','No legacy implementation of the output feedback LPI exists.'); % MMP, 11/06/2026
        end                                                                 % MMP, 11/06/2026
        switch ser      % CALL_EXEC returns prog and, for a container executive, info % MMP, 11/06/2026
            case 'legacy'
%               if nv==1,   fh = @() legacy_PDEstability(PIEo,stk);         % MMP, 11/06/2026 (was)
%               else,       fh = @() legacy_PDEstability_2D(PIEo,stk);      % MMP, 11/06/2026 (was)
                if nv==1,   fh = @() call_exec(@legacy_PDEstability,1,PIEo,stk); % MMP, 11/06/2026
                else,       fh = @() call_exec(@legacy_PDEstability_2D,1,PIEo,stk); % MMP, 11/06/2026
                end
            otherwise
%               if strcmp(problem,'of'),    fh = @() PIETOOLS_Hinf_OFcontrol_sop(PIE,stk); % MMP, 11/06/2026 (was)
%               else,                       fh = @() PIETOOLS_stability_sop(PIE,stk,'P'); % MMP, 11/06/2026 (was)
                if strcmp(problem,'of'),    fh = @() call_exec(@PIETOOLS_Hinf_OFcontrol_sop,4,PIE,stk); % MMP, 11/06/2026
                else,                       fh = @() call_exec(@PIETOOLS_stability_sop,3,PIE,stk,'P'); % MMP, 11/06/2026
                end
        end
        row = run_one(fh);
        row.problem = problem;  row.parse_only = parse_only;                % MMP, 11/06/2026
        row.nv = nv;    row.d = d;      row.series = ser;
        R = [R; orderfields(row)]; %#ok<AGROW>
%       fprintf('##nv %d d %d %-9s total %8.2f s  solve %8.2f s  parse %8.2f s  %-10s ndv %8d  m %7d  %s %s\n', ...
%           nv,d,ser,row.t_total,row.t_solve,row.t_parse,row.status,row.ndv,row.m,datestr(now,'HH:MM:SS'),row.err); % MMP, 11/06/2026 (was)
        fprintf('##%s%s nv %d d %d %-9s total %8.2f s  solve %8.2f s  parse %8.2f s  build %8.2f s  %-10s ndv %8d  m %7d  %s %s\n', ...
            problem,repmat(' (parse only)',1,parse_only),nv,d,ser,row.t_total,row.t_solve,row.t_parse,row.t_build,row.status,row.ndv,row.m,datestr(now,'HH:MM:SS'),row.err); % MMP, 11/06/2026
    end
end
end

function row = run_one(fh)
% Time FH: total wall time, MOSEK time, status and SDP shape of its program
% row = struct('t_total',NaN,'t_solve',NaN,'t_parse',NaN,'nrep',0,'parse_reps','','status','', ...
%              'ndv',NaN,'m',NaN,'Ks','','sumKs2',NaN,'err','');            % MMP, 11/06/2026 (was)
row = struct('t_total',NaN,'t_solve',NaN,'t_parse',NaN,'t_build',NaN,'t_lpisolve',NaN, ...
             'nrep',0,'parse_reps','','build_reps','','status','', ...
             'ndv',NaN,'m',NaN,'Ks','','sumKs2',NaN,'err','');              % MMP, 11/06/2026
try
    tt = [];    ts = [];    tb = [];    tl = [];                            % MMP, 11/06/2026
    for r = 1:3
%       t0 = tic;   prog = fh();    t = toc(t0);    % no evalc: R2025b crashes on evalc around handles in a loop % MMP, 11/06/2026 (was)
        t0 = tic;   [prog,info] = fh();     t = toc(t0);    % no evalc: R2025b crashes on evalc around handles in a loop % MMP, 11/06/2026
        tt(end+1) = t;  ts(end+1) = prog.solinfo.info.cpusec; %#ok<AGROW>
        % Build and LPISOLVE times of the passes (one pass at max_raise = 0), % MMP, 11/06/2026
        % timed by LPI_SOLVE_LOOP_SOP: what the executive spent before and in % MMP, 11/06/2026
        % the solver call, apart from the work after it                     % MMP, 11/06/2026
        b = NaN;    l = NaN;                                                % MMP, 11/06/2026
        if isstruct(info) && isfield(info,'hist') && ~isempty(info.hist)    % MMP, 11/06/2026
            b = sum([info.hist.t_build]);   l = sum([info.hist.t_solve]);   % MMP, 11/06/2026
        end                                                                 % MMP, 11/06/2026
        tb(end+1) = b;  tl(end+1) = l; %#ok<AGROW>                          % MMP, 11/06/2026
        if t>=10,   break,  end
    end
    row.t_total = median(tt);   row.t_solve = median(ts);   row.t_parse = median(tt-ts);
    row.t_build = median(tb);   row.t_lpisolve = median(tl);                % MMP, 11/06/2026
    row.nrep = numel(tt);
    row.parse_reps = mat2str(tt-ts,4);      % every repetition, first call first
    row.build_reps = mat2str(tb,4);                                         % MMP, 11/06/2026
    I = prog.solinfo.info;
    if I.pinf,          row.status = 'infeasible';
    elseif I.dinf,      row.status = 'unbounded';
    elseif I.numerr==2, row.status = 'unknown';
    elseif I.numerr==1, row.status = 'inaccurate';
    else,               row.status = 'optimal';
    end
    try
        S = lpi_shape_sop(prog);
        row.ndv = S.ndv;    row.m = S.m;    row.Ks = mat2str(S.Ks);    row.sumKs2 = sum(S.Ks.^2);
    catch
    end
catch ME
    row.err = ME.message;
end
end

function [prog,info] = call_exec(f,nout,varargin)
% MMP, 11/06/2026: PROG, the first output of the executive F called with
% NOUT outputs, and INFO, its last output when NOUT > 1 (the container
% executives, whose info.hist holds the build and solve times), else [].
out = cell(1,nout);     [out{:}] = f(varargin{:});
prog = out{1};  info = [];  if nout>1,  info = out{end};    end
end

function st = settings_for(nv,d)
% The settings of degree d (see the header)
if nv==1
    st = lpisettings('veryheavy');
    st.sosineq_on = 0;
    st.dd1  = {d,[d d 2*d],[d d 2*d]};
    st.dd12 = st.dd1;
    st.dd2  = {d+1,[d+2 d+2 2*d+2],[d+2 d+2 2*d+2]};
    st.dd3  = {d,[d+1 d+1 2*d+1],[d+1 d+1 2*d+1]};
    st.ddZ  = [2*d 2*d 2*d];
else
    st = lpisettings('light');
    st.settings_2d.LF_deg = scale_deg(st.settings_2d.LF_deg,d);
end
end

function sn = native_from_stock(st)                                         % MMP, 11/12/2026
% SN = NATIVE_FROM_STOCK(ST) returns the stock-format 1-D settings ST in the
% native form the container executives read since 11/12/2026: lpi.pos from
% dd1 = {m,[b b j],[b b j]} and options1.sep, lpi.var.deg{1} from ddZ,
% sop.stock from dd2, options2, override2, dd3, options3; the stock-format
% fields removed. INPUT st: 'struct', SETTINGS_FOR(1,d). OUTPUT sn: 'struct'.
f = {'dd1','dd12','options1','options12','override1','ddZ','dd2','dd3', ...
     'options2','options3','override2','sosineq_on','opts'};                % MMP, 11/12/2026
sn = rmfield(st,f(isfield(st,f)));                                          % MMP, 11/12/2026
d1 = st.dd1;                                                                % MMP, 11/12/2026
sn.lpi.pos = struct('mult',[d1{1} 0 0],'int',[d1{2}(1) 1 1], ...
                    'joint',[d1{2}(3) 2 2],'sep',st.options1.sep);          % MMP, 11/12/2026
sn.lpi.var = struct('deg',{{st.ddZ,1,1}});                                  % MMP, 11/12/2026
sn.sop.stock = struct('override2',st.override2,'options2',st.options2, ...
    'dd2',{st.dd2},'options3',st.options3,'dd3',{st.dd3});                  % MMP, 11/12/2026
end                                                                         % MMP, 11/12/2026

function D = scale_deg(D,d)
% Every numeric array inside the struct/cell D multiplied by d
if isnumeric(D),        D = D*d;
elseif iscell(D),       D = cellfun(@(x) scale_deg(x,d),D,'UniformOutput',false);
elseif isstruct(D)
    f = fieldnames(D);
    for i = 1:numel(D),  for j = 1:numel(f),  D(i).(f{j}) = scale_deg(D(i).(f{j}),d);  end,  end
end
end
