function R = bench_lpi_degree(nv,dlist,series)
% R = BENCH_LPI_DEGREE(NV,DLIST,SERIES) times the parsing and the solving of
% the stability LPI in the direct form, P > 0 and T'PA + A'PT < 0, for the
% heat equation
%     x_t = sum_{i<=NV} x_{s_i s_i} + lam x   on [0,1]^NV, x = 0 on the boundary
% (lam = pi^2/2 for NV = 1 and pi^2 for NV = 2, both stable), as a function
% of the degree d of the Lyapunov variable P, with the legacy classes and
% with the container classes. Parsing time = wall time of the executive
% minus the SDP solver time reported by MOSEK.
% INPUT
% - nv:     1 or 2, the number of spatial variables;
% - dlist:  vector of positive integers, the Lyapunov degrees d;
% - series: cellstr, any of
%             'legacy'    LEGACY_PDESTABILITY / LEGACY_PDESTABILITY_2D (the
%                         stock executive) on PIE2OPVAR(PIE): opvar, dopvar
%                         (opvar2d, dopvar2d) and lpi_programming;
%             'sop_stock' PIETOOLS_STABILITY_SOP(PIE,settings,'P') with
%                         settings.sop.slack = 'stock', the stock SDP (1-D);
%             'sop'       PIETOOLS_STABILITY_SOP(PIE,settings,'P') with its
%                         default slack.
% OUTPUT
% - R: struct array, one row per (d, series), fields nv, d, series,
%      t_total, t_solve, t_parse (median seconds over the repetitions),
%      nrep, status, ndv (decision variables), m (equality rows), Ks (PSD
%      block sizes, as text), sumKs2 (sum of squared PSD sizes), err.
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

% The PDE and its PIE (copvar operators), and the PIE with opvar operators
pvar s1 s2
if nv==1
    x = pde_var(1,s1,[0,1]);
    PDE = [diff(x,'t')==diff(x,s1,2)+pi^2/2*x; subs(x,s1,0)==0; subs(x,s1,1)==0];
else
    x = pde_var(1,[s1;s2],[0,1;0,1]);
    PDE = [diff(x,'t')==diff(x,s1,2)+diff(x,s2,2)+pi^2*x;
           subs(x,s1,0)==0; subs(x,s1,1)==0; subs(x,s2,0)==0; subs(x,s2,1)==0];
end
PIE = convert(PDE,'pie');
PIEo = pie2opvar(PIE);

R = [];
for d = dlist
    st = settings_for(nv,d);
    for k = 1:numel(series)
        ser = series{k};
        stk = st;
        stk.sop = struct('max_raise',0,'max_raise_feas',0,'witness',false,'keepdual',false,'verbose',false);
        if strcmp(ser,'sop_stock'),     stk.sop.slack = 'stock';    end
        switch ser
            case 'legacy'
                if nv==1,   fh = @() legacy_PDEstability(PIEo,stk);
                else,       fh = @() legacy_PDEstability_2D(PIEo,stk);
                end
            otherwise
                fh = @() PIETOOLS_stability_sop(PIE,stk,'P');
        end
        row = run_one(fh);
        row.nv = nv;    row.d = d;      row.series = ser;
        R = [R; orderfields(row)]; %#ok<AGROW>
        fprintf('##nv %d d %d %-9s total %8.2f s  solve %8.2f s  parse %8.2f s  %-10s ndv %8d  m %7d  %s %s\n', ...
            nv,d,ser,row.t_total,row.t_solve,row.t_parse,row.status,row.ndv,row.m,datestr(now,'HH:MM:SS'),row.err);
    end
end
end

function row = run_one(fh)
% Time FH: total wall time, MOSEK time, status and SDP shape of its program
row = struct('t_total',NaN,'t_solve',NaN,'t_parse',NaN,'nrep',0,'parse_reps','','status','', ...
             'ndv',NaN,'m',NaN,'Ks','','sumKs2',NaN,'err','');
try
    tt = [];    ts = [];
    for r = 1:3
        t0 = tic;   prog = fh();    t = toc(t0);    % no evalc: R2025b crashes on evalc around handles in a loop
        tt(end+1) = t;  ts(end+1) = prog.solinfo.info.cpusec; %#ok<AGROW>
        if t>=10,   break,  end
    end
    row.t_total = median(tt);   row.t_solve = median(ts);   row.t_parse = median(tt-ts);
    row.nrep = numel(tt);
    row.parse_reps = mat2str(tt-ts,4);      % every repetition, first call first
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

function D = scale_deg(D,d)
% Every numeric array inside the struct/cell D multiplied by d
if isnumeric(D),        D = D*d;
elseif iscell(D),       D = cellfun(@(x) scale_deg(x,d),D,'UniformOutput',false);
elseif isstruct(D)
    f = fieldnames(D);
    for i = 1:numel(D),  for j = 1:numel(f),  D(i).(f{j}) = scale_deg(D(i).(f{j}),d);  end,  end
end
end
