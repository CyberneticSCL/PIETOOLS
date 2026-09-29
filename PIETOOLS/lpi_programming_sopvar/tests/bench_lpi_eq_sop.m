function T = bench_lpi_eq_sop(qs,Ns,reps)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% T = BENCH_LPI_EQ_SOP(QS,NS,REPS) measures the cost the dispatcher
% 'lpi_eq_sop' adds to the routine it calls, over a sweep of the number of
% decision variables q and of spatial variables N, and the cost of
% 'lpiprogram_sop' against 'lpiprogram'.
%
% For each (q, N): a square 'cdopvar' decision operator P on L2^n[s1..sN]
% from lpivar_cdopvar at degree d(N) (2, 1, 1 for N = 1, 2, 3), with the
% component count n chosen so that numel(P.Zd) is about q, on an
% lpiprogram_sop program; then lpi_eq_sop(prog,P) and lpi_eq_cdopvar(prog,P)
% alternately, REPS times each, minimum wall time kept; the two programs
% must be equal (expr At and b). At the largest q per N, the profiler's
% PeakMem of each call (memory allocated above the call's entry, the
% function and its callees).
%
% QS default [1e4 1e5 1e6], NS default 1:3, REPS default 3.
% OUTPUT table: N, n, q, t_sop, t_direct, ratio, rows, peak_sop_MB,
% peak_direct_MB.
%
% Initial coding MMP, 09/29/2026
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if nargin<1 || isempty(qs),     qs = [1e4 1e5 1e6];     end
if nargin<2 || isempty(Ns),     Ns = 1:3;               end
if nargin<3 || isempty(reps),   reps = 3;               end
warning('off','sopvar:noncanonicalMultiplier');
warning('off','sdopvar:noncanonicalMultiplier');
dg = [2 1 1 1];
rows = [];
for N = Ns
    vars = arrayfun(@(i) sprintf('s%d',i),1:N,'uni',0);
    D = repmat([0 1],N,1);
    % Decision variables per component pair at degree d: count once at n = 1.
    [~,P1] = lpivar_cdopvar(lpiprogram_sop(vars,D),1,{vars},D,dg(N));
    c1 = numel(P1.Zd);
    for q = qs
        n = max(1,round(sqrt(q/c1)));
        [prog,P] = lpivar_cdopvar(lpiprogram_sop(vars,D),n,{vars},D,dg(N));
        tq = inf(1,2);
        for r = 1:reps
            t0 = tic;   pa = lpi_eq_sop(prog,P);        tq(1) = min(tq(1),toc(t0));
            t0 = tic;   pb = lpi_eq_cdopvar(prog,P);    tq(2) = min(tq(2),toc(t0));
        end
        same = isequal(pa.expr.At,pb.expr.At) && isequal(pa.expr.b,pb.expr.b);
        assert(same,'bench_lpi_eq_sop: programs differ at N = %d, q = %d',N,numel(P.Zd))
        nr = 0;     for i = 1:pa.expr.num,  nr = nr + numel(pa.expr.b{i});   end
        pk = [NaN NaN];
        if q==max(qs)
            pk(1) = peak_mb(@() lpi_eq_sop(prog,P),'lpi_eq_sop');
            pk(2) = peak_mb(@() lpi_eq_cdopvar(prog,P),'lpi_eq_cdopvar');
        end
        rows = [rows; N, n, numel(P.Zd), tq, tq(1)/tq(2), nr, pk];                  %#ok<AGROW>
        fprintf('N=%d n=%4d q=%8d: lpi_eq_sop %.4f s, lpi_eq_cdopvar %.4f s (x%.3f), rows %d, peak %.1f / %.1f MB\n',...
            N,n,numel(P.Zd),tq(1),tq(2),tq(1)/tq(2),nr,pk(1),pk(2));
        clear prog P pa pb
    end
end
T = array2table(rows,'VariableNames',{'N','n','q','t_sop','t_direct','ratio','rows','peak_sop_MB','peak_direct_MB'});

% lpiprogram_sop against lpiprogram (N <= 2) and alone (N up to 9).
pvar s1 s2
t = inf(1,2);
for r = 1:20
    t0 = tic;   lpiprogram_sop([s1;s2],[0 1]);  t(1) = min(t(1),toc(t0));
    t0 = tic;   lpiprogram([s1;s2],[0 1]);      t(2) = min(t(2),toc(t0));
end
fprintf('lpiprogram_sop %.2f ms, lpiprogram %.2f ms (N = 2, min of 20)\n',1e3*t(1),1e3*t(2));
for N = [1 3 5 9]
    vars = arrayfun(@(i) sprintf('s%d',i),1:N,'uni',0);     tn = inf;
    for r = 1:20,   t0 = tic;   lpiprogram_sop(vars,[0 1]);   tn = min(tn,toc(t0));   end
    fprintf('lpiprogram_sop N = %d: %.2f ms\n',N,1e3*tn);
end
end


function mb = peak_mb(f,name)
% Profiler PeakMem (MB) of the function NAME called by f.
profile clear;  profile('-memory','on');
f();
profile off;    p = profile('info');
k = find(strcmp({p.FunctionTable.FunctionName},name),1);
mb = NaN;   if ~isempty(k),     mb = p.FunctionTable(k).PeakMem/2^20;   end
end
