function T = bench_lpi_eq_sop(qs,Ns,reps)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% T = BENCH_LPI_EQ_SOP(QS,NS,REPS) measures the wall time and peak memory
% of lpi_eq_sop(prog,P) and of lpi_eq_cdopvar(prog,P), the routine it
% dispatches to, for a square 'cdopvar' decision operator P, over a sweep
% of the number q of decision variables and N of spatial variables, and
% the wall time of lpiprogram_sop and of lpiprogram. For each N in NS and
% q in QS, P maps L2^n[s1,...,sN] to itself, [0,1] in each variable, and
% is declared by lpivar_cdopvar at degree 2 (N = 1) or 1 (N >= 2) on the
% program lpiprogram_sop({'s1',...,'sN'},D), with
% n = max(1,round(sqrt(q/c1))), c1 the number of decision variables at
% n = 1. The two routines run alternately, REPS times each, the least wall
% time of each is kept, and the two programs must have equal expr.At and
% expr.b. At q = max(QS), the MATLAB profiler records the PeakMem of each
% call (the memory allocated above the entry of the call, by the function
% and its callees). Last, the least of 20 wall times of lpiprogram_sop and
% of lpiprogram on [s1;s2] with domain [0,1], and of lpiprogram_sop alone
% on 1, 3, 5 and 9 variables, are printed.
%
% INPUT
% - qs:     (optional) 'double' vector, the target numbers of decision
%           variables; default [1e4 1e5 1e6];
% - Ns:     (optional) 'double' vector, the numbers of spatial variables,
%           each in 1..4; default 1:3;
% - reps:   (optional) 'double', the number of timed calls of each
%           routine; default 3.
%
% OUTPUT
% - T:      'table' with one row per (N,q) and columns N, n, q (the number
%           of decision variables of P, numel(P.Zd)), t_sop and t_direct
%           (seconds, lpi_eq_sop and lpi_eq_cdopvar), ratio
%           (t_sop/t_direct), rows (the equality rows of the program),
%           peak_sop_MB and peak_direct_MB (PeakMem in MB, NaN for q below
%           max(QS)); one printed line per row and per lpiprogram timing.
%
% NOTES
% Two programs that differ are an error.
%
% For support, contact M. Peet, Arizona State University at mpeet@asu.edu

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PIETOOLS - bench_lpi_eq_sop
%
% Copyright (C) 2026 PIETOOLS Team
%
% This program is free software; you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation; either version 2 of the License, or
% (at your option) any later version.
%
% This program is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
%
% You should have received a copy of the GNU General Public License
% along with this program; if not, write to the Free Software
% Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% If you modify this code, document all changes carefully and include date
% authorship, and a brief description of modifications
%
% Initial coding MMP, 09/29/2026
% MMP, 10/28/2026: Header to the header protocol (what the function does,
%                inputs, outputs), the licence block added, step comments
%                in the body and the header of the local function; code
%                unchanged. No header text was deleted.

% Defaults; no warnings for non-canonical multipliers
if nargin<1 || isempty(qs),     qs = [1e4 1e5 1e6];     end
if nargin<2 || isempty(Ns),     Ns = 1:3;               end
if nargin<3 || isempty(reps),   reps = 3;               end
warning('off','sopvar:noncanonicalMultiplier');
warning('off','sdopvar:noncanonicalMultiplier');
% Degree of P for N = 1, 2, 3, 4
dg = [2 1 1 1];
rows = [];
% Sweep over (N, q)
for N = Ns
    vars = arrayfun(@(i) sprintf('s%d',i),1:N,'uni',0);
    D = repmat([0 1],N,1);
    % Decision variables per component pair at degree d: count once at n = 1.
    [~,P1] = lpivar_cdopvar(lpiprogram_sop(vars,D),1,{vars},D,dg(N));
    c1 = numel(P1.Zd);
    for q = qs
        % Component count for about q decision variables; the program and P
        n = max(1,round(sqrt(q/c1)));
        [prog,P] = lpivar_cdopvar(lpiprogram_sop(vars,D),n,{vars},D,dg(N));
        % Least wall time of each routine over REPS alternating calls
        tq = inf(1,2);
        for r = 1:reps
            t0 = tic;   pa = lpi_eq_sop(prog,P);        tq(1) = min(tq(1),toc(t0));
            t0 = tic;   pb = lpi_eq_cdopvar(prog,P);    tq(2) = min(tq(2),toc(t0));
        end
        % The two programs must be equal; count the equality rows
        same = isequal(pa.expr.At,pb.expr.At) && isequal(pa.expr.b,pb.expr.b);
        assert(same,'bench_lpi_eq_sop: programs differ at N = %d, q = %d',N,numel(P.Zd))
        nr = 0;     for i = 1:pa.expr.num,  nr = nr + numel(pa.expr.b{i});   end
        % Peak memory of each call, at the largest q only
        pk = [NaN NaN];
        if q==max(qs)
            pk(1) = peak_mb(@() lpi_eq_sop(prog,P),'lpi_eq_sop');
            pk(2) = peak_mb(@() lpi_eq_cdopvar(prog,P),'lpi_eq_cdopvar');
        end
        % Record and print the row
        rows = [rows; N, n, numel(P.Zd), tq, tq(1)/tq(2), nr, pk];                  %#ok<AGROW>
        fprintf('N=%d n=%4d q=%8d: lpi_eq_sop %.4f s, lpi_eq_cdopvar %.4f s (x%.3f), rows %d, peak %.1f / %.1f MB\n',...
            N,n,numel(P.Zd),tq(1),tq(2),tq(1)/tq(2),nr,pk(1),pk(2));
        clear prog P pa pb
    end
end
% The table of the sweep
T = array2table(rows,'VariableNames',{'N','n','q','t_sop','t_direct','ratio','rows','peak_sop_MB','peak_direct_MB'});

% lpiprogram_sop against lpiprogram (N <= 2) and alone (N up to 9).
pvar s1 s2
% N = 2: least of 20 wall times of each
t = inf(1,2);
for r = 1:20
    t0 = tic;   lpiprogram_sop([s1;s2],[0 1]);  t(1) = min(t(1),toc(t0));
    t0 = tic;   lpiprogram([s1;s2],[0 1]);      t(2) = min(t(2),toc(t0));
end
fprintf('lpiprogram_sop %.2f ms, lpiprogram %.2f ms (N = 2, min of 20)\n',1e3*t(1),1e3*t(2));
% lpiprogram_sop alone, N = 1, 3, 5, 9
for N = [1 3 5 9]
    vars = arrayfun(@(i) sprintf('s%d',i),1:N,'uni',0);     tn = inf;
    for r = 1:20,   t0 = tic;   lpiprogram_sop(vars,[0 1]);   tn = min(tn,toc(t0));   end
    fprintf('lpiprogram_sop N = %d: %.2f ms\n',N,1e3*tn);
end
end


function mb = peak_mb(f,name)
% MB = PEAK_MB(F,NAME) returns the PeakMem, in MB, that the MATLAB profiler
% records for the function NAME during the call F().
% INPUT  f: 'function_handle', called with no argument; name: 'char', the
%        name of the profiled function.
% OUTPUT mb: 'double', NaN when NAME is not in the profile.
% Profile the call with memory statistics
profile clear;  profile('-memory','on');
f();
profile off;    p = profile('info');
% PeakMem of NAME
k = find(strcmp({p.FunctionTable.FunctionName},name),1);
mb = NaN;   if ~isempty(k),     mb = p.FunctionTable(k).PeakMem/2^20;   end
end
