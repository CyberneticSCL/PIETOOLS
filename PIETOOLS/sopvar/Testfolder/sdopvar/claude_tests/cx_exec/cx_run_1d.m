function R = cx_run_1d(families)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% R = CX_RUN_1D(FAMILIES) runs every 1-D executive comparison serially:
% stock opvar executive against its container transcription cx_<exec>, on
% the baseline plants, via cx_compare. FAMILIES is a cellstr subset of
% {'stability','hinf','h2'} (default: all three).
%
% Serial on purpose: overlapping runs corrupt the timings (the baseline
% suite's concurrency lesson), and timings are part of the comparison.
% Measured 2026-09-25: 154 s for all 21 cases (Mosek, i9-14900KF).
%
% Initial coding MMP, 09/25/2026
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if nargin<1 || isempty(families),   families = {'stability','hinf','h2'};   end
C = {};
if any(strcmp(families,'stability')),   C = [C, reshape(num2cell(cx_cases_stability('1d')),1,[])];  end
if any(strcmp(families,'hinf')),        C = [C, reshape(num2cell(cx_cases_hinf()),1,[])];            end
if any(strcmp(families,'h2')),          C = [C, reshape(num2cell(cx_cases_h2()),1,[])];             end
% The families' case structs carry the same fields in different orders.
fl = {'id','exec','setname','solver','kind','plant','thresh'};
for k = 1:numel(C)
    s = struct();
    for f = fl,     if isfield(C{k},f{1}), s.(f{1}) = C{k}.(f{1}); else, s.(f{1}) = []; end,  end
    C{k} = s;
end
t0 = tic;
fprintf('cx_run_1d: %d cases, shapes as stock/container\n',numel(C));
R = cx_compare([C{:}]);
fprintf('cx_run_1d: %.1f s total\n',toc(t0));
end
