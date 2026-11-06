function export_scaling_csv(f23,f4,outdir,fof)
% EXPORT_SCALING_CSV(F23,F4,OUTDIR[,FOF]) writes the data of the scaling
% plots as CSV files in OUTDIR:
%   ops_vs_q.csv                   one row per (sweep, q, operation, class):
%                                  sweep ('nnz_fixed', 'nnz_fixed_shared',
%                                  'nnz_4q'), n, d, q, nnz, op, class, t
%                                  (median s), nrep, bytes;
%   stability_parse_vs_degree.csv  one row per (nv, d, series, parse_only)
%                                  of the stability LPI;
%   of_parse_vs_degree.csv         the same for the output feedback LPI (only
%                                  when FOF is given).
% The two LPI files have the columns problem, nv, d, series, parse_only
% (1: the solver was the stub MOSEKOPT, t_solve = 0), t_total, t_solve,
% t_parse (s), pct_parse (100 t_parse/t_total; NaN when parse_only),
% t_build and t_lpisolve (s, the LPI construction and the LPISOLVE call of
% a container executive; NaN for legacy and for rows measured before these
% fields), nrep, parse_reps (parsing time of every call, first call first),
% build_reps, status, ndv, m, Ks, sumKs2, err.
% INPUT
% - f23:    char, MAT-file with Rq2 (dvarQ), Rs2 (dvarS), Rn2 (dvarN) of
%           BENCH_CLASS_OPS;
% - f4:     char, MAT-file whose struct arrays (any names) hold the rows of
%           BENCH_LPI_DEGREE for the stability LPI;
% - outdir: char, output folder;
% - fof:    (optional) char, MAT-file whose struct arrays hold the rows of
%           BENCH_LPI_DEGREE for the output feedback LPI.
% OUTPUT
% none; the files.

% Initial coding MMP, 11/06/2026: data of the scaling plots (README.md).
% MMP, 11/06/2026: the LPI tables by one local function LPI_CSV, which
%                replaces the inline stability block; the rows carry problem,
%                parse_only, t_build and t_lpisolve (rows saved before those
%                fields existed are stability rows with a full solve, without
%                build times); the output feedback LPI to
%                of_parse_vs_degree.csv.

S = load(f23);
map = {'Rq2','nnz_fixed'; 'Rs2','nnz_fixed_shared'; 'Rn2','nnz_4q'};
C = {};
for k = 1:size(map,1)
    if ~isfield(S,map{k,1}),    continue,   end
    R = S.(map{k,1});
    for i = 1:numel(R)
        r = R(i);
        C(end+1,:) = {map{k,2},r.n,r.d,r.q,r.nnz,r.op,r.cls,r.t,r.nrep,r.bytes}; %#ok<AGROW>
    end
end
T = cell2table(C,'VariableNames',{'sweep','n','d','q','nnz','op','class','t','nrep','bytes'});
writetable(T,fullfile(outdir,'ops_vs_q.csv'));

% MMP, 11/06/2026: start of the LPI tables. Replaces the inline block of
% 7688dcf2, which wrote R1 and R2 of F4 by STRUCT2TABLE (stability, full
% solves only).
lpi_csv(f4,'stability',fullfile(outdir,'stability_parse_vs_degree.csv'));
if nargin>=4 && ~isempty(fof)
    lpi_csv(fof,'of',fullfile(outdir,'of_parse_vs_degree.csv'));
end
fprintf('##wrote the CSV files in %s\n',outdir);
end

function lpi_csv(f,problem,out)
% Every struct array of the MAT-file F as one table, sorted by nv, d, series
cols = {'problem','nv','d','series','parse_only','t_total','t_solve','t_parse','pct_parse', ...
        't_build','t_lpisolve','nrep','parse_reps','build_reps','status','ndv','m','Ks','sumKs2','err'};
S = load(f);    C = {};
for fn = fieldnames(S)'
    R = S.(fn{1});
    if ~isstruct(R),    continue,   end
    for i = 1:numel(R)
        r = R(i);
        if ~isfield(r,'problem'),       r.problem = problem;    end
        if ~isfield(r,'parse_only'),    r.parse_only = false;   end
        if ~isfield(r,'parse_reps'),    r.parse_reps = '';      end
        if ~isfield(r,'t_build'),       r.t_build = NaN;        end
        if ~isfield(r,'t_lpisolve'),    r.t_lpisolve = NaN;     end
        if ~isfield(r,'build_reps'),    r.build_reps = '';      end
        pct = NaN;  if ~r.parse_only,   pct = 100*r.t_parse/r.t_total;  end
        C(end+1,:) = {r.problem,r.nv,r.d,r.series,r.parse_only,r.t_total,r.t_solve,r.t_parse,pct, ...
                      r.t_build,r.t_lpisolve,r.nrep,r.parse_reps,r.build_reps,r.status,r.ndv,r.m,r.Ks,r.sumKs2,r.err}; %#ok<AGROW>
    end
end
T = cell2table(C,'VariableNames',cols);
T = sortrows(T,{'nv','d','series','parse_only'});
writetable(T,out);
end
% MMP, 11/06/2026: end of the LPI tables.
