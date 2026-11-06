function export_scaling_csv(f23,f4,outdir)
% EXPORT_SCALING_CSV(F23,F4,OUTDIR) writes the data of the scaling plots as
% two CSV files in OUTDIR:
%   ops_vs_q.csv                 one row per (sweep, q, operation, class):
%                                sweep ('nnz_fixed', 'nnz_fixed_shared',
%                                'nnz_4q'), n, d, q, nnz, op, class, t
%                                (median s), nrep, bytes;
%   stability_parse_vs_degree.csv  one row per (nv, d, series): nv, d,
%                                series, t_total, t_solve, t_parse (s),
%                                pct_parse (100 t_parse/t_total), nrep,
%                                parse_reps (parsing time of every call,
%                                first call first), status, ndv, m, Ks,
%                                sumKs2, err.
% INPUT
% - f23:    char, MAT-file with Rq2 (dvarQ), Rs2 (dvarS), Rn2 (dvarN) of
%           BENCH_CLASS_OPS;
% - f4:     char, MAT-file with R1 (1-D) and R2 (2-D) of BENCH_LPI_DEGREE;
% - outdir: char, output folder.
% OUTPUT
% none; the two files.

% Initial coding MMP, 11/06/2026: data of the scaling plots (README.md).

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

S = load(f4);
R = [];
for fn = {'R1','R2'}
    if isfield(S,fn{1}),    R = [R; S.(fn{1})(:)]; end %#ok<AGROW>
end
T = struct2table(R);
T.pct_parse = 100*T.t_parse./T.t_total;
T = T(:,{'nv','d','series','t_total','t_solve','t_parse','pct_parse','nrep','parse_reps','status','ndv','m','Ks','sumKs2','err'});
writetable(T,fullfile(outdir,'stability_parse_vs_degree.csv'));
fprintf('##wrote %s and %s\n',fullfile(outdir,'ops_vs_q.csv'),fullfile(outdir,'stability_parse_vs_degree.csv'));
end
