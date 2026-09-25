function bl_verify()
% bl_verify -- score every cuADMM certificate against the ORIGINAL data.
%
% The GPU arm's TSV records what cuADMM says about itself. That is the solver's
% self-report, computed from the same transformed system it was handed, so it
% cannot see an error in the transform. This pass rebuilds each certificate in
% SeDuMi cone coordinates and recomputes the residual, the primal PSD margin and
% the dual PSD margin from At/b/c as they were before the dump -- which is what
% makes the GPU numbers evidence rather than testimony.
%
% Reported per (case, tolerance):
%   rel_b_v   verified row residual. Expected ~2x cuADMM's own primal
%             infeasibility, because cuADMM's stopping rule divides by
%             (1+||b||) while this divides by ||b||, and ||b||=1 after the
%             dump's normalisation. A ratio far from 2 is a red flag.
%   psd_v     primal PSD margin over the rebuilt blocks
%   dual_v    PSD margin of the dual slack C - sum y_k A_k
%   asym      symmetry defect of the rebuilt blocks; must be ~0 by construction,
%             so anything else means the svec inverse is wrong

cuadmm_path;
HERE = fileparts(mfilename('fullpath'));
DMP  = fullfile(cuadmm_outdir(),'baseline','dumps');
TSV  = fullfile(cuadmm_outdir(),'baseline','bl_verify.tsv');
COLS = {'id','tol','status','rel_b_v','psd_v','psd_relv','dual_v','asym','obj_v','normx_v','note'};
if ~exist(TSV,'file')
    fid=fopen(TSV,'w'); fprintf(fid,'%s\n',strjoin(COLS,sprintf('\t'))); fclose(fid);
end
done = readdone(TSV);

d = dir(DMP);
for i = 1:numel(d)
    if ~d(i).isdir || any(strcmp(d(i).name,{'.','..'})), continue; end
    id = d(i).name;
    mf = fullfile(DMP,[id '.mat']);
    if ~exist(mf,'file'), continue; end
    for tol = {'1e-4','1e-6'}
        key = [id '|' tol{1}];
        if isKey(done,key), continue; end
        xf = fullfile(DMP,id,['X_opt_' tol{1} '.txt']);
        if ~exist(xf,'file'), continue; end
        row = struct('id',id,'tol',tol{1},'status','','rel_b_v',NaN,'psd_v',NaN, ...
                     'psd_relv',NaN,'dual_v',NaN,'asym',NaN,'obj_v',NaN, ...
                     'normx_v',NaN,'note','');
        try
            r = cuimport(fullfile(DMP,id),mf,['_' tol{1}]);
            row.rel_b_v=r.rel_b_norm; row.psd_v=r.psd_min; row.psd_relv=r.psd_relmin;
            row.dual_v=r.dual_min;    row.asym=r.max_asym; row.obj_v=r.obj_norm;
            row.normx_v=r.normx;      row.status='ok';
        catch ME
            row.status='ERR'; row.note=regexprep(ME.message,'[\t\r\n]+',' ');
        end
        appendrow(TSV,COLS,row);
        fprintf('VF %-18s %-5s %-4s rel=%-11s psd=%-11s dual=%-11s %s\n', ...
            row.id,row.tol,row.status,num2str(row.rel_b_v),num2str(row.psd_v), ...
            num2str(row.dual_v),row.note);
    end
end
fprintf('VFDONE\n');
end

function appendrow(TSV,COLS,row)
fid=fopen(TSV,'a'); s=cell(1,numel(COLS));
for i=1:numel(COLS)
    v=row.(COLS{i});
    if ischar(v), s{i}=v; elseif isempty(v), s{i}=''; else, s{i}=num2str(v,'%.6g'); end
end
fprintf(fid,'%s\n',strjoin(s,sprintf('\t'))); fclose(fid);
end

function d = readdone(TSV)
d = containers.Map('KeyType','char','ValueType','logical');
fid=fopen(TSV,'r'); fgetl(fid);
while true
    l=fgetl(fid); if ~ischar(l), break; end
    p=strsplit(l,sprintf('\t'),'CollapseDelimiters',false);
    if numel(p)>=3 && strcmp(p{3},'ok'), d([p{1} '|' p{2}])=true; end
end
fclose(fid);
end
