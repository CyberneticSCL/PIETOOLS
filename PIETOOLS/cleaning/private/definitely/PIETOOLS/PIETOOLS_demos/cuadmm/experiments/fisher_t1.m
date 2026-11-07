% fisher_t1 -- does the function reproduce the shipped script, and is
% use_bnd=false coherent?
cuadmm_path;
fprintf('FT variant|R|m|Kf|nblk|Ns|c_nnz|types|ok|rel_b|feasratio|numerr|t\n');
for ub = [true false]
    try
        [~,I] = cuadmm_private('fisher_prog',4.01,ub,true,'mosek');
        fprintf('FT use_bnd=%d|%.2f|%d|%d|%d|[%s]|%d|%s|%d|%.3e|%+.4f|%d|%.1f\n', ...
            ub,I.R,I.m,I.Kf,I.nblk,strtrim(num2str(I.Ns)),I.c_nnz,I.types, ...
            I.ok,I.rel_b,I.feasratio,I.numerr,I.t);
    catch ME
        fprintf('FT use_bnd=%d|FAILED|%s\n',ub,strrep(ME.message,newline,' '));
    end
end
fprintf('FTDONE\n');
