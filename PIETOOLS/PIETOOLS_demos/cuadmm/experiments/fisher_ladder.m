% fisher_ladder -- where does feasibility break in R?  Paper and the script's
% own comment both say ~4.048.  Three modes:
%   'opt'    minimise gam at fixed R           (c != 0, the shipped default)
%   'gamfix' bound kept with gam NUMERIC       (c  = 0, faithful to Thm 1)
%   'nobnd'  bound dropped                     (c  = 0, strictly weaker)
cuadmm_path;
RL = [1 2 3 4 4.05 4.2 5];
fprintf('FL mode|R|m|Kf|c_nnz|ok|rel_b|feasratio|numerr|t\n');
MODES = {'opt',true; 'gamfix',2.0; 'nobnd',false};
for k = 1:size(MODES,1)
    lastok = NaN;
    for R = RL
        try
            [~,I] = cuadmm_private('fisher_prog',R,MODES{k,2},true,'mosek');
            fprintf('FL %s|%.2f|%d|%d|%d|%d|%.3e|%+.4f|%d|%.1f\n', ...
                MODES{k,1},R,I.m,I.Kf,I.c_nnz,I.ok,I.rel_b,I.feasratio,I.numerr,I.t);
            if I.ok, lastok = R; end
        catch ME
            fprintf('FL %s|%.2f|ERR|%s\n',MODES{k,1},R,strrep(ME.message,newline,' '));
        end
    end
    fprintf('FLMAX %s|last_feasible_R=%.2f\n',MODES{k,1},lastok);
end
fprintf('FLDONE\n');
