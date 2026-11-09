function hd_l1v()
% HD_L1V, 10/28/2026: build-only check of the one-variable-space
% translation (poslpivar_settings_2d_sop, working tree) against stock, on
% the 2-D stability LPI of PDE 42 (R^n, L2[s1], L2[s2], L2[s1,s2]) and PDE
% 43 (L2[s1], L2[s1,s2]): stock PIETOOLS_PIE2PDEstability (-> _2D) against
% PIETOOLS_stability_sop 'Q' (posed as 'P' in 2-D) with the stock slack, on
% the same settings_2d:
%   'stock stripped'   lpisettings stripped unchanged (all caps x = theta);
%   'asym caps'        the same with LF_deg and eq_deg dx{2}, dx{3}, dy{2},
%                      dy{3} given unequal caps in the variable and in theta,
%                      so a swap of the two in the translation shows.
here = fileparts(mfilename('fullpath'));
lbdir = fullfile(fileparts(here),'library_baseline');
addpath(lbdir);     lb_setup();
addpath(fullfile(here,'shadow'),'-begin');
w = which('Sedumi2Mosek','-all');
assert(startsWith(w{1},fullfile(here,'shadow')),'hd_l1v: shadow Sedumi2Mosek is not first on the path');
fprintf('translator: %s\n', which('poslpivar_settings_2d_sop'));
global HD_FP
S = load(fullfile(lbdir,'catalog.mat'));    C = S.C;
cat_of = @(k) C(find(strcmp({C.key},k),1));

S0 = lpisettings('stripped',0,'','mosek');
A0 = S0;                                % asymmetric caps: [var; theta; joint]
A0.settings_2d.LF_deg.dx{2} = [2;1;2];  A0.settings_2d.LF_deg.dx{3} = [1;2;2];
A0.settings_2d.LF_deg.dy{2} = [2,1,2];  A0.settings_2d.LF_deg.dy{3} = [1,2,2];
A0.settings_2d.eq_deg.dx{2} = [4;3;4];  A0.settings_2d.eq_deg.dx{3} = [3;4;4];
A0.settings_2d.eq_deg.dy{2} = [4,3,4];  A0.settings_2d.eq_deg.dy{3} = [3,4,4];
sets = {'stock stripped', S0; 'asym caps', A0};
for id = {'PDE:42','PDE:43'}
    PIE = lb_build_system(cat_of(id{1}), [], true);
    fprintf('%s T.dim %s\n', id{1}, mat2str(PIE.T.dim));
    for k = 1:size(sets,1)
        Ss = sets{k,2};
        Sp = lpisettings_sop('stripped',0,'','mosek');
        Sp.settings_2d = Ss.settings_2d;    Sp.sos_opts = Ss.sos_opts;
        Sp.sop.slack = 'stock';     Sp.sop.max_raise = 0;
        cases = {'stock', @() PIETOOLS_PDEstability(PIE,Ss);       % P form (_2D); PIE2PDEstability_2D is the Q form
                 'sop',   @() PIETOOLS_stability_sop(PIE,Sp,'Q')};
        for c = 1:2
            HD_FP = [];     t0 = tic;
            lbl = sprintf('%s %s %s', id{1}, sets{k,1}, cases{c,1});
            try
                cases{c,2}();
                fprintf('##FP %-30s returned without reaching the solver\n', lbl);
            catch ME
                t = toc(t0);
                if ~strcmp(ME.identifier,'hdshadow:stop')
                    fprintf('##FP %-30s ERROR %s (%s)\n', lbl, ME.message, ME.stack(1).name);
                    continue
                end
                f = HD_FP;
                fprintf(['##FP %-30s build %5.1f s | m %6d | scal %8d | Kf %d | Ks [%s] | nnzA %d | ' ...
                         'sumabsA %.10g | sumsqA %.10g | rownnz2 %d | normb %.10g\n'], ...
                        lbl, t, f.m, f.scal, f.Kf, num2str(f.Ks), f.nnzA, f.sumabsA, f.sumsqA, f.rownnz2, f.normb);
            end
        end
    end
end
end
