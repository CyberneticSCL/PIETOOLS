function hd_build()
% HD_BUILD, 10/28/2026: build-only comparison of the 2-D H-infinity LPI of
% PDE 45 (reaction-diffusion, exact 0.11101) on the stock and sopvar paths
% at matched degrees. The shadow Sedumi2Mosek records the SDP that sossolve
% hands to MOSEK and stops before the solve.
%   1. stock P form   PIETOOLS_Hinf_gain (-> _2D), lpisettings stripped
%   2. sop P form     stock slack, stock stripped settings_2d and sos_opts
%   3. stock Q form   PIETOOLS_Hinf_gain_2D_non_coercive, stripped
%   4. sop Q form     stock slack, as 2
%   5. sop P form     'like' slack, lpisettings_sop light
%   6. sop P form     'like' slack, lpisettings_sop stripped
here = fileparts(mfilename('fullpath'));
lbdir = fullfile(fileparts(here),'library_baseline');
addpath(lbdir);     lb_setup();
addpath(fullfile(here,'shadow'),'-begin');
w = which('Sedumi2Mosek','-all');
assert(startsWith(w{1},fullfile(here,'shadow')),'hd_build: shadow Sedumi2Mosek is not first on the path');
global HD_FP
S = load(fullfile(lbdir,'catalog.mat'));    C = S.C;
cat_of = @(k) C(find(strcmp({C.key},k),1));
PIE = lb_build_system(cat_of('PDE:45'), [], true);

% Nonzero parameter cells of T and A on L2[s1,s2] (3x3 grid: 1 multiplier,
% 2 lower, 3 upper integral; rows s1, columns s2)
nzc = @(x) ~(isempty(x) || (isnumeric(x) && all(x(:)==0)) || ...
             (isa(x,'polynomial') && all(x.coefficient(:)==0)));
fprintf('T.dim %s\n', mat2str(PIE.T.dim));
fprintf('T.R22 nonzero cells %s\n', mat2str(cellfun(nzc,PIE.T.R22)));
fprintf('A.R22 nonzero cells %s\n', mat2str(cellfun(nzc,PIE.A.R22)));

S0 = lpisettings('stripped',0,'','mosek');
Sp = lpisettings_sop('stripped',0,'','mosek');
Sp.settings_2d = S0.settings_2d;    Sp.sos_opts = S0.sos_opts;
Sp.sop.slack = 'stock';     Sp.sop.max_raise = 0;
L  = lpisettings_sop('light',0,'','mosek');
St = lpisettings_sop('stripped',0,'','mosek');
cases = {'stock P stripped',            @() PIETOOLS_Hinf_gain(PIE,S0);
         'sop P stock-slack stripped',  @() PIETOOLS_Hinf_gain_sop(PIE,Sp,'P');
         'stock Q stripped',            @() PIETOOLS_Hinf_gain_2D_non_coercive(PIE,S0);
         'sop Q stock-slack stripped',  @() PIETOOLS_Hinf_gain_sop(PIE,Sp,'Q');
         'sop P like light',            @() PIETOOLS_Hinf_gain_sop(PIE,L,'P');
         'sop P like stripped',         @() PIETOOLS_Hinf_gain_sop(PIE,St,'P')};
for k = 1:size(cases,1)
    HD_FP = [];     t0 = tic;
    try
        cases{k,2}();
        fprintf('##FP %-28s returned without reaching the solver\n', cases{k,1});
    catch ME
        t = toc(t0);
        if ~strcmp(ME.identifier,'hdshadow:stop')
            fprintf('##FP %-28s ERROR %s (%s)\n', cases{k,1}, ME.message, ME.stack(1).name);
            continue
        end
        f = HD_FP;
        fprintf(['##FP %-28s build %6.1f s | m %6d | n %8d | scal %8d | Kf %d | Ks [%s] | nnzA %d | ' ...
                 'sumabsA %.10g | sumsqA %.10g | rownnz2 %d | normb %.10g nnzb %d | nnzc %d sumc %.6g\n'], ...
                cases{k,1}, t, f.m, f.n, f.scal, f.Kf, num2str(f.Ks), f.nnzA, ...
                f.sumabsA, f.sumsqA, f.rownnz2, f.normb, f.nnzb, f.nnzc, f.sumc);
    end
end
end
