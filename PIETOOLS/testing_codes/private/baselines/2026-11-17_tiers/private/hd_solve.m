function hd_solve()
% HD_SOLVE, 10/28/2026: solves of the sopvar 2-D H-infinity LPI in the P
% form (the form stock PIETOOLS_Hinf_gain poses on a 2-D PIE) for PDE
% 45-47 (exact 0.11101, 0.17110, 0.20668):
%   1. 'like' slack, lpisettings_sop stripped;
%   2. 'like' slack, lpisettings_sop light;
%   3. stock slack at the stock stripped settings_2d, PDE 47 only (stock
%      stripped certified 0.20707431 on an SDP that differs in 59 rows).
here = fileparts(mfilename('fullpath'));
lbdir = fullfile(fileparts(here),'library_baseline');
addpath(lbdir);     lb_setup();
assert(~contains(which('Sedumi2Mosek'),'hinf2d'),'hd_solve: the shadow Sedumi2Mosek is on the path');
S = load(fullfile(lbdir,'catalog.mat'));    C = S.C;
cat_of = @(k) C(find(strcmp({C.key},k),1));
ex = struct('id',{'45','46','47'},'exact',{0.11101,0.17110,0.20668});
PIEs = cell(1,3);
for i = 1:3,    PIEs{i} = lb_build_system(cat_of(['PDE:' ex(i).id]), [], true);    end
St = lpisettings_sop('stripped',0,'','mosek');
L  = lpisettings_sop('light',0,'','mosek');
S0 = lpisettings('stripped',0,'','mosek');
Sp = lpisettings_sop('stripped',0,'','mosek');
Sp.settings_2d = S0.settings_2d;    Sp.sos_opts = S0.sos_opts;
Sp.sop.slack = 'stock';     Sp.sop.max_raise = 0;
runs = {'P like stripped', St, 1:3;
        'P like light', L, 1:3;
        'P stock-slack stripped', Sp, 3};
for r = 1:size(runs,1)
    for i = runs{r,3}
        t0 = tic;
        try
            [~,~,gam,info] = PIETOOLS_Hinf_gain_sop(PIEs{i}, runs{r,2}, 'P');
            h = info.hist;
            fprintf('##SOL PDE %s %-24s %-10s gam %-12.8g exact %.5f %7.1f s | passes %d | dw %s | m %s | build %s | solve %s\n', ...
                ex(i).id, runs{r,1}, info.status, gam, ex(i).exact, toc(t0), numel(h), mat2str([h.dw]), ...
                mat2str([h.m]), mat2str(round([h.t_build],1)), mat2str(round([h.t_solve],1)));
        catch ME
            fprintf('##SOL PDE %s %-24s ERROR %s (%s)\n', ex(i).id, runs{r,1}, ME.message, ME.stack(1).name);
        end
    end
end
end
