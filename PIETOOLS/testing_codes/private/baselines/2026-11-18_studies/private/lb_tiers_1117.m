function lb_tiers_1117(part)
% LB_TIERS_1117(PART) the library pass of 11/17/2026 under the tiers of
% settings_sopvar on the sopvar path (LB_RUN, path 'sop'): every catalog row
% (PDE 1-72 with the 3-D cases 53-55 and the mixed-domain cases 56-72 added
% today, DDE, NDSDDF), lpisettings_sop(tier,0,'','mosek'), the _sop
% executive and form PIETOOLS_auto_execute_sop uses for the library's
% executive. Two processes with their own results files, so that the 2-D
% projection of LB_RUN reads the lighter tiers it needs:
%   'A': light, heavy, veryheavy into results_A.csv (cost order);
%   'B': extreme, stripped into results_B.csv.
% Any other PART is a tier name (or cellstr), run into results_<PART>.csv.
% custom equals light in settings_sopvar (checked by lb_catalog_add's
% caller) and is not run. Size gate on every solve: LB_NDV_CAP 2e6,
% LB_M_CAP 45000 (gate/lpisolve.m); 2-D rows above a projected 600 s per
% solve are skipped (lb_run max2d). MMP, 11/17/2026.
here = fileparts(mfilename('fullpath'));
cd(here);
setappdata(0,'LB_NDV_CAP',2e6);   setappdata(0,'LB_M_CAP',45000);   % 64 GiB workstation, other sessions running MATLAB: a 4M-variable, m 41k program peaked at 39.7 GiB
switch upper(char(part))
    case 'A',   tiers = {'light','heavy','veryheavy'};   tag = 'A';
    case 'B',   tiers = {'extreme','stripped'};          tag = 'B';
    otherwise,  tiers = cellstr(part);                   tag = strjoin(tiers,'_');
end
opts = struct('path','sop','results',fullfile(here,['results_',tag,'.csv']),'max2d',600);
fprintf('## LB_TIERS part %s: %s -> %s\n',tag,strjoin(tiers,' '),opts.results);
for t = 1:numel(tiers)
    if exist(fullfile(here,'PAUSE_sop'),'file'),    disp('PAUSE_sop present: stopped');   return;     end
    lb_run('all', tiers(t), opts);
end
fprintf('## LB_TIERS part %s DONE %s\n',tag,char(datetime('now')));
end
