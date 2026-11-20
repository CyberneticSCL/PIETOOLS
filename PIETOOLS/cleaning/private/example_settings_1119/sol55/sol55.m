% SOL55 runs PDE 55 (3-D, two coupled heat fields, lam = 1) in the stability
% form 'P' at the tier stripped: one build and one solve on a Sol node
% (maintainer approval 11/19/2026). The start is row 1 of EXU3D2
% (cleaning/private/example_settings_1119/phaseB2_3d2_starts.csv), which the
% workstation size gate refused at 7,340,048 variables, m 105,717. The caps
% sit just above that size, so a larger program is built but not solved.
% Result: exs_SOL55.csv in this folder (ex_search log: status, MOSEK
% status, eta, feasibility ratio, variables, m, seconds).
% MMP, 11/19/2026.
H = '/scratch/mpeet/pietools/pde55/repo/testing_codes/private/baselines/2026-11-18_studies/private';
cd(H);
info = lb_setup();
fprintf('## root %s\n## mosek %s\n## MATLAB threads %d\n',info.root,info.mosek,maxNumCompThreads);
if isempty(info.mosek),  error('sol55:mosek','mosekopt not on the path');  end
ex_search(fullfile(H,'sol55_starts.csv'),'SOL55',struct('keys',{{'PDE:55'}}, ...
    'mcap',106000,'ndvcap',7.4e6,'rowmin',1e4,'maxbuild',1,'stopmin',1e4, ...
    'statedir','_sol55'));
