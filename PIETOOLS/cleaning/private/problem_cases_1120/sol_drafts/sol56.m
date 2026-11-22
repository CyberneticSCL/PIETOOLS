% SOL56 runs PDE 56 (2-D heat with an interior delay, 3-D delay state; 0.9 tau*) in the direct stability form 'Q'
% at light with dD 0, the cap and deg_Q{3} = [1 2 2] – the degree that certifies its 1-D twin PDE 43 – one build
% and one solve on a Sol node. DRAFT (problem-case run of 11/20/2026), not submitted. Result: exs_SOL56.csv in this
% folder. MMP, 11/20/2026.
H = '/scratch/mpeet/pietools/pde55/repo/testing_codes/private/baselines/2026-11-18_studies/private';
cd(H);
info = lb_setup();
fprintf('## root %s\n## mosek %s\n## MATLAB threads %d\n',info.root,info.mosek,maxNumCompThreads);
if isempty(info.mosek),  error('sol56:mosek','mosekopt not on the path');  end
ex_search(fullfile(H,'sol56_starts.csv'),'SOL56',struct('keys',{{'PDE:56'}}, ...
    'mcap',136000,'ndvcap',7.9e6,'rowmin',1e4,'maxbuild',1,'stopmin',1e4, ...
    'statedir','_sol56'));
