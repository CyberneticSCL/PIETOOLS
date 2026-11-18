function S = pf_setup(route)
% S = PF_SETUP(ROUTE) puts the scratch shadows of the psatz_faces run first
% on the path (poscopvar_direct with per-direction row offsets, lpi_ineq_sop
% reading the tier's rows, heatNd_lpi with the 'rows' term case), adds
% testing_codes/private and MOSEK, removes cleaning/ and orphaned_file/,
% and asserts the resolution of every shadowed name. ROUTE true also puts
% the route shadow of exec_tools_sop (the N-D default storage in 2-D) first.
% Returns the scratch root. SCRATCH, MMP, 11/15/2026.
if nargin<1 || isempty(route),  route = false;  end
S = 'C:/Users/mpeet/AppData/Local/Temp/claude/363c499c6332a503/af90e0a1-bcff-4e8d-9a24-9aa4c52491cc/scratchpad/psatz_faces';
R = 'C:/Users/mpeet/OneDrive/Documents/GitHub/PIETOOLS/PIETOOLS';
w = warning('off','MATLAB:mpath:privateDirectoriesNotAllowedOnPath');
addpath(genpath([R '/testing_codes/private']));    warning(w);
p = strsplit(path,pathsep);
bad = p(contains(p,[filesep 'cleaning']) | contains(p,[filesep 'orphaned_file']) | contains(p,'psatz_mix') | contains(p,'psatz_faces'));
if ~isempty(bad),   rmpath(strjoin(bad,pathsep));   end
d = dir('C:/Program Files/Mosek/*/toolbox/r20*');   addpath(fullfile(d(end).folder,d(end).name));
addpath([S '/shadow']);     addpath([S '/heat']);   addpath(S);
if route,   addpath([S '/route']);  end
names = {'poscopvar_direct','lpi_ineq_sop','heatNd_lpi','pf_terms'};
if route,   names{end+1} = 'exec_tools_sop';    end
for f = names
    wf = strrep(which(f{1}),'\','/');
    assert(startsWith(wf,S),"shadow of "+f{1}+" not first on the path: "+wf);
end
if ~route
    wf = strrep(which('exec_tools_sop'),'\','/');
    assert(startsWith(wf,R),"exec_tools_sop should be the repo's: "+wf);
end
for f = {'copvar','heatNd_pie','heatNd_solve','heatNd_bisect','lpiprogram','PIETOOLS_stability_sop','mixed_domain_systems','get_lift_degs'}
    wf = strrep(which(f{1}),'\','/');
    assert(startsWith(wf,R),f{1}+" does not resolve inside the repo: "+wf);
end
fprintf('pf_setup: shadows on (route %d), repo %s, MOSEK %s\n',route,R,d(end).name);
end
