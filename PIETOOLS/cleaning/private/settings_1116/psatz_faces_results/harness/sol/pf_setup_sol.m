function S = pf_setup(route)
% S = PF_SETUP(ROUTE) Sol copy of the psatz_faces path setup: the git
% archive of HEAD faca4e51 at REPO on the path (cleaning/ and orphaned_file/
% removed), the testing_codes/private folders, the MOSEK 11 toolbox, then
% the scratch shadows first (poscopvar_direct with per-direction row
% offsets, lpi_ineq_sop reading the tier's rows, heatNd_lpi with the 'rows'
% term case, heatNd_solve with the echo option). ROUTE true is refused: the
% executive route shadow is not shipped. Returns the Sol scratch root.
% SCRATCH, MMP, 11/15/2026 (Sol copy of the workstation pf_setup).
if nargin<1 || isempty(route),  route = false;  end
assert(~route,'pf_setup:route','the executive route shadow is not on Sol');
S = '/scratch/mpeet/pietools/psatz_faces/psatz_faces';
R = '/scratch/mpeet/pietools/psatz_faces/repo';
MSK = '/scratch/mpeet/pietools/solvers/mosek/11.0/toolbox/r2019b';
addpath(genpath(R));
w = warning('off','MATLAB:mpath:privateDirectoriesNotAllowedOnPath');
addpath(genpath([R '/testing_codes/private']));    warning(w);
p = strsplit(path,pathsep);
bad = p(contains(p,'/cleaning') | contains(p,'/orphaned_file'));
if ~isempty(bad),   rmpath(strjoin(bad,pathsep));   end
addpath(MSK);
addpath([S '/shadow']);     addpath([S '/heat']);   addpath(S);
for f = {'poscopvar_direct','lpi_ineq_sop','heatNd_lpi','heatNd_solve','pf_terms'}
    wf = which(f{1});
    assert(startsWith(wf,S),"shadow of "+f{1}+" not first on the path: "+wf);
end
for f = {'copvar','heatNd_pie','heatNd_sdp','heatNd_lindep','lpiprogram','get_lift_degs','exec_tools_sop','Sedumi2Mosek'}
    wf = which(f{1});
    assert(startsWith(wf,R),f{1}+" does not resolve inside the repo: "+wf);
end
assert(startsWith(which('mosekopt'),MSK),"mosekopt not from "+MSK);
fprintf('pf_setup (Sol): shadows on, repo %s, MOSEK %s, MATLAB threads %d\n',R,MSK,maxNumCompThreads);
end
