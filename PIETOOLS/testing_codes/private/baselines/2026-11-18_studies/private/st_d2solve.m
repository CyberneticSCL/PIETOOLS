function st_d2solve(cd3,frac,capinc,dD,degQ3,form)
% ST_D2SOLVE(CD3,FRAC,CAPINC,DD,DEGQ3,FORM)  Study A3/A4 of PLAN_1118_studies.md: solve the
% direct stability form 'Q' at the light tier on D2, the 2-D delayed heat
% equation with its 3-D delay state (library case 56), at tau = FRAC*tau*,
% tau* = 1.034726, with the slack cap cap_delta(3) = CD3; the other knobs as
% in the Sol jobs and the size sweeps (3-D psatz 'faces1', wcross 0, Q at
% the tier degree, no raise). The cap-1 rung (2,721,692 variables, m
% 82,551) was found to factor sparsely on the workstation (Fable session,
% 11/18 05:50: dense dim 3, about 190 s per iteration within 64 GB, 4
% iterations run), so it is solved here rather than on Sol. CAPINC is
% lpi.ineq.cap_inc ([] none; 'auto' the triangular cap, which by Lemma 1.2
% of the Newton polytope note keeps every certificate on D2, whose targets
% are purely integral). Prints '@@ A3' lines; the MOSEK log goes into the
% session log. MMP, 11/18/2026.
if nargin<1,  cd3 = 1;  end
if nargin<2,  frac = 0.5;  end
if nargin<3,  capinc = [];  end       % lpi.ineq.cap_inc: [] none, 'auto' the triangular cap
if nargin<4,  dD = 1;  end           % lpi.ineq.dD, the lift raise (light: 1)
if nargin<5,  degQ3 = 1;  end        % sop.deg_Q{3}, the degree of Q in three variables (light: 1)
if nargin<6,  form = 'Q';  end       % the stability form: 'Q' direct, 'P' storage (stronger)
R = 'C:/Users/mpeet/OneDrive/Documents/GitHub/PIETOOLS/PIETOOLS';
t0 = tic;   cd(R);  pietools_path_update;
fprintf('@@ A3 start %s: cap_delta(3) %g, cap_inc %s, tau = %g tau*\n',char(datetime('now')),cd3,mat2str(capinc),frac);
fprintf('@@ A3 lpisolve %s\n',which('lpisolve'));
PDE_t = PIETOOLS_PDE_Ex_2D_Delayed_Heat_Eq_Fridman(0,{sprintf('tau=%.12g*1.034726;',frac)});
PDE = initialize_PIETOOLS_PDE(PDE_t,true);
PIE = initialize(convert(PDE,{'pie','silent'}));
fprintf('@@ A3 pie %.0f s, dim %d\n',toc(t0),PIE.dim);
st = lpisettings_sop('light');
st.sos_opts.solver = 'mosek';   st.sop.verbose = true;
st.lpi.ineq.psatz = {'auto','faces1','faces1'};
st.sop.cap_delta = [NaN NaN cd3];   st.sop.max_raise_Q(3) = 0;   st.lpi.ineq.cap_inc = capinc;
st.lpi.ineq.dD = dD;   st.sop.deg_Q{3} = degQ3;
t1 = tic;       % sossolve prints the MOSEK log into the session log
[prog,~,info] = PIETOOLS_stability_sop(PIE,st,form);
fprintf('@@ A3 form %s\n',form);
ts = toc(t1);
sh = struct('ndv',NaN,'m',NaN);   if isfield(info,'shape'),  sh = info.shape;  end
eta = NaN;  if isfield(info,'hist') && ~isempty(info.hist) && isfield(info.hist,'eta'),  eta = info.hist(end).eta;  end
fr = NaN;   if isfield(prog,'solinfo') && isfield(prog.solinfo,'info') && isfield(prog.solinfo.info,'feasratio'),  fr = prog.solinfo.info.feasratio;  end
fprintf('@@ A3 result: status %s, certified %d, eta %.3g, optimality measure %.3g, ndv %d, m %d, executive %.0f s, total %.0f s\n', ...
    info.status,info.certified,eta,fr,sh.ndv,sh.m,ts,toc(t0));
end
