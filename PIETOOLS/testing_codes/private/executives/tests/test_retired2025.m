function R = test_retired2025()
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% R = TEST_RETIRED2025 checks the PIETOOLS 2025 names and inputs that stop
% with an error since 11/28/2026, and records what a numeric opts.psatz
% of lpi_ineq selects; errors at the first failed check:
%   (1) each of the 21 executives of PIETOOLS 2025 in
%       executives/retired_2025 resolves there once, and a call raises
%       <name>:retired with a message holding the call of this release
%       (the executive and the form of the table below);
%   (2) getController(P,Z,TOL) and getObserver(P,Z,TOL), TOL a number,
%       raise getController:call2025 and getObserver:call2025; with
%       OPTS = [] or a 'struct' the guard raises nothing;
%   (3) lpi_ineq(prog,P,opts) with opts.pure or opts.sep raises
%       lpi_ineq:call2025; with opts.psatz = 'none' it runs;
%   (4) on the 1-D operator T'*T of a heat equation on [0,1], the term
%       codes and Gram dimensions of the certificate of lpi_ineq for
%       opts.psatz = 0, 1, [0 1], 'none' and 'auto' (recorded, R.psatz):
%       a number is read as term codes (0 plain, 1 box product), so the
%       2025 call opts.psatz = 1 (the plain and the psatz term) poses the
%       box product alone.
%
% OUTPUT
% - R:  'struct' with fields n (checks passed) and psatz (5 x 1
%       'struct': the value, info.terms.codes, info.N, info.ndec).
%
% PIETOOLS - test_retired2025
%
% If you modify this code, document all changes carefully and include date
% authorship, and a brief description of modifications
%
% Initial coding MMP, 11/28/2026
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
n = 0;
% (1) the 2025 executive and the call that poses its LPI now
T = {'PIETOOLS_PDEstability',              'PIETOOLS_stability(PIE,settings,''P'')';
     'PIETOOLS_PDEstability_dual',         'PIETOOLS_stability(PIE,settings,''Pd'')';
     'PIETOOLS_PIE2PDEstability',          'PIETOOLS_stability(PIE,settings,''Q'')';
     'PIETOOLS_PIE2PDEstability_dual',     'PIETOOLS_stability(PIE,settings,''Qd'')';
     'PIETOOLS_Hinf_gain_coercive',        'PIETOOLS_Hinf_gain(PIE,settings,''P'')';
     'PIETOOLS_Hinf_gain_dual',            'PIETOOLS_Hinf_gain(PIE,settings,''Qd'')';
     'PIETOOLS_Hinf_gain_dual_coercive',   'PIETOOLS_Hinf_gain(PIE,settings,''Pd'')';
     'PIETOOLS_H2_norm_o',                 'PIETOOLS_H2_norm(PIE,settings,''o'')';
     'PIETOOLS_H2_norm_c',                 'PIETOOLS_H2_norm(PIE,settings,''c'')';
     'PIETOOLS_H2_norm_o_coercive',        'PIETOOLS_H2_norm(PIE,settings,''oco'')';
     'PIETOOLS_H2_norm_c_coercive',        'PIETOOLS_H2_norm(PIE,settings,''cco'')';
     'PIETOOLS_stability_2D',              'PIETOOLS_stability(PIE,settings,''P'')';
     'PIETOOLS_stability_dual_2D',         'PIETOOLS_stability(PIE,settings,''Pd'')';
     'PIETOOLS_Hinf_gain_2D',              'PIETOOLS_Hinf_gain(PIE,settings,''P'',gain)';
     'PIETOOLS_Hinf_gain_2D_non_coercive', 'PIETOOLS_Hinf_gain(PIE,settings,''Q'',gain)';
     'PIETOOLS_Hinf_gain_dual_2D',         'PIETOOLS_Hinf_gain(PIE,settings,''Pd'',gain)';
     'PIETOOLS_H2_norm_2D_c',              'PIETOOLS_H2_norm(PIE,settings,''cco'')';
     'PIETOOLS_H2_norm_2D_c_non_coercive', 'PIETOOLS_H2_norm(PIE,settings,''c'')';
     'PIETOOLS_H2_norm_2D_o',              'PIETOOLS_H2_norm(PIE,settings,''oco'')';
     'PIETOOLS_H2_norm_2D_o_non_coercive', 'PIETOOLS_H2_norm(PIE,settings,''o'')';
     'PIETOOLS_Hinf_estimator_2D',         'PIETOOLS_Hinf_estimator(PIE,settings,gain)'};
for k = 1:size(T,1)
    w = which(T{k,1},'-all');
    chk(isscalar(w) && contains(strrep(w{1},'\','/'),'/executives/retired_2025/'), ...
        sprintf('(1) %s resolves %d time(s)',T{k,1},numel(w)));
    e = raised(@() feval(T{k,1},[],[]));
    chk(strcmp(e.identifier,[T{k,1} ':retired']) && contains(e.message,T{k,2}), ...
        sprintf('(1) %s: "%s" "%s"',T{k,1},e.identifier,e.message));
    n = n+1;
end
fprintf('  passed: (1) %d executives of PIETOOLS 2025 stop with the call of this release\n',size(T,1));
% (2) the TOL of getController, getObserver
for f = {'getController','getObserver'}
    e = raised(@() feval(f{1},1,1,1e-3));
    chk(strcmp(e.identifier,[f{1} ':call2025']),sprintf('(2) %s(P,Z,1e-3): "%s"',f{1},e.identifier));
    for o = {[], struct()}
        e = raised(@() feval(f{1},1,1,o{1}));
        chk(~strcmp(e.identifier,[f{1} ':call2025']),sprintf('(2) %s with OPTS %s raised call2025',f{1},class(o{1})));
    end
    n = n+1;
end
fprintf('  passed: (2) getController, getObserver: a numeric third input stops; [] and a struct pass the guard\n');
% (3), (4) on T'*T of the heat equation on [0,1]
pvar s t
x = pde_var(s,[0,1]);
PIE = convert(initialize([diff(x,t) == diff(x,s,2); subs(x,s,0) == 0; subs(x,s,1) == 0]));
TT = PIE.T'*PIE.T;
prog = lpiprogram(PIE.vars,PIE.dom);
for f = {'pure','sep'}
    e = raised(@() lpi_ineq(prog,TT,struct(f{1},1)));
    chk(strcmp(e.identifier,'lpi_ineq:call2025'),sprintf('(3) opts.%s: "%s"',f{1},e.identifier));
    n = n+1;
end
[~,~,i0] = lpi_ineq(prog,TT,struct('psatz','none'));
chk(isequal(i0.terms.codes,0),'(3) opts.psatz = ''none'' did not pose the plain term alone');
n = n+1;
fprintf('  passed: (3) lpi_ineq: opts.pure, opts.sep stop; opts.psatz = ''none'' runs\n');
vals = {0, 1, [0 1], 'none', 'auto'};
R.psatz = struct('value',vals,'codes',[],'N',[],'ndec',[]);
for k = 1:numel(vals)
    [~,~,ik] = lpi_ineq(prog,TT,struct('psatz',vals{k}));
    R.psatz(k).codes = ik.terms.codes;  R.psatz(k).N = ik.N;  R.psatz(k).ndec = ik.ndec;
    fprintf('  (4) psatz %-6s codes %-8s Gram %-10s ndec %d\n',mat2str(vals{k}),mat2str(ik.terms.codes),mat2str(ik.N),ik.ndec);
end
chk(isequal(R.psatz(2).codes,1),'(4) opts.psatz = 1 is not read as the code of the box product');
n = n+1;
R.n = n;
fprintf('test_retired2025: %d checks passed\n',n);
end

function e = raised(f)
% The error F raises, or a struct with empty identifier and message.
e = struct('identifier','','message','');
try
    f();
catch err
    e = struct('identifier',err.identifier,'message',err.message);
end
end

function chk(c,msg)
if ~c, error('test_retired2025:fail','test_retired2025: %s',msg); end
end
