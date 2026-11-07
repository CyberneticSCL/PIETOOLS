% t7_sepdiag.m -- whose bug is sep=1?  Stage the pipeline so the failure point
% is identified rather than guessed: stock executive, then mirror build, then
% solve, then each part of the residual panel.
cuadmm_path;
HERE = fileparts(mfilename('fullpath'));
global CENSUS_PROG
pvar s t
x   = pde_var('state',1,s,[0,1]);
PIE = convert([diff(x,t,1)==diff(x,s,2)+2*x; subs(x,s,0)==0; subs(x,s,1)==0]);

st = lpisettings('light');
st.options1.sep=1; st.options12.sep=1; st.options2.sep=1; st.options3.sep=1;

stage = @(n,f) run_stage(n,f);

% 1. stock executive with the capture stub (build only)
cuadmm_shadow('lpisolve');
stage('stock-build', @() evalc('PIETOOLS_PIE2PDEstability(PIE,st);'));
cuadmm_shadow('off');

% 2. mirror build
Mv = [];
stage('mirror-build', @() assignin('caller','dummy',0));
try, Mv = cuadmm_private('stab_mirror',PIE,st); fprintf('T7 mirror-build OK\n');
catch ME, fprintf('T7 mirror-build FAIL %s\n',ME.message); end

% 3. stock executive, REAL solve (no stub)
stage('stock-solve', @() evalc('PIETOOLS_PIE2PDEstability(PIE,st);'));

% 4. mirror solve + panel
if ~isempty(Mv)
    ok = true;
    try
        opts.solver='sedumi'; opts.params.fid=0;
        evalc('sol = lpisolve(Mv.prog,opts);');
        fprintf('T7 mirror-solve OK\n');
    catch ME, fprintf('T7 mirror-solve FAIL %s\n',ME.message); ok=false; end
    if ok
        for nm = {'Pop','Dop','Deop'}
            try
                o = getsol_lpivar(sol, Mv.(nm{1}));
                fprintf('T7 getsol(%s) OK dim=[%s]\n', nm{1}, num2str(o.dim(:)'));
            catch ME
                fprintf('T7 getsol(%s) FAIL %s\n', nm{1}, ME.message);
            end
        end
    end
end
fprintf('T7DONE\n');

function run_stage(n,f)
try, f(); fprintf('T7 %s OK\n',n);
catch ME, fprintf('T7 %s FAIL %s\n',n,strrep(ME.message,newline,' ')); end
end
