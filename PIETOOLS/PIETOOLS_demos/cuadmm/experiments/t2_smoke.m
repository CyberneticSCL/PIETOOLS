% t2_smoke.m -- does the lpisolve capture stub work, and does sdpshape agree
% with what the solver would actually see?  Three executives: one feasibility
% (c=0) and one objective-bearing (min gam), plus the dual form.

cuadmm_path;
HERE = fileparts(mfilename('fullpath'));
cuadmm_shadow('lpisolve');     % stub must SHADOW lpisolve
fprintf('T2 lpisolve resolves to %s\n', which('lpisolve'));

global CENSUS_PROG

% ---- a plain 1-D reaction-diffusion, Dirichlet, scalar state
pvar s t
phi = pde_var('state',1,s,[0,1]);
sys = [diff(phi,t,1)==diff(phi,s,2)+2*phi;
       subs(phi,s,0)==0;
       subs(phi,s,1)==0];
PIE = convert(sys);
fprintf('T2 PIE dim %d\n', PIE.dim);

st = lpisettings('light');

runs = {'PIETOOLS_PIE2PDEstability','PIETOOLS_PDEstability'};
for i = 1:numel(runs)
    CENSUS_PROG = [];
    try
        f = str2func(runs{i});
        evalc('f(PIE,st);');                 % suppress executive chatter
        S = cuadmm_private('sdpshape',CENSUS_PROG);
        fprintf(['T2 %-30s m=%-6d nvar=%-8d Kf=%-4d nblk=%-2d Ns=[%s] ' ...
                 'nnz=%-8d |b|=%.3e c_nnz=%d dim_ok=%d\n'], ...
                runs{i}, S.m, S.nvar, S.Kf, S.nblk, num2str(S.Ks), ...
                S.nnzAt, S.normb, S.c_nnz, S.dim_ok);
    catch ME
        fprintf('T2 %-30s FAILED %s\n', runs{i}, ME.message);
    end
end
fprintf('T2DONE\n');
cuadmm_shadow('off');   % remove the shadow: it stays active for the whole session otherwise
