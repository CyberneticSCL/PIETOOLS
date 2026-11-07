% t4_mirror.m -- is stab_mirror faithful to the stock executive?
% Compares the SDP each poses, field by field.  A mirror that differs in m,
% Kf, block sizes or nnz is not measuring the same problem and must not be
% used to judge residuals.
cuadmm_path;
HERE = fileparts(mfilename('fullpath'));
cuadmm_shadow('lpisolve');
global CENSUS_PROG

pvar s t
x = pde_var('state',1,s,[0,1]);
PIE = convert([diff(x,t,1)==diff(x,s,2)+2*x; subs(x,s,0)==0; subs(x,s,1)==0]);

for pre = {'stripped','light','heavy'}
    st = lpisettings(pre{1});
    CENSUS_PROG = [];
    evalc('PIETOOLS_PIE2PDEstability(PIE,st);');
    A = cuadmm_private('sdpshape',CENSUS_PROG);
    M = cuadmm_private('sdpshape',getfield(cuadmm_private('stab_mirror',PIE,st),'prog'));                 %#ok<GFLD>
    same = A.m==M.m && A.Kf==M.Kf && isequal(A.Ks,M.Ks) && A.nvar==M.nvar ...
           && A.nnzAt==M.nnzAt && abs(A.normb-M.normb)<=1e-14*max(1,A.normb);
    fprintf('T4 %-10s stock m=%d Kf=%d Ns=[%s] nnz=%d | mirror m=%d Kf=%d Ns=[%s] nnz=%d | MATCH=%d\n',...
        pre{1},A.m,A.Kf,strtrim(num2str(A.Ks)),A.nnzAt, ...
               M.m,M.Kf,strtrim(num2str(M.Ks)),M.nnzAt, same);
end
fprintf('T4DONE\n');
cuadmm_shadow('off');   % remove the shadow: it stays active for the whole session otherwise
