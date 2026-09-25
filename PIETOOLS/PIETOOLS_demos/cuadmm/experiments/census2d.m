% census2d.m -- PHASE 0 continued: the 2-D executives, same recipe as
% census1d.  2-D is one class among several here, NOT a privileged target;
% it is censused because breadth requires it.

cuadmm_path;
HERE = fileparts(mfilename('fullpath'));
cuadmm_shadow('lpisolve');
global CENSUS_PROG

lam = 2;
mkA = @() build2(lam,false,false);
mkB = @() build2(lam,true ,false);   % + disturbance w, regulated z
mkD = @() build2(lam,true ,true );   % + sensed output y

st = lpisettings('light');

J = { 'PIETOOLS_stability_2D'              , mkA
      'PIETOOLS_stability_dual_2D'         , mkA
      'PIETOOLS_Hinf_gain_2D'              , mkB
      'PIETOOLS_Hinf_gain_2D_non_coercive' , mkB
      'PIETOOLS_Hinf_gain_dual_2D'         , mkB
      'PIETOOLS_Hinf_estimator_2D'         , mkD
      'PIETOOLS_H2_norm_2D_c'              , mkB
      'PIETOOLS_H2_norm_2D_o'              , mkB
      'PIETOOLS_H2_norm_2D_c_non_coercive' , mkB
      'PIETOOLS_H2_norm_2D_o_non_coercive' , mkB };

fprintf('CEN2 header exec,m,nvar,Kf,nblk,Ns,eigcost,svec,nnzAt,normb,c_nnz,feas,dim_ok,t_build\n');
R2 = struct([]);
for i = 1:size(J,1)
    name = J{i,1};
    CENSUS_PROG = [];
    try
        PIE = J{i,2}();
        f  = str2func(name);
        t0 = tic;  evalc('f(PIE,st);');  tb = toc(t0);
        S  = cuadmm_private('sdpshape',CENSUS_PROG);
        fprintf('CEN2 %s,%d,%d,%d,%d,[%s],%g,%d,%d,%.4e,%d,%d,%d,%.1f\n', ...
            name,S.m,S.nvar,S.Kf,S.nblk,strtrim(num2str(S.Ks)),S.eigcost, ...
            S.svec_len,S.nnzAt,S.normb,S.c_nnz,S.feas,S.dim_ok,tb);
        S.name = name;  R2(end+1).S = S;                             %#ok<SAGROW>
    catch ME
        fprintf('CEN2 %s,FAILED,%s\n', name, strrep(ME.message,newline,' '));
    end
end
save(fullfile(cuadmm_outdir(),'census2d.mat'),'R2');
fprintf('CEN2DONE\n');


cuadmm_shadow('off');   % remove the shadow: it stays active for the whole session otherwise

function PIE = build2(lam,with_w,with_y)
% Scalar reaction-diffusion on the unit square, Dirichlet on all four edges.
pvar s1 s2 t
x = pde_var('state',1,[s1;s2],[0,1;0,1]);

rhs = diff(x,s1,2) + diff(x,s2,2) + lam*x;
if with_w
    w = pde_var('input',1);
    rhs = rhs + s1*w;
end
sys = diff(x,t,1) == rhs;
if with_w
    z = pde_var('output',1);
    sys = [sys; z == int(int(x,s1,[0,1]),s2,[0,1])];
end
if with_y
    y = pde_var('sense',1);
    sys = [sys; y == int(int(x,s1,[0,1]),s2,[0,1])];
end
sys = [sys; subs(x,s1,0)==0; subs(x,s1,1)==0;
            subs(x,s2,0)==0; subs(x,s2,1)==0];
PIE = convert(sys);
end
