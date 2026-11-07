% census1d.m -- PHASE 0: what SDP does each 1-D PIETOOLS executive actually
% pose?  No solve, no GPU.  The stub captures the assembled program and
% sdpshape reports the cone the solver would see.
%
% Four systems of increasing I/O structure, all the SAME underlying physics
% (scalar reaction-diffusion, Dirichlet, lam=2 well inside the pi^2 limit) so
% that differences between rows are attributable to the EXECUTIVE, not to the
% plant.  This is the whole point of the census: isolate what each analysis
% question costs.

cuadmm_path;
HERE = fileparts(mfilename('fullpath'));
cuadmm_shadow('lpisolve');
global CENSUS_PROG

lam = 2;
mkA = @() build(lam,false,false,false);
mkB = @() build(lam,true ,false,false);   % + disturbance w, regulated z
mkC = @() build(lam,true ,true ,false);   % + control u
mkD = @() build(lam,true ,false,true );   % + sensed output y

st = lpisettings('light');

% executive , system
J = { 'PIETOOLS_PDEstability'            , mkA
      'PIETOOLS_PDEstability_dual'       , mkA
      'PIETOOLS_PIE2PDEstability'        , mkA
      'PIETOOLS_PIE2PDEstability_dual'   , mkA
      'PIETOOLS_well_posedness'          , mkA
      'PIETOOLS_Hinf_gain'               , mkB
      'PIETOOLS_Hinf_gain_coercive'      , mkB
      'PIETOOLS_Hinf_gain_dual'          , mkB
      'PIETOOLS_Hinf_gain_dual_coercive' , mkB
      'PIETOOLS_H2_norm_c'               , mkB
      'PIETOOLS_H2_norm_o'               , mkB
      'PIETOOLS_H2_norm_c_coercive'      , mkB
      'PIETOOLS_H2_norm_o_coercive'      , mkB
      'PIETOOLS_Hinf_control'            , mkC
      'PIETOOLS_H2_control'              , mkC
      'PIETOOLS_Hinf_estimator'          , mkD
      'PIETOOLS_H2_estimator'            , mkD };

fprintf('CEN header exec,m,nvar,Kf,nblk,Ns,eigcost,svec,nnzAt,normb,c_nnz,feas,dim_ok,t_build\n');
R = struct([]);
for i = 1:size(J,1)
    name = J{i,1};
    CENSUS_PROG = [];
    try
        PIE = J{i,2}();
        f  = str2func(name);
        t0 = tic;  evalc('f(PIE,st);');  tb = toc(t0);
        S  = cuadmm_private('sdpshape',CENSUS_PROG);
        fprintf('CEN %s,%d,%d,%d,%d,[%s],%g,%d,%d,%.4e,%d,%d,%d,%.1f\n', ...
            name,S.m,S.nvar,S.Kf,S.nblk,strtrim(num2str(S.Ks)),S.eigcost, ...
            S.svec_len,S.nnzAt,S.normb,S.c_nnz,S.feas,S.dim_ok,tb);
        S.name = name;  R(end+1).S = S;                              %#ok<SAGROW>
    catch ME
        fprintf('CEN %s,FAILED,%s\n', name, strrep(ME.message,newline,' '));
    end
end
save(fullfile(cuadmm_outdir(),'census1d.mat'),'R');
fprintf('CENDONE\n');


cuadmm_shadow('off');   % remove the shadow: it stays active for the whole session otherwise

function PIE = build(lam,with_w,with_u,with_y)
% Scalar reaction-diffusion on [0,1], Dirichlet, optionally with a distributed
% disturbance, a distributed control, and a distributed measurement.
pvar s t
x = pde_var('state',1,s,[0,1]);

% Right-hand side, assembled before any '==' so the dynamics equation is built
% once.  Never concatenate [] into a pde_struct list -- vertcat of a double
% with a pde_struct is what broke the first attempt at this census.
rhs = diff(x,s,2) + lam*x;
if with_w
    w = pde_var('input',1);
    rhs = rhs + s*w;
end
if with_u
    u = pde_var('control',1);
    rhs = rhs + u;
end
sys = diff(x,t,1) == rhs;

if with_w                                   % regulated output accompanies w
    z = pde_var('output',1);
    sys = [sys; z == int(x,s,[0,1])];
end
if with_y
    y = pde_var('sense',1);
    sys = [sys; y == int(x,s,[0,1])];
end
sys = [sys; subs(x,s,0)==0; subs(x,s,1)==0];
PIE = convert(sys);
end
