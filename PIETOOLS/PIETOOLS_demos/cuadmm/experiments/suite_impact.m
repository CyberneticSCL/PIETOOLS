% suite_impact.m -- what does the 1-D settings change do to the BENCHMARK SET?
%
% The 1-D recommendation is Dup = 2 on the negativity degrees (dd2/dd3),
% keeping the stock product psatz: measured reach 0.50 -> 0.9999 with P=I.
% Every 1-D benchmark's SDP shape therefore changes, and since cuADMM's cost
% is t/iter ~ m^1.10 * sum(N^3)^0.035 (measured), the relevant question is what
% happens to m across the whole suite -- not just on the stability executive
% the reach was measured on.
cuadmm_path;
HERE = fileparts(mfilename('fullpath'));
cuadmm_shadow('lpisolve');
global CENSUS_PROG

lam = 2;
mkA = @() build(lam,false,false,false);
mkB = @() build(lam,true ,false,false);
mkC = @() build(lam,true ,true ,false);
mkD = @() build(lam,true ,false,true );

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

fprintf('SI exec|m1|m2|m_ratio|eig1|eig2|nnz1|nnz2|Nmax1|Nmax2|pred_titer_ratio\n');
for i = 1:size(J,1)
    S = cell(1,2);
    ok = true;
    for k = 1:2
        CENSUS_PROG = [];
        try
            PIE = J{i,2}();
            st  = lpisettings('light');  st.eppos2 = 1e-2;  st.eppos = 1e-2;
            if k==2, st = dup(st,2); end
            f = str2func(J{i,1});  evalc('f(PIE,st);');
            S{k} = cuadmm_private('sdpshape',CENSUS_PROG);
        catch ME
            fprintf('SI %s|FAILED(Dup=%d)|%s\n',J{i,1},k,strrep(ME.message,newline,' '));
            ok = false; break
        end
    end
    if ~ok, continue; end
    pr = (S{2}.m/S{1}.m)^1.10 * (S{2}.eigcost/max(S{1}.eigcost,eps))^0.035;
    fprintf('SI %s|%d|%d|%.2f|%g|%g|%d|%d|%d|%d|%.2f\n', J{i,1}, ...
        S{1}.m,S{2}.m,S{2}.m/S{1}.m, S{1}.eigcost,S{2}.eigcost, ...
        S{1}.nnzAt,S{2}.nnzAt, S{1}.Nmax,S{2}.Nmax, pr);
end
fprintf('SIDONE\n');

cuadmm_shadow('off');   % remove the shadow: it stays active for the whole session otherwise

function st = dup(st,Dup)
n1=1;n2=1;n3=1;n4=n2+n3;
st.dd2 = {n1+Dup,   [n2+Dup-1, n3+Dup,   n4+Dup  ], [n2+Dup-1, n3+Dup,   n4+Dup  ]};
st.dd3 = {n1+Dup-1, [n2+Dup-2, n3+Dup-1, n4+Dup-1], [n2+Dup-2, n3+Dup-1, n4+Dup-1]};
end

function PIE = build(lam,with_w,with_u,with_y)
pvar s t
x = pde_var('state',1,s,[0,1]);
rhs = diff(x,s,2) + lam*x;
if with_w, w = pde_var('input',1); rhs = rhs + s*w; end
if with_u, u = pde_var('control',1); rhs = rhs + u; end
sys = diff(x,t,1) == rhs;
if with_w, z = pde_var('output',1); sys = [sys; z == int(x,s,[0,1])]; end
if with_y, y = pde_var('sense',1);  sys = [sys; y == int(x,s,[0,1])]; end
sys = [sys; subs(x,s,0)==0; subs(x,s,1)==0];
PIE = convert(sys);
end
