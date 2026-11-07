% probe_psatz.m -- price the Positivstellensatz variants by SHAPE before
% spending solve time.  eq_use_psatz is [0;0] in BOTH light and heavy, and a
% psatz term on the negativity constraint is the standard way to strengthen
% this relaxation -- the one untested axis that could raise reach rather than
% trade it.  It costs m, which the measured cost model says is what matters.
cuadmm_path;
HERE = fileparts(mfilename('fullpath'));
cuadmm_shadow('lpisolve');
global CENSUS_PROG
V = { 'heavy'                 , @() bh()
      'heavy+eqpsatz1'        , @() eqp(bh(),1)
      'heavy+eqpsatz2'        , @() eqp(bh(),2)
      'heavy+eqpsatz12'       , @() eqp(bh(),[1;2])
      'heavy+LFpsatz1'        , @() lfp(bh(),1)
      'heavy+Dup2'            , @() dup(bh(),2)
      'heavy+Dup2+eqpsatz1'   , @() eqp(dup(bh(),2),1)
      'heavy+Dup2+eqpsatz12'  , @() eqp(dup(bh(),2),[1;2])
      'heavy+Dup1+eqpsatz12'  , @() eqp(dup(bh(),1),[1;2]) };
fprintf('PP label|m|nblk|Nmax|Ns|eigcost|nnzAt|pred_titer_vs_heavy\n');
m0 = [];
for v = 1:size(V,1)
    CENSUS_PROG = [];
    try
        PIE = mk2d(2); st = V{v,2}();
        f = @PIETOOLS_stability_2D; evalc('f(PIE,st);');
        S = cuadmm_private('sdpshape',CENSUS_PROG);
        if isempty(m0), m0 = S.m; e0 = S.eigcost; end
        pr = (S.m/m0)^1.10 * (S.eigcost/e0)^0.035;
        fprintf('PP %s|%d|%d|%d|[%s]|%g|%d|%.3f\n', V{v,1},S.m,S.nblk,S.Nmax, ...
                strtrim(num2str(S.Ks)),S.eigcost,S.nnzAt,pr);
    catch ME
        fprintf('PP %s|FAILED|%s\n',V{v,1},strrep(ME.message,newline,' '));
    end
end
fprintf('PPDONE\n');

cuadmm_shadow('off');   % remove the shadow: it stays active for the whole session otherwise

function st = bh()
st = lpisettings('heavy');
st.settings_2d.eppos = 1e-2*[1;1;1;1]; st.eppos=1e-2; st.eppos2=1e-2;
st.sos_opts.solver = 'mosek';
end
function st = eqp(st,v)
% Set eq_use_psatz AND sync eq_opts_psatz{j}.psatz, which the settings file
% builds as eq_opts_psatz{j}.psatz = eq_use_psatz(j).
v = v(:);
st.settings_2d.eq_use_psatz = v;
for j = 1:numel(v)
    st.settings_2d.eq_opts_psatz{j}.psatz = v(j);
end
end
function st = lfp(st,v)
st.settings_2d.LF_use_psatz = v;
for j = 1:numel(v)
    st.settings_2d.LF_opts_psatz{j}.psatz = v(j);
end
end
function st = dup(st,Dup)
s2 = st.settings_2d;
dx=s2.LF_deg.dx; dy=s2.LF_deg.dy; d2=s2.LF_deg.d2;
s2.eq_deg.dx = {Dup+dx{1};  Dup+dx{2};  Dup+dx{3}};
s2.eq_deg.dy = {Dup+dy{1},  Dup+dy{2},  Dup+dy{3}};
s2.eq_deg.d2 = cellfun(@(c) Dup+c, d2, 'UniformOutput', false);
if isfield(s2,'eq_deg_psatz')
  for j=1:numel(s2.eq_deg_psatz)
    s2.eq_deg_psatz{j}.dx = {Dup-1+dx{1}; Dup-1+dx{2}; Dup-1+dx{3}};
    s2.eq_deg_psatz{j}.dy = {Dup-1+dy{1}, Dup-1+dy{2}, Dup-1+dy{3}};
    s2.eq_deg_psatz{j}.d2 = cellfun(@(c) Dup-1+c, d2, 'UniformOutput', false);
  end
end
st.settings_2d = s2;
end
function PIE = mk2d(lam)
pvar s1 s2 t
x = pde_var('state',1,[s1;s2],[0,1;0,1]);
PIE = convert([diff(x,t,1)==diff(x,s1,2)+diff(x,s2,2)+lam*x;
               subs(x,s1,0)==0; subs(x,s1,1)==0;
               subs(x,s2,0)==0; subs(x,s2,1)==0]);
end
