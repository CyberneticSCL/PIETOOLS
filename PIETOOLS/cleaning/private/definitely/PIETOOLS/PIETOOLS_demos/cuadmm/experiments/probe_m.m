% probe_m.m -- shape-only pricing of heavy-family variants, ranked by m.
% The cost model measurement says cuADMM's t/iter ~ m^1.10 and is nearly flat
% in sum(N^3), so the tuning objective is MINIMUM m SUBJECT TO CERTIFYING.
% heavy is the only configuration that certifies so far, at m=8950; heavy turns
% LF_opts.sep OFF, so turning it back on is the main m lever available.
cuadmm_path;
HERE = fileparts(mfilename('fullpath'));
cuadmm_shadow('lpisolve');
global CENSUS_PROG
V = { 'heavy'                , @()                      bh()
      'heavy+LFsep1'         , @()            sep1(     bh())
      'heavy+Dup2'           , @()        dup(          bh() ,2)
      'heavy+Dup1'           , @()        dup(          bh() ,1)
      'heavy+LFsep1+Dup2'    , @()        dup(sep1(     bh()),2)
      'heavy+LFsep1+Dup1'    , @()        dup(sep1(     bh()),1)
      'light'                , @()                      bl() };
fprintf('PM label|m|nblk|Nmax|Ns|eigcost|nnzAt|pred_rel_titer\n');
base_m = [];
for v = 1:size(V,1)
    CENSUS_PROG = [];
    try
        PIE = mk2d(2);  st = V{v,2}();
        f = @PIETOOLS_stability_2D;  evalc('f(PIE,st);');
        S = cuadmm_private('sdpshape',CENSUS_PROG);
        if isempty(base_m), base_m = S.m; end
        % predicted per-iteration cost RELATIVE to heavy, from the measured
        % model t/iter ~ m^1.10 * (sumN3)^0.035
        pr = (S.m/base_m)^1.10;
        fprintf('PM %s|%d|%d|%d|[%s]|%g|%d|%.3f\n', V{v,1},S.m,S.nblk,S.Nmax, ...
                strtrim(num2str(S.Ks)),S.eigcost,S.nnzAt,pr);
    catch ME
        fprintf('PM %s|FAILED|%s\n',V{v,1},strrep(ME.message,newline,' '));
    end
end
fprintf('PMDONE\n');

cuadmm_shadow('off');   % remove the shadow: it stays active for the whole session otherwise

function st = bh(), st = lpisettings('heavy'); st = epp(st); end
function st = bl(), st = lpisettings('light'); st = epp(st); end
function st = epp(st)
st.settings_2d.eppos = 1e-2*[1;1;1;1]; st.eppos=1e-2; st.eppos2=1e-2;
st.sos_opts.solver = 'mosek';
end
function st = sep1(st), st.settings_2d.LF_opts.sep = ones(6,1); end
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
