% shapemap2d.m -- settings sweep for the 2-D executives (maintainer request).
% Same method as shapemap1d: CPU only, no solves, shape reported against both
% cost models.
%
% What the 2-D defaults already do, and what they do not:
%   LF_opts.sep  = ones(6,1)   -- separability on the LYAPUNOV operator is ON
%   eq_opts.sep  = zeros(1,6)  -- and OFF on the negativity operator
% which is exactly the split that the 1-D solves showed to be the only feasible
% one (sep on the negativity operator gives pinf=1 even at lam=0).
%
%   Dupx = Dupy = Dup2 = 3     -- HARDCODED, against Dup=1 in 1-D.  The 1-D map
% measured Dup=3 at 17x the eig work of Dup=1, so this is the first thing to
% vary here.

cuadmm_path;
HERE = fileparts(mfilename('fullpath'));
cuadmm_shadow('lpisolve');
global CENSUS_PROG

EX = {'PIETOOLS_stability_2D', @sysA2
      'PIETOOLS_Hinf_gain_2D', @sysB2};

V = { 'default(Dup=3)'   , @() lpisettings('light')
      'Dup=1'            , @() setDup2(lpisettings('light'),1)
      'Dup=2'            , @() setDup2(lpisettings('light'),2)
      'LFsep=0'          , @() setLFsep(lpisettings('light'),0)
      'ismultiplier=1'   , @() setZop(lpisettings('light'),'ismultiplier',1)
      'isscalar=1'       , @() setZop(lpisettings('light'),'isscalar',1)
      'LFpsatz=1'        , @() setfld2(lpisettings('light'),'LF_use_psatz',1)
      'eppos=1'          , @() setepp(lpisettings('light'),1) };

fprintf('SM2 header exec,variant,m,nvar,Kf,nblk,Nmax,Ns,eigcost,nnzAt,svec,normb,t\n');
for e = 1:size(EX,1)
    for v = 1:size(V,1)
        CENSUS_PROG = [];
        try
            PIE = EX{e,2}();
            st  = V{v,2}();
            f   = str2func(EX{e,1});
            t0 = tic; evalc('f(PIE,st);'); tb = toc(t0);
            S = cuadmm_private('sdpshape',CENSUS_PROG);
            fprintf('SM2 %s,%s,%d,%d,%d,%d,%d,[%s],%g,%d,%d,%.3e,%.1f\n', ...
                EX{e,1},V{v,1},S.m,S.nvar,S.Kf,S.nblk,S.Nmax, ...
                strtrim(num2str(S.Ks)),S.eigcost,S.nnzAt,S.svec_len,S.normb,tb);
        catch ME
            fprintf('SM2 %s,%s,FAILED,%s\n',EX{e,1},V{v,1}, ...
                    strrep(ME.message,newline,' '));
        end
    end
end
fprintf('SM2DONE\n');


cuadmm_shadow('off');   % remove the shadow: it stays active for the whole session otherwise

function st = setDup2(st,Dup)
% Rebuild the negativity degrees for a different bump.  Dupx/Dupy/Dup2 are
% locals in settings_PIETOOLS_light_2D and are baked into eq_deg, so they can
% only be changed by regenerating it from LF_deg.
s2 = st.settings_2d;
dx = s2.LF_deg.dx;  dy = s2.LF_deg.dy;  d2 = s2.LF_deg.d2;
s2.eq_deg.dx = {Dup+dx{1};  Dup+dx{2};  Dup+dx{3}};
s2.eq_deg.dy = {Dup+dy{1},  Dup+dy{2},  Dup+dy{3}};
s2.eq_deg.d2 = cellfun(@(c) Dup+c, d2, 'UniformOutput', false);
if isfield(s2,'eq_deg_psatz')
    for j = 1:numel(s2.eq_deg_psatz)
        s2.eq_deg_psatz{j}.dx = {Dup-1+dx{1};  Dup-1+dx{2};  Dup-1+dx{3}};
        s2.eq_deg_psatz{j}.dy = {Dup-1+dy{1},  Dup-1+dy{2},  Dup-1+dy{3}};
        s2.eq_deg_psatz{j}.d2 = cellfun(@(c) Dup-1+c, d2, 'UniformOutput', false);
    end
end
st.settings_2d = s2;
end

function st = setLFsep(st,v)
st.settings_2d.LF_opts.sep = v*ones(6,1);
end

function st = setZop(st,f,v)
st.settings_2d.Zop_opts.(f) = v;
end

function st = setfld2(st,f,v)
st.settings_2d.(f) = v;
end

function st = setepp(st,v)
st.settings_2d.eppos = v*[1;1;1;1];
st.eppos = v;  st.eppos2 = v;
end

function PIE = sysA2()
pvar s1 s2 t
x = pde_var('state',1,[s1;s2],[0,1;0,1]);
sys = [diff(x,t,1) == diff(x,s1,2)+diff(x,s2,2)+2*x;
       subs(x,s1,0)==0; subs(x,s1,1)==0;
       subs(x,s2,0)==0; subs(x,s2,1)==0];
PIE = convert(sys);
end

function PIE = sysB2()
pvar s1 s2 t
x = pde_var('state',1,[s1;s2],[0,1;0,1]);
w = pde_var('input',1);
z = pde_var('output',1);
sys = [diff(x,t,1) == diff(x,s1,2)+diff(x,s2,2)+2*x + s1*w;
       z == int(int(x,s1,[0,1]),s2,[0,1]);
       subs(x,s1,0)==0; subs(x,s1,1)==0;
       subs(x,s2,0)==0; subs(x,s2,1)==0];
PIE = convert(sys);
end
