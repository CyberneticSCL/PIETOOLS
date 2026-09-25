% dumpcfgs.m -- build a 2-FACTOR experiment that separates cuADMM's two
% candidate cost drivers, and write each cell as a cuADMM problem directory.
%
% The claim under test is my own, and so far only asserted: that cuADMM's
% per-iteration cost is driven by sum(N^3) (a DENSE eigendecomposition per
% block per iteration) and not by the constraint count m.  The settings shape
% map happens to supply a near-orthogonal design:
%
%   group A  (m nearly FIXED, sum(N^3) spanning 95x)
%     2d_sep0_Dup1   m=3451  sumN3=8.0e5   nnz=203k
%     2d_sep0_Dup2   m=3615  sumN3=1.0e7   nnz=320k
%     2d_sep0_Dup3   m=4283  sumN3=7.6e7   nnz=731k
%
%   group B  (sum(N^3) FIXED, m halved)
%     2d_sep0_Dup2 vs 2d_Dup2      sumN3 1.0e7, m 3615 vs 1834
%     2d_sep0_Dup3 vs 2d_base      sumN3 7.6e7, m 4283 vs 3456
%
% Feasibility is irrelevant here -- per-iteration cost is a property of the
% data shape, not of whether a certificate exists -- but it is recorded so no
% one later mistakes these for certificates.

cuadmm_path;
HERE = fileparts(mfilename('fullpath'));
cuadmm_shadow('lpisolve');     % capture, never solve
global CENSUS_PROG

ROOT = fullfile(cuadmm_outdir(),'cu_cfg');
if ~exist(ROOT,'dir'), mkdir(ROOT); end

C = { '2d_sep0_Dup1', @() dup(sep0(b2()),1)
      '2d_sep0_Dup2', @() dup(sep0(b2()),2)
      '2d_sep0_Dup3', @()     sep0(b2())
      '2d_Dup2'     , @() dup(      b2() ,2)
      '2d_base'     , @()           b2()   };

fprintf('DC label|m|nvar|Kf|nblk|Nmax|Ns|eigcost|nnzAt|vec_len\n');
for c = 1:size(C,1)
    lab = C{c,1};
    CENSUS_PROG = [];
    PIE = mk2d(2);
    st  = C{c,2}();
    f   = @PIETOOLS_stability_2D;
    evalc('f(PIE,st);');
    prog = CENSUS_PROG;
    S = cuadmm_private('sdpshape',prog);

    Atf=[]; bf=[];
    for i=1:prog.expr.num, Atf=[Atf,prog.expr.At{i}]; bf=[bf;prog.expr.b{i}]; end

    RR = mkRR(prog);                     % processvars' permutation, rebuilt
    D.At = RR'*Atf;
    D.b  = bf/max(norm(full(bf)),eps);   % NORMALISE b, as sossolve does, so
                                         % cuADMM's /(1+||b||) means what a
                                         % reader expects (see the eppos work)
    D.c  = sparse(size(Atf,1),1);        % feasibility: c = 0
    D.K  = struct('f',S.Kf,'l',0,'q',[],'s',S.Ks);
    D.Ns = S.Ks;  D.Kf = S.Kf;
    dfile = fullfile(ROOT,[lab '.mat']);
    save(dfile,'-struct','D');
    info = dump2cuadmm(dfile, fullfile(ROOT,lab));

    fprintf('DC %s|%d|%d|%d|%d|%d|[%s]|%g|%d|%d\n', lab,S.m,S.nvar,S.Kf, ...
        S.nblk,S.Nmax,strtrim(num2str(S.Ks)),S.eigcost,info.nnz_At,info.vec_len);
end
fprintf('DCDONE\n');


cuadmm_shadow('off');   % remove the shadow: it stays active for the whole session otherwise

function RR = mkRR(prog)
% Rebuild sossolve/processvars' RR: identity on the decision variables, then
% ANTI-block-diagonal for each 'poly' variable and block-diagonal for each
% 'sos'/extra variable.  Assumed identity elsewhere in this work; built
% properly here so the dumped At is the one the solver would actually see.
RR = speye(prog.var.idx{1}-1);
for i = 1:prog.var.num
    sz = prog.var.idx{i+1}-prog.var.idx{i};
    switch prog.var.type{i}
        case 'poly', RR = spantiblkdiag(RR,speye(sz));
        case 'sos',  RR = spblkdiag(RR,speye(sz));
    end
end
for i = 1:prog.extravar.num
    RR = spblkdiag(RR,speye(prog.extravar.idx{i+1}-prog.extravar.idx{i}));
end
end

function st = b2()
st = lpisettings('light');
st.settings_2d.eppos = 1e-2*[1;1;1;1];  st.eppos = 1e-2;  st.eppos2 = 1e-2;
end

function st = sep0(st)
st.settings_2d.LF_opts.sep = zeros(6,1);
end

function st = dup(st,Dup)
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

function PIE = mk2d(lam)
pvar s1 s2 t
x = pde_var('state',1,[s1;s2],[0,1;0,1]);
PIE = convert([diff(x,t,1)==diff(x,s1,2)+diff(x,s2,2)+lam*x;
               subs(x,s1,0)==0; subs(x,s1,1)==0;
               subs(x,s2,0)==0; subs(x,s2,1)==0]);
end
