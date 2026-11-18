function R = pf_poinc(tag,name,PIE,degs,optl,nbis,hiMul)
% R = PF_POINC(TAG,NAME,PIE,DEGS,OPTL,NBIS,HIMUL) the mixed-domain analogue
% of the Poincare operator: X(gam) = gam I - T'*T on the fundamental spaces
% of the PIE, strictly positive for gam > ||T||^2. lam = a Rayleigh lower
% bound of ||T||^2 over Legendre products orthonormal on the PIE's own
% domain (degrees DEGS, one entry per space); then for each option struct
% in OPTL a bisection of NBIS steps on gam in [lam, HIMUL lam] (HIMUL
% default 4) of the certificate X(gam) >= 0 by LPI_INEQ_SOP (MOSEK optimal
% = certified). Prints lam, the certified reach gam_min/lam and the last
% failing status. Results in results/poinc_<TAG>.mat.
% SCRATCH, MMP, 11/16/2026; domains of the PIE (not [0,1]) the same day.
if nargin<7 || isempty(hiMul),  hiMul = 4;  end
S = pf_setup(false);
T = PIE.T;
st = exec_sop_settings(lpisettings_sop('light'));   tl = exec_tools_sop(PIE,st);
I = tl.eye(2,1);                                     % identity on T's input spaces
[sp,~] = copvar_space_list(T,'in');
vars = T.vars;  dom = T.dom;                          % registry names and intervals
pv = struct();  lo = struct();  hi = struct();
for k = 1:numel(vars)
    pv.(vars{k}) = polynomial(1,1,vars(k),[1 1]);     % pvar cannot assign inside a function with a nested function
    lo.(vars{k}) = dom(k,1);    hi.(vars{k}) = dom(k,2);
end
E = {};
for k = 1:numel(sp)
    v = sp{k};
    if isempty(v),  e = {polynomial(1)};
    elseif numel(v)==1
        e = cell(1,degs(k)+1);
        for m = 0:degs(k),  e{m+1} = shleg(m,pv.(v{1}),lo.(v{1}),hi.(v{1}));   end
    else
        e = {};
        for m = 0:degs(k), for n = 0:degs(k)
            e{end+1} = shleg(m,pv.(v{1}),lo.(v{1}),hi.(v{1}))*shleg(n,pv.(v{2}),lo.(v{2}),hi.(v{2})); %#ok<AGROW>
        end, end
    end
    for i = 1:numel(e)
        c = repmat({polynomial(0)},1,numel(sp));   c{k} = e{i};   E{end+1} = c;    %#ok<AGROW>
    end
end
nE = numel(E);  Y = cell(1,nE);
for j = 1:nE,   Y{j} = apply_cop(T,E{j});   end
Mg = zeros(nE);
for i = 1:nE, for j = i:nE,  Mg(i,j) = ipc(Y{i},Y{j},sp,lo,hi);  Mg(j,i) = Mg(i,j);  end, end
lam = max(eig((Mg+Mg')/2));
fprintf('@@ POINC %s lam=%.6f (Rayleigh over %d functions; domain %s)\n',name,lam,nE,mat2str(dom));
out = sprintf('%s/results/poinc_%s.mat',S,tag);
R = struct('name',{},'opts',{},'lam',{},'gmin',{},'reach',{},'last_fail',{},'solves',{},'t',{});
for j = 1:numel(optl)
    o = optl{j};    olab = jsonencode(o);   t0 = tic;
    lo_g = lam;   hi_g = hiMul*lam;   nsol = 0;  lastfail = '';
    [ok,stat] = cert(gamop(hi_g),o,st);   nsol = nsol+1;
    if ~ok,     lastfail = stat;    hi_g = NaN;
    else
        for b = 1:nbis
            mid = (lo_g+hi_g)/2;
            [ok,stat] = cert(gamop(mid),o,st);    nsol = nsol+1;
            if ok,  hi_g = mid;   else,   lo_g = mid;   lastfail = stat;    end
        end
    end
    r = struct('name',name,'opts',olab,'lam',lam,'gmin',hi_g,'reach',hi_g/lam,'last_fail',lastfail,'solves',nsol,'t',toc(t0));
    R(end+1) = r;                                                           %#ok<AGROW>
    fprintf('@@ REACH %s opts=%s lam=%.6f gmin=%.6f reach=%.5f last_fail=%s solves=%d t=%.0f\n',name,olab,lam,hi_g,hi_g/lam,lastfail,nsol,toc(t0));
    save(out,'R');
end
fprintf('## POINC DONE %s\n',tag);
    function X = gamop(g)
        X = g*I - T'*T;
    end
end

function [ok,stat] = cert(X,o,st)
% certificate of X >= 0 by lpi_ineq_sop with options o; ok = MOSEK optimal
try
    prog = lpiprogram(X.vars(:),[],X.dom,st);
    [prog,~,~] = lpi_ineq_sop(prog,X,o);
    prog = lpisolve(prog,st.sos_opts);
    cl = lpi_classify_sop(prog,[]);
    stat = cl.status;   ok = strcmp(stat,'optimal');
catch ME
    stat = ['error: ' strrep(ME.message,newline,' ')];    ok = false;
end
end

function L = shleg(n,s,a,b)
% Legendre polynomial of degree n on [a,b], orthonormal there
u = (2*s-a-b)/(b-a);  L0 = polynomial(1);  L1 = u;
if n==0,    L = L0/sqrt(b-a);   return;     end
for m = 1:n-1,  L2 = ((2*m+1)*u*L1 - m*L0)/(m+1);   L0 = L1;   L1 = L2;    end
L = L1*sqrt((2*n+1)/(b-a));
end

function y = apply_cop(P,v)
y = cell(size(P.C,1),1);
for i = 1:size(P.C,1)
    acc = polynomial(zeros(P.dim_out(i),1));
    for j = 1:size(P.C,2)
        if isempty(P.C{i,j}),   continue,   end
        acc = acc + apply_sopvar(P.C{i,j},polynomial(v{j}));
    end
    y{i} = acc;
end
end

function v = ipc(x,y,sp,lo,hi)
% L2 inner product of the state cells x, y over the spaces sp on the PIE's intervals
v = 0;
for k = 1:numel(sp)
    p = polynomial(x{k})*polynomial(y{k});
    for d = 1:numel(sp{k})
        nm = sp{k}{d};
        p = int(p,polynomial(1,1,{nm},[1 1]),lo.(nm),hi.(nm));
    end
    v = v + double(p);
end
end
