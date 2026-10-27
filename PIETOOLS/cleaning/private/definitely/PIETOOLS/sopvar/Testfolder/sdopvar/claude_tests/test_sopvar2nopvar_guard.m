%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% TEST_SOPVAR2NOPVAR_GUARD checks 'sopvar2nopvar' after the fix of
% 09/29/2026: degR is read from ZR (it was read from ZL, so the degree check
% compared degL with itself and never fired), and a basis other than the
% complete set 0:deg is rejected, as 'sdopvar2ndopvar' does since
% 08/29/2026.
%
% A 'nopvar' carries one degree per variable and the complete basis 0:deg
% on both sides, with no dummy monomial in a multiplier direction. The
% converter picks the columns of each parameter by POSITION in that basis,
% so any other basis puts coefficients on the wrong monomials or drops them.
%
%   (1) ZL and ZR of different degree in some direction -> error. ZL = 0:1,
%       ZR = 0:2 used to return a degree-1 nopvar without the s'^2 column;
%   (2) a gapped basis on either side or both            -> error. ZL = [0 2],
%       ZR = 0:2 used to return a degree-2 nopvar with the s^2 row read as
%       s^1;
%   (2b) content on a dummy monomial of a multiplier direction, possible
%       only when params are assigned past the canonical form -> error. It
%       used to be dropped (kernel error 1.45 against the definition);
%   (3) complete bases in 1, 2 and 3 variables, m ~= n included:
%       (a) the operator is preserved. Both objects are applied to the same
%           random polynomial function: the sopvar by 'apply_sopvar'
%           (Testfolder/sdopvar), the nopvar by 'apply_nopvar' below, from
%           the kernel definition in the nopvar.m header with kernels from
%           'pi_ndopvar_kernels'. The kernels are also compared cell by
%           cell against 'pi_blk_kernels' (sopvar.m header). No converter
%           is in the loop (CLAUDE.md S4);
%       (b) 'copvar2nopvar' of the 1 x 1 container returns the same object;
%       (c) the output is identical to that of the 0a66f793 version, when a
%           renamed copy of it, 'sopvar2nopvar_0a66f793', is on the path
%           (not in the repository); skipped otherwise.
%
% Initial coding MMP, 09/29/2026
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

rng(20260929);
warning('off','sopvar:noncanonicalMultiplier');
tol = 1e-9;
npass = 0;
has_ref = exist('sopvar2nopvar_0a66f793','file')==2;

%% (1), (2) bases a nopvar cannot represent
fprintf('=== rejected bases ===\n');
DEG = 'different degrees';      GAP = 'complete set 0:deg';
probes = {
    {'ZL 0:1, ZR 0:2',            {[0;1]},           {[0;1;2]},         DEG}
    {'ZL 0:2, ZR 0:1',            {[0;1;2]},         {[0;1]},           DEG}
    {'2 var, ZR 0:2 in s2',       {[0;1],[0;1]},     {[0;1],[0;1;2]},   DEG}
    {'3 var, ZL 0:2 in s3',       {[0;1],[0;1],(0:2)'}, {[0;1],[0;1],[0;1]}, DEG}
    {'ZL = ZR = [0 2]',           {[0;2]},           {[0;2]},           GAP}
    {'ZL [0 2], ZR 0:2',          {[0;2]},           {[0;1;2]},         GAP}
    {'ZL 0:2, ZR [0 2]',          {[0;1;2]},         {[0;2]},           GAP}
    {'2 var, [0 2] in s2',        {[0;1],[0;2]},     {[0;1],[0;2]},     GAP}
    {'3 var, [1] in s2',          {[0;1],1,[0;1]},   {[0;1],1,[0;1]},   GAP}
    };
for k = 1:size(probes,1)
    pr = probes{k};
    P = build_sopvar([1,1],pr{2},pr{3},0.8);
    msg = '';
    try
        sopvar2nopvar(P);
    catch ME
        msg = ME.message;
    end
    if isempty(msg)
        error('test_sopvar2nopvar_guard: %s was accepted.',pr{1});
    elseif ~contains(msg,pr{4})
        error('test_sopvar2nopvar_guard: %s raised the wrong error: %s',pr{1},msg);
    end
    fprintf('  rejects : %-22s %s\n',pr{1},shorten(msg));
    npass = npass+1;
end

%% (2b) multiplier content on a dummy monomial (params assigned directly)
for nv = 1:2
    ZZ = repmat({(0:2)'},1,nv);
    P = build_sopvar([1,1],ZZ,ZZ,0.8);
    P.params{1}(:,end) = 1;     % cell 1: multiplier in every direction
    msg = '';
    try
        sopvar2nopvar(P);
    catch ME
        msg = ME.message;
    end
    if ~contains(msg,'multiplier direction')
        error('test_sopvar2nopvar_guard: multiplier content on a dummy monomial (%d var) was not refused: %s',nv,msg);
    end
    fprintf('  rejects : %-22s %s\n',sprintf('%d var, mult. s''^2',nv),shorten(msg));
    npass = npass+1;
end

%% (3) complete bases: the operator is preserved
fprintf('\n=== complete bases: operator preserved ===\n');
% {label, [m n], degree per variable}
cases = {
    {'1D m1 n1 deg0',     [1 1], 0}
    {'1D m1 n1 deg2',     [1 1], 2}
    {'1D m2 n3 deg1',     [2 3], 1}
    {'1D m3 n1 deg3',     [3 1], 3}
    {'2D m1 n1 deg[1 1]', [1 1], [1 1]}
    {'2D m2 n2 deg[2 1]', [2 2], [2 1]}
    {'2D m1 n2 deg[0 2]', [1 2], [0 2]}
    {'3D m1 n1 deg[1 1 1]', [1 1], [1 1 1]}
    {'3D m2 n1 deg[1 2 1]', [2 1], [1 2 1]}
    };
for ic = 1:size(cases,1)
    lbl = cases{ic}{1};     dims = cases{ic}{2};    degs = cases{ic}{3};
    N = numel(degs);
    Z = arrayfun(@(d) (0:d)',degs,'UniformOutput',false);
    P = build_sopvar(dims,Z,Z,0.6);
    if ~isequal(P.ZL(:),Z(:)) || ~isequal(P.ZR(:),Z(:))
        error('test_sopvar2nopvar_guard: the constructor changed the bases (%s).',lbl);
    end
    Pn = sopvar2nopvar(P);

    names = P.vars.in(:).';
    dumnames = strcat(names,'_dum');
    pts = rand(4,2*N);

    % (a) kernels, cell by cell, each from its own class definition
    Ks = pi_blk_kernels(P,[]);
    Kn = pi_ndopvar_kernels(Pn,zeros(0,1),dumnames);
    if numel(Ks)~=numel(Kn)
        error('test_sopvar2nopvar_guard: %d sopvar cells against %d nopvar cells (%s).',...
              numel(Ks),numel(Kn),lbl);
    end
    eK = 0;     sK = 0;
    for j = 1:numel(Ks)
        eK = max(eK,max_at(Ks{j}-Kn{j},[names,dumnames],pts));
        sK = max(sK,max_at(Ks{j},[names,dumnames],pts));
    end
    % (a) the action on a random function, degree 2 in every variable
    x = rand_function(dims(2),names);
    ys = apply_sopvar(P,x);
    yn = apply_nopvar(Pn,x);
    eY = max_at(ys-yn,names,pts(:,1:N));
    sY = max_at(ys,names,pts(:,1:N));
    if ~(eK<=tol*max(1,sK)) || ~(eY<=tol*max(1,sY))
        error('test_sopvar2nopvar_guard: operator changed, kernels %.3g, action %.3g (%s).',...
              eK,eY,lbl);
    end
    npass = npass+2;

    % (b) the container wrapper
    Pn2 = copvar2nopvar(copvar({P}));
    if ~same_nopvar(Pn2,Pn)
        error('test_sopvar2nopvar_guard: copvar2nopvar differs from sopvar2nopvar (%s).',lbl);
    end
    npass = npass+1;

    % (c) unchanged from 0a66f793 on a basis both versions accept
    rtxt = 'ref not on path';
    if has_ref
        if ~same_nopvar(sopvar2nopvar_0a66f793(P),Pn)
            error('test_sopvar2nopvar_guard: output differs from 0a66f793 (%s).',lbl);
        end
        rtxt = 'identical to 0a66f793';
        npass = npass+1;
    end
    fprintf('  passed: %-21s kernels %.0e | action %.0e | copvar ok | %s\n',...
            lbl,eK,eY,rtxt);
end

if has_ref,     rsum = 'comparison with 0a66f793 run';
else,           rsum = 'comparison with 0a66f793 SKIPPED (sopvar2nopvar_0a66f793 not on the path)';
end
fprintf('\ntest_sopvar2nopvar_guard passed (%d checks; %s).\n',npass,rsum);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function P = build_sopvar(dims,ZL,ZR,dnsty)
% A random sopvar on L2[s1,...,sN] -> L2[s1,...,sN], [0,1] in every
% direction, with the given bases, in canonical multiplier form (a
% multiplier direction has content only on ZR monomials of degree 0 there),
% so the constructor leaves it as given.

N = numel(ZL);
names = arrayfun(@(k) sprintf('s%d',k),1:N,'UniformOutput',false);
vars = struct('in',{names},'out',{names});
dom = struct('in',repmat([0,1],N,1),'out',repmat([0,1],N,1));
ZL = cellfun(@(z) reshape(z,[],1),ZL,'UniformOutput',false);
ZR = cellfun(@(z) reshape(z,[],1),ZR,'UniformOutput',false);
% Degrees of the ZR monomials, first variable slowest (pi_monom_vector)
dm = zeros(1,0);
for i = 1:N
    dm = [kron(dm,ones(numel(ZR{i}),1)), kron(ones(size(dm,1),1),ZR{i})];  %#ok<AGROW>
end
NL = prod(cellfun(@numel,ZL));     NR = size(dm,1);
params = cell([3*ones(1,N),1,1]);
for k = 1:numel(params)
    g = cell(1,N);
    [g{:}] = ind2sub([3*ones(1,N),1],k);
    g = cell2mat(g);
    keep = all(dm(:,g==1)==0,2);
    keep = repmat(keep,dims(2),1);          % matrix column outer, monomial inner
    C = sprandn(dims(1)*NL,dims(2)*NR,dnsty);
    C(:,~keep) = 0;
    params{k} = C;
end
P = sopvar(params,vars,ZL,ZR,dom,dims);

end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function y = apply_nopvar(Pn,x)
% (Pn x)(s) = sum_j int_[a,b] I_j(s,t) R_j(s,t) x(t) dt, from the nopvar.m
% header: cell index 1 in a direction is the multiplier, 2 the integral
% over [a_i,s_i], 3 over [s_i,b_i]. A multiplier direction carries no dummy
% monomial, so x keeps s_i there.

N = numel(Pn.deg);
s = Pn.vars(:,1);   t = Pn.vars(:,2);
K = pi_ndopvar_kernels(Pn,zeros(0,1));
sz = [size(Pn.C),1];
y = polynomial(zeros(Pn.dim(1),size(x,2)));
for j = 1:numel(Pn.C)
    g = cell(1,N);
    [g{:}] = ind2sub(sz,j);
    g = cell2mat(g);
    xj = x;
    for i = find(g>1)
        xj = subs(xj,s(i),t(i));
    end
    trm = K{j}*xj;
    for i = find(g>1)
        if g(i)==2
            trm = int(trm,t(i),Pn.dom(i,1),s(i));
        else
            trm = int(trm,t(i),s(i),Pn.dom(i,2));
        end
    end
    y = y + trm;
end

end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function x = rand_function(n,names)
% n x 1 polynomial with random coefficients, degree <= 2 in each variable.

Zx = pi_monom_vector(repmat({(0:2)'},1,numel(names)),names);
x = randn(n,size(Zx,1))*Zx;

end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function e = max_at(p,names,pts)
% Largest absolute entry of the polynomial matrix p over the rows of pts,
% pts(r,i) being the value of variable names{i}.

e = 0;
for r = 1:size(pts,1)
    pr = subs(polynomial(p),names(:),pts(r,:).');
    e = max(e,max(abs(double(pr(:)))));
end

end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function tf = same_nopvar(A,B)
% Identical nopvar objects: coefficients, degrees, domain and variables.

tf = isequal(size(A.C),size(B.C)) && isequal(A.deg,B.deg) && isequal(A.dom,B.dom) ...
     && isequal(A.vars.varname,B.vars.varname) && isequal(A.vars.degmat,B.vars.degmat) ...
     && isequal(A.vars.coef,B.vars.coef) && isequal(A.vars.matdim,B.vars.matdim);
for k = 1:numel(A.C)
    if ~tf,     return,     end
    tf = isequal(A.C{k},B.C{k});
end

end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function s = shorten(s)
s = regexprep(strtrim(s),'\s+',' ');
if numel(s)>60,     s = [s(1:60),'...'];    end
end
