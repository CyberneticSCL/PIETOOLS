% Test_quadPoly_combine.m -- verify quadPoly/combine, whose local function
% combine_one_side builds its index maps (old full tensor index -> sorted
% -> deduplicated) for all indices at once since 11/25/2026, instead of
% one index per call of lin2subv/subv2lin, and treats an empty exponent
% vector as [0] (the constant monomial, as basisSize does).
%
% Each result is compared with (1) the code before 11/25/2026 (local
% function combine_old) on the same object with empty exponent vectors
% replaced by 0 (the same object if there are none), whose output must be
% the same object, bit for bit; (2) the polynomial: F(s,t) =
% (I_m x Z(s))'*C*(I_n x Z(t)) is the same before and after at random
% points (names all nonempty; an empty exponent vector as [0]); (3) the
% canonical form: names sorted (unnamed last), degrees sorted and unique,
% no variable of degree 0 only, C of size m*d_s x n*d_t. Covers 0 to 6
% variables per side, repeated and unsorted degrees, constant-only
% variables, empty exponent vectors, repeated, empty and missing names,
% row and column cells, 1 to 3 rows and columns, zero and empty C
% (Part 1); 4 variables of 18 degrees per side (104976 tensor indices),
% against the old code, and the time of the new code, required < 0.5 s
% (Part 2); and the two failures of the old code on empty exponent
% vectors: an index error if the other variables give more than one
% monomial, and all variables of the side dropped otherwise (Part 3).
%
% DJ (CC), 11/25/2026: Initial coding;

rng(20261125);

%% Part 1: against the old code, the polynomial, and the canonical form
nobj = 3000;
neval = 0;
for c=1:nobj
    is_named = rand<0.5;
    P = rand_quadPoly_test(is_named);
    label = sprintf('object %d',c);
    E = combine(P);
    R = combine_old(empty2zero(P));
    chk(isequal(getByteStreamFromArray({E.C,E.Zs,E.Zt,E.dim,E.ns,E.nt}), ...
        getByteStreamFromArray({R.C,R.Zs,R.Zt,R.dim,R.ns,R.nt})),[label,', vs the old code']);
    check_canonical(E,label);
    if is_named
        [svals,tvals] = deal(0.5+rand(1,5),0.5+rand(1,5));
        F_in = eval_quadPoly(P,svals,tvals);
        F_out = eval_quadPoly(E,svals,tvals);
        chk(norm(F_in-F_out,inf)<=1e-12*max(1,norm(F_in,inf)),[label,', polynomial']);
        neval = neval+1;
    end
end
fprintf('Test_quadPoly_combine Part 1 (%d objects, %d evaluated) passed.\n',nobj,neval)

%% Part 2: 104976 tensor indices per side
Zs = {randi([0,9],18,1),randi([0,9],18,1),(17:-1:0)',randi([0,3],18,1)};
Zt = {(0:17)',randi([0,9],18,1),zeros(18,1),randi([0,9],18,1)};
nZ = 18^4;
C = sparse(randi(2*nZ,500,1),randi(nZ,500,1),randn(500,1),2*nZ,nZ);
P = quadPoly(C,Zs,Zt,[2,1],{'s3','s1','','s1'},{'t2','t1','t3',''});
tic;    E = combine(P);     t = toc;
R = combine_old(P);
chk(isequal(getByteStreamFromArray({E.C,E.Zs,E.Zt,E.dim,E.ns,E.nt}), ...
    getByteStreamFromArray({R.C,R.Zs,R.Zt,R.dim,R.ns,R.nt})),'104976 indices, vs the old code');
check_canonical(E,'104976 indices');
chk(t<0.5,sprintf('combine with 104976 indices per side took %.2f s, expected < 0.5 s',t));
fprintf('Test_quadPoly_combine Part 2 (104976 indices per side, %.3f s) passed.\n',t)

%% Part 3: empty exponent vectors
% Other variables with more than one monomial (old code: index error)
P = quadPoly(sparse([1;2]),{[0;1],[]},{},[1,1],{'a','b'},{});
E = combine(P);
chk(isequal(E.Zs,{[0;1]}) && isequal(E.ns,{'a'}) && isequal(E.C,sparse([1;2])),'{[0;1],[]}');
% One monomial in all (old code: Zs = {}, the polynomial b^5 became 1)
P = quadPoly(sparse(3),{[],5},{},[1,1],{'a','b'},{});
E = combine(P);
chk(isequal(E.Zs,{5}) && isequal(E.ns,{'b'}) && isequal(E.C,sparse(3)),'{[],5}');
svals = 0.5+rand(1,5);
chk(abs(eval_quadPoly(P,svals,svals)-eval_quadPoly(E,svals,svals))<=1e-12*3*1.5^5,'{[],5}, polynomial');
fprintf('Test_quadPoly_combine Part 3 (empty exponent vectors) passed.\n')

disp('Test_quadPoly_combine passed.')



%% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %%
function check_canonical(E,label)
% Check the canonical form of a combined quadPoly.
% INPUTS:
% - E:      'quadPoly' object returned by combine;
% - label:  char, description of the case for error messages.
for side = {{E.Zs,E.ns},{E.Zt,E.nt}}
    [Z,names] = deal(side{1}{:});
    chk(numel(names)==numel(Z),[label,', number of names']);
    for i=1:numel(Z)
        z = Z{i};
        chk(iscolumn(z) && isequal(z,unique(z)) && ~isequal(z,0),[label,', degrees']);
    end
    if ~isempty(names)
        key = names;
        key(cellfun(@isempty,key)) = {char(127)};
        chk(issorted(key),[label,', names sorted']);
    end
end
d_s = prod(cellfun(@numel,E.Zs));
d_t = prod(cellfun(@numel,E.Zt));
chk(issparse(E.C) && isequal(size(E.C),[E.dim(1)*d_s,E.dim(2)*d_t]),[label,', size of C']);
end


function F = eval_quadPoly(P,svals,tvals)
% Evaluate F(s,t) = (I_m x Z(s))'*C*(I_n x Z(t)) from the definition of
% quadPoly, Z(s) = Zs1(s1) x ... x Zsk(sk), Zsi(si) = si.^Zs{i}(:), with
% an empty exponent vector as [0], as in quadPoly/combine.
% INPUTS:
% - P:              'quadPoly' object, with names from rand_quadPoly_test;
% - svals, tvals:   1 x 5 arrays, the values of the names of name_pool.
% OUTPUTS:
% - F:      dim(1) x dim(2) array.
P = empty2zero(P);
Zs = 1;
for i=1:numel(P.Zs)
    Zs = kron(Zs,svals(strcmp(name_pool,P.ns{i})).^P.Zs{i}(:));
end
Zt = 1;
for i=1:numel(P.Zt)
    Zt = kron(Zt,tvals(strcmp(name_pool,P.nt{i})).^P.Zt{i}(:));
end
m = P.dim(1);   n = P.dim(2);
C = P.C;
if isempty(C)
    C = sparse(m*numel(Zs),n*numel(Zt));
end
F = full(kron(eye(m),Zs)'*C*kron(eye(n),Zt));
end


function P = empty2zero(P)
% Replace the empty exponent vectors of P by 0 (the constant monomial).
% INPUTS:
% - P:      'quadPoly' object.
% OUTPUTS:
% - P:      'quadPoly' object, the same polynomial.
if iscell(P.Zs)
    P.Zs(cellfun(@isempty,P.Zs)) = {0};
end
if iscell(P.Zt)
    P.Zt(cellfun(@isempty,P.Zt)) = {0};
end
end


function names = name_pool()
% Names of the random variables (5, as for svals and tvals).
names = {'s1','s2','a','b','x<1>'};
end


function P = rand_quadPoly_test(is_named)
% Random quadPoly: per side 0 to 6 variables with 1 to 5 random degrees in
% 0:4 (repeated, unsorted; some all zero; some none, an empty exponent
% vector, counted as one monomial), names from name_pool, repeated;
% if not is_named, also empty and missing names; row or column cells;
% 1 to 3 rows and columns; C random sparse, zero, or empty.
% INPUTS:
% - is_named:   logical, true for nonempty names only.
% OUTPUTS:
% - P:      'quadPoly' object.
[Zs,ns] = rand_side(is_named);
[Zt,nt] = rand_side(is_named);
m = randi(3);   n = randi(3);
nr = m*prod(max(cellfun(@numel,Zs),1));
nc = n*prod(max(cellfun(@numel,Zt),1));
r = rand;
if r<0.05
    C = sparse(nr,nc);
elseif r<0.08
    C = [];
else
    nz = randi([0,min(200,nr*nc)]);
    C = sparse(randi(nr,nz,1),randi(nc,nz,1),round(10*randn(nz,1))/4,nr,nc);  % some cancel when summed
end
P = quadPoly(C,Zs,Zt,[m,n],ns,nt);
end


function [Z,names] = rand_side(is_named)
% Random exponent cell and names for one side of rand_quadPoly_test.
k = randi([0,6]);
Z = cell(1,k);
for i=1:k
    len = randi(5);
    r = rand;
    if r<0.04
        z = [];
    elseif r<0.19
        z = zeros(len,1);
    else
        z = randi([0,4],len,1);
    end
    if rand<0.3
        z = z.';
    end
    Z{i} = z;
end
pool = name_pool;
if is_named
    names = pool(randi(numel(pool),1,k));
else
    pool = [pool,{'',''}];
    r = rand;
    if r<0.2
        names = {};
    elseif r<0.3
        names = cell(1,k);
    else
        names = pool(randi(numel(pool),1,k));
    end
end
if rand<0.2
    Z = Z(:);
    names = names(:);
end
end


function obj = combine_old(obj)
% quadPoly/combine before 11/25/2026 (index maps one element at a time),
% its code verbatim without the comments, as reference.
% INPUTS:
% - obj:    'quadPoly' object.
% OUTPUTS:
% - obj:    'quadPoly' object with sorted names and sorted, unique degrees.
if isempty(obj.Zs), obj.Zs = {}; end
if isempty(obj.Zt), obj.Zt = {}; end
if isempty(obj.ns), obj.ns = cell(1, numel(obj.Zs)); end
if isempty(obj.nt), obj.nt = cell(1, numel(obj.Zt)); end
m = obj.dim(1);
n = obj.dim(2);
Zs = obj.Zs;
Zt = obj.Zt;
ds = basisSize_old(Zs);
dt = basisSize_old(Zt);
[Zs_new, ns_new, mapS, ds_new] = combine_one_side_old(obj.Zs, obj.ns);
[Zt_new, nt_new, mapT, dt_new] = combine_one_side_old(obj.Zt, obj.nt);
if isempty(obj.C) || nnz(obj.C) == 0
    obj.Zs = Zs_new; obj.ns = ns_new;
    obj.Zt = Zt_new; obj.nt = nt_new;
    obj.C  = sparse(m*ds_new, n*dt_new);
    return;
end
C = sparse(obj.C);
[r, c, v] = find(C);
rIn  = ceil(r(:) / ds);
sMon = mod(r(:) - 1, ds) + 1;
sMon2  = mapS(sMon);
cIn  = ceil(c(:) / dt);
tMon = mod(c(:) - 1, dt) + 1;
tMon2  = mapT(tMon);
r2 = sMon2 + (rIn - 1)*ds_new;
c2 = tMon2 + (cIn - 1)*dt_new;
C2 = sparse(r2, c2, v(:), m*ds_new, n*dt_new);
obj.Zs = Zs_new; obj.ns = ns_new;
obj.Zt = Zt_new; obj.nt = nt_new;
obj.C  = C2;
end


function [Z_new, names_new, mapFull, d_new] = combine_one_side_old(Z, names)
% combine_one_side of quadPoly/combine before 11/25/2026, its code
% verbatim without the comments.
% INPUTS:
% - Z:      1 x k cell of exponent vectors;
% - names:  1 x k cellstr of variable names, or empty.
% OUTPUTS:
% - Z_new, names_new, mapFull, d_new: see quadPoly/combine.
if isempty(Z)
    Z_new = {};
    names_new = {};
    mapFull = 1;
    d_new = 1;
    return;
end
k = numel(Z);
if isempty(names)
    names = cell(1,k);
    for i = 1:k, names{i} = ''; end
end
key = names;
for i = 1:k
    if isempty(key{i})
        key{i} = char(127);
    end
end
[~, perm] = sort(key);
Zsrt = Z(perm);
nsrt = names(perm);
sizes_srt = cellfun(@numel, Zsrt);
d_srt = prod(sizes_srt);
sizes_old = cellfun(@numel, Z);
mapPerm = zeros(prod(sizes_old), 1);
if isempty(mapPerm)
    Z_new = {};
    names_new = {};
    mapFull = 1;
    d_new = 1;
    return;
end
for idx = 1:prod(sizes_old)
    subs_old = lin2subv(idx, sizes_old);
    subs_srt = subs_old(perm);
    mapPerm(idx) = subv2lin(subs_srt, sizes_srt);
end
Z_new = cell(1,k);
names_new = nsrt;
posMap = cell(1,k);
sizes_new = zeros(1,k);
for i = 1:k
    zi = Zsrt{i}(:);
    if isempty(zi)
        zi = 0;
    end
    [zu, ~, ic] = unique(zi, 'sorted');
    Z_new{i} = zu;
    sizes_new(i) = numel(zu);
    pm = ic(:);
    posMap{i} = pm;
end
d_new = prod(sizes_new);
mapDedup = zeros(d_srt, 1);
for idx = 1:d_srt
    subs = lin2subv(idx, sizes_srt);
    subs2 = zeros(1,k);
    for i = 1:k
        subs2(i) = posMap{i}(subs(i));
    end
    mapDedup(idx) = subv2lin(subs2, sizes_new);
end
mapFull = mapDedup(mapPerm);
keep = true(1,k);
for i = 1:k
    zi = Z_new{i}(:);
    if numel(zi) == 1 && zi(1) == 0
        keep(i) = false;
    end
end
Z_new = Z_new(keep);
names_new = names_new(keep);
end


function subs = lin2subv(idx, dims)
% Subscripts of the linear index idx in an array of size dims, the last
% dimension fastest (quadPoly/combine before 11/25/2026, verbatim).
k = numel(dims);
subs = zeros(1,k);
r = idx - 1;
for i = 1:k-1
    stride = prod(dims(i+1:end));
    subs(i) = floor(r / stride) + 1;
    r = r - (subs(i)-1)*stride;
end
subs(k) = r + 1;
end


function idx = subv2lin(subs, dims)
% Linear index of the subscripts subs in an array of size dims, the last
% dimension fastest (quadPoly/combine before 11/25/2026, verbatim).
k = numel(dims);
idx = 1;
for i = 1:k-1
    stride = prod(dims(i+1:end));
    idx = idx + (subs(i)-1)*stride;
end
idx = idx + (subs(k)-1);
end


function d = basisSize_old(Zcell)
% Number of monomials of the tensor basis of the exponent vectors Zcell
% (quadPoly/combine, verbatim).
if isempty(Zcell), d = 1; return; end
d = 1;
for i = 1:numel(Zcell)
    zi = Zcell{i};
    if isempty(zi), zi = 0; end
    d = d * numel(zi);
end
end


function chk(cond,msg)
% The condition must hold
if ~cond
    error("Test_quadPoly_combine failed: %s.",msg)
end
end
