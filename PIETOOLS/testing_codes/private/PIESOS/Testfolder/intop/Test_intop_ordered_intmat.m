% Test_intop_ordered_intmat.m -- Tests the invariant that the intmat of an
% 'intop' is ordered, and canonicalize_intmat, which splits rows that are not.
%
% Oracle: the semantics of the intervals, chi_{[x,y]}(z)=1 if x<=z<=y, and 0
% otherwise (delta for equal indices). The integral over a row is evaluated
% by brute force: for each of the n! orderings of the (non-multiplier)
% variables, check directly that every variable lies between its limits, and
% integrate the kernel over the cell of the ordering with 'int' of polynomial
% objects. This uses neither canonicalize_intmat/enumerate_valid_orders nor
% the signed integration in apply_op.
%
% Parts
%   1. ordered rows are detected, and returned unchanged;
%   2. canonicalize_intmat on hand-checked rows, also with kernels (repeated
%      for split rows, removed for rows with an empty region);
%   3. apply_op of the functional declared from canonicalize_intmat against
%      the oracle, for hand-built rows and random valid rows, d=2..4,
%      N=1,2,3 spatial variables, kernels of several terms; the constructor
%      rejects the rows as they are;
%   4. isvalid and the constructor reject a row that is not ordered, and
%      refer to canonicalize_intmat; chkval=true skips the check;
%   5. mtimes_tensopvar(K,R) against apply_op(K,apply_op(R,w)), for a non-
%      ordered K (split with canonicalize_intmat) and 3-PI operators R whose
%      limits involve the variable.
%
% DJ (CC), 10/22/2026: Initial coding;
%
% NOTES
% - Part 3 compares per direction: for the product kernel
%   k_1(s_1)*...*k_N(s_N) and the product state v_1*...*v_N, the value is
%   the product over the directions of the sum over terms of the integrals.

% Run pietools_path_update first.
rng(20261022);
nfail = 0;
chk = @(c,msg) fprintf('%s %s\n',ternary(c,'ok  ','FAIL'),msg);

%% Part 1: ordered rows
rows_ord = {
    [0,3,0,1], [0,3,1,1], [0,4,0,3,1,4], [0,4,0,3,0,1], [3,4,0,4,2,4,0,5], ...
    [0,4,1,3,1,1], [0,5,0,3,1,1,0,5], [0,2,0,3,0,4], [2,2,0,3]};
for k=1:numel(rows_ord)
    [im,~,io] = intop.canonicalize_intmat(rows_ord{k});
    ok = io && isequal(im,rows_ord{k});
    nfail = nfail + ~ok;
    chk(ok,sprintf('ordered, unchanged: [%s]',num2str(rows_ord{k})));
end
rows_unord = {
    [0,4,1,3,0,4], [2,3,0,4,0,4], [0,4,0,4,1,2], [3,4,1,3,0,4], ...
    [3,4,3,4,0,5,0,5], [0,5,4,3,0,5,1,1], [0,5,1,3,4,5,0,5]};
for k=1:numel(rows_unord)
    [~,~,io] = intop.canonicalize_intmat(rows_unord{k},[],true);
    nfail = nfail + io;
    chk(~io,sprintf('not ordered: [%s]',num2str(rows_unord{k})));
end
% x1,x2 in [x3,x4]: x3<=x4 is not implied by the limits of the other variables,
% so that the rule using all other variables at once (not recursively) would
% call this row ordered, although the signed integral is positive for x4<x3
[~,~,io] = intop.canonicalize_intmat([3,4,3,4,0,5,0,5],[],true);
nfail = nfail + io;
chk(~io,'x1,x2 in [x3,x4] is not ordered');
% A cell input gives a cell output, for N=2
[im,~,io] = intop.canonicalize_intmat({[0,4,0,3,1,4],[0,4,1,3,0,4]},[],true);
ok = iscell(im) && iscell(io) && all(io{1}) && ~all(io{2});
nfail = nfail + ~ok;
chk(ok,'cell input: one is_ord per direction');

%% Part 2: hand-checked splits
cases = {
    [0,4,1,3,0,4],          [0,2,0,3,0,4];
    [0,4,0,4,1,2],          [0,3,0,4,0,2];
    [3,4,3,4,0,5,0,5],      [0,2,0,4,0,1,0,5; 0,4,0,1,0,2,0,5];
    [0,5,4,3,0,5,1,1],      [0,2,0,3,0,5,1,1];
    [0,5,1,3,4,5,0,5],      [0,2,0,3,0,5,0,1; 0,4,0,3,0,5,0,2; 0,2,0,4,0,5,0,3]};
for k=1:size(cases,1)
    [im,pp] = intop.canonicalize_intmat(cases{k,1},{7});
    ok = isequal(sortrows(im),sortrows(cases{k,2})) && all(cell2mat(pp)==7) && size(pp,1)==size(im,1);
    nfail = nfail + ~ok;
    chk(ok,sprintf('split of [%s] into %d row(s), kernel repeated',num2str(cases{k,1}),size(im,1)));
end
% Several rows at once: unchanged rows keep their position, split rows are
% replaced in place, rows with an empty region are removed with their kernel
M = [0,4,0,3,1,4; 0,4,1,3,0,4; 3,4,1,3,0,4; 0,4,0,3,0,1];
[im,pp] = intop.canonicalize_intmat({M},{1;2;3;4});
ok = isequal(im{1},[0,4,0,3,1,4; 0,2,0,3,0,4; 0,4,0,3,0,1]) && isequal(cell2mat(pp)',[1,2,4]);
nfail = nfail + ~ok;
chk(ok,'rows with and without a split, and with an empty region: rows and kernels');
% N=2: the kernels are indexed along the direction of the split row
[im,pp] = intop.canonicalize_intmat({[0,3,0,1; 0,3,1,1],[3,4,3,4,0,5,0,5]},{10;20});
ok = isequal(size(pp),[2,2]) && isequal(cell2mat(pp),[10,10;20,20]) && size(im{2},1)==2 && size(im{1},1)==2;
nfail = nfail + ~ok;
chk(ok,'N=2: kernels repeated along the second direction');
[im,pp] = intop.canonicalize_intmat({[3,4,3,4,0,5,0,5],[0,3,0,1; 0,3,1,1]},{10,20});
ok = isequal(size(pp),[2,2]) && isequal(cell2mat(pp),[10,20;10,20]) && size(im{1},1)==2;
nfail = nfail + ~ok;
chk(ok,'N=2: kernels repeated along the first direction');
% All integrals over empty regions: the functional is zero
try
    intop.canonicalize_intmat({[3,4,1,3,0,4]},{1});     ok = false;
catch
    ok = true;
end
nfail = nfail + ~ok;
chk(ok,'error if all rows of a direction have an empty region');

%% Part 3: apply_op against the zero-convention oracle
nmis = 0;      ntrial = 0;      nsplit = 0;     nrej = 0;       nzero = 0;
for trial=1:60
    N = randi(3);
    d = randi([2,4],1,N);
    intmat = cell(1,N);     nrw = zeros(1,N);
    for i=1:N
        % Up to 3 random valid rows per direction, most of them not ordered
        rr = zeros(0,2*d(i));
        for att=1:300
            % each limit is a variable (60%), else a or b, and sometimes a multiplier
            r = zeros(1,2*d(i));
            for j=1:d(i)
                oth = setdiff(1:d(i),j);
                if rand<0.6,    r(2*j-1) = oth(randi(numel(oth)));  end
                r(2*j) = d(i)+1;
                if rand<0.6,    r(2*j) = oth(randi(numel(oth)));    end
                if rand<0.1 && r(2*j-1)>0,  r(2*j) = r(2*j-1);  end
            end
            if is_valid_row(r,d(i))
                rr(end+1,:) = r;                                            %#ok<SAGROW>
            end
            if size(rr,1)>=3,   break,  end
        end
        if isempty(rr),     rr = repmat([0,d(i)+1],1,d(i));     end
        intmat{i} = unique(rr,'rows');
        nrw(i) = size(intmat{i},1);
    end
    % Kernels per direction and term, state per direction
    vars = arrayfun(@(i) sprintf('s%d',i),1:N,'UniformOutput',false);
    pv = cell(1,N);     k_i = cell(1,N);    v_i = cell(1,N);
    for i=1:N
        pv{i} = polynomial(arrayfun(@(k) sprintf('s%d_dum%d',i,k),(1:d(i))','UniformOutput',false));
        k_i{i} = cell(nrw(i),1);
        for r=1:nrw(i)
            k_i{i}{r} = rand_poly([1,1],pv{i},2*ones(d(i),1),6);
        end
        v_i{i} = rand_poly([1,1],pv{i},ones(d(i),1),6);
    end
    params = cell([nrw,1]);
    for lin=1:numel(params)
        sub = cell(1,N);        [sub{:}] = ind2sub([nrw,1],lin);
        kk = polynomial(1);
        for i=1:N,  kk = kk*k_i{i}{sub{i}};  end
        params{lin} = kk;
    end
    % Oracle: product over directions of the sum over terms
    ref = 1;
    for i=1:N
        Ii = 0;
        for r=1:nrw(i)
            Ii = Ii + zero_integral(intmat{i}(r,:),k_i{i}{r}*v_i{i},pv{i});
        end
        ref = ref*Ii;
    end
    % The constructor rejects rows that are not ordered, pointing to canonicalize_intmat
    [~,~,io] = intop.canonicalize_intmat(intmat,[],true);
    is_unord = any(cellfun(@(a) ~all(a),io));
    try
        intop(params,intmat,vars,[0,1]);        rejected = false;
    catch ME
        rejected = contains(ME.message,'not ordered') && contains(ME.message,'canonicalize_intmat');
    end
    if is_unord
        nrej = nrej + 1;
        if ~rejected,   nmis = nmis+1;  fprintf('FAIL not rejected: %s\n',mat2str_cell(intmat));  end
    elseif rejected
        nmis = nmis+1;  fprintf('FAIL rejected, but ordered: %s\n',mat2str_cell(intmat));
    end
    % Declare the functional from canonicalize_intmat
    try
        [intmat_c,params_c] = intop.canonicalize_intmat(intmat,params);
    catch
        % A direction in which all regions are empty: the functional is zero
        nzero = nzero + 1;
        if abs(ref)>1e-12,  nmis = nmis+1;  fprintf('FAIL zero functional, oracle %.6g: %s\n',ref,mat2str_cell(intmat));  end
        continue
    end
    ntrial = ntrial+1;
    nsplit = nsplit + any(cellfun(@(a,b) ~isequal(a,b),intmat_c,intmat));
    Kop = intop(params_c,intmat_c,vars,[0,1]);
    v = polynomial(1);
    for i=1:N,  v = v*v_i{i};  end
    val = double(apply_op(Kop,v));
    if abs(val-ref) > 1e-9*max(1,abs(ref))
        nmis = nmis+1;
        fprintf('FAIL mismatch: N=%d intmat=%s: apply_op %.6g, oracle %.6g\n',N,mat2str_cell(intmat),val,ref);
    end
end
nfail = nfail + nmis;
chk(nmis==0,sprintf('random rows: %d declared (%d with a split, %d zero functionals), %d rejected as given, %d mismatches',ntrial,nsplit,nzero,nrej,nmis));
if nsplit<5
    nfail = nfail+1;
    disp('FAIL too few random trials exercised a split');
end

% Hand-built rows, one direction, with a fixed kernel and state
vars = {'s'};
for k=1:numel(rows_unord)
    r = rows_unord{k};      d = numel(r)/2;
    pvk = polynomial(arrayfun(@(j) sprintf('s_dum%d',j),(1:d)','UniformOutput',false));
    kern = 1 + pvk(1) + 2*pvk(d) + pvk(2)*pvk(d);
    vst = 1 + pvk(1)*pvk(2);
    ref = zero_integral(r,kern*vst,pvk);
    try
        [im,pp] = intop.canonicalize_intmat({r},{kern});
        val = double(apply_op(intop(pp,im,vars,[0,1]),vst));
    catch
        val = 0;            % empty region
    end
    ok = abs(val-ref)<1e-10;
    nfail = nfail + ~ok;
    chk(ok,sprintf('hand-built [%s]: apply_op %.6g, oracle %.6g',num2str(r),val,ref));
end
% Example of the documentation: int int int_{s1}^{s3} s1 ds2 ds3 ds1 = 1/24
% (indicator), where reversed limits counted negatively give -1/12
pvk = polynomial({'s_dum1';'s_dum2';'s_dum3'});
[im,pp] = intop.canonicalize_intmat({[0,4,1,3,0,4]},{pvk(1)});
Kop = intop(pp,im,vars,[0,1]);
ok = abs(double(apply_op(Kop,polynomial(1)))-1/24)<1e-12;
nfail = nfail + ~ok;
chk(ok,'example [0,4,1,3,0,4], kernel s1: 1/24');

%% Part 4: isvalid and the constructor
Kbad = Kop;     Kbad.intmat{1} = [0,4,1,3,0,4];
[lv,em] = isvalid(Kbad);
ok = ~lv && contains(em,'not ordered') && contains(em,'canonicalize_intmat');
nfail = nfail + ~ok;
chk(ok,'isvalid reports a manually assigned row that is not ordered, and refers to canonicalize_intmat');
Kok = Kop;      Kok.intmat{1} = [0,4,0,3,1,4];
ok = isvalid(Kok);
nfail = nfail + ~ok;
chk(ok,'isvalid accepts an ordered row');
try
    intop({pvk(1)},{[0,4,1,3,0,4]},{'s'},[0,1]);        ok = false;
catch ME
    ok = contains(ME.message,'canonicalize_intmat');
end
nfail = nfail + ~ok;
chk(ok,'constructor (chkval=false) rejects a row that is not ordered');
try
    Kskip = intop({pvk(1)},{[0,4,1,3,0,4]},{'s'},[0,1],true);
    ok = isequal(Kskip.intmat{1},[0,4,1,3,0,4]);
catch
    ok = false;
end
nfail = nfail + ~ok;
chk(ok,'constructor with chkval=true does not check the order');

%% Part 5: composition with tensor-PI operators
dd = 3;
vars = struct('in',{{'s'}},'out',{{'s'}});
dom = struct('in',[0,1],'out',[0,1]);
depmat = struct('in',true(dd,1),'out',true(dd,1));
nbad = 0;       nrun = 0;
for seed=1:6
    rng(seed);
    Bop = rand_tensopvar(ones(dd,2),vars,dom,depmat,1,[2,2],struct('in',1,'out',1));
    vars_A = vars;      vars_A.in = vars.out;   vars_A.out = cell(1,0);
    dom_A = dom;        dom_A.in = dom.out;     dom_A.out = zeros(0,2);
    pv = polynomial({'s_dum1';'s_dum2';'s_dum3'});
    par = rand_poly([1,1],pv,ones(dd,1),10);
    [im,pp] = intop.canonicalize_intmat({[0,4,1,3,0,4]},{par});
    Aop = intop(pp,im,vars_A,dom_A);
    pvars = polynomial(Bop.vars.in(:));
    v = cell(1,dd);     v_full = 1;
    for j=1:dd
        v{j} = rand_poly([size(Bop.ops{1,j},2),1],pvars,2,10);
        v_full = kron(v_full,subs(v{j},pvars,polynomial({['s_dum',num2str(j)]})));
    end
    Cop = mtimes_tensopvar(Aop,Bop);
    ABv = apply_op(Aop,apply_op(Bop,v));
    Cv = apply_op(Cop,v_full);
    err = cleanpoly(ABv-Cv,1e-10);
    e = full(max(max(abs(err.C))));
    nrun = nrun+1;
    nbad = nbad + (e > 1e-8*max(1,full(max(max(abs(double(ABv)))))));
end
nfail = nfail + nbad;
chk(nbad==0,sprintf('mtimes_tensopvar of a split functional with %d random 3-PI operators: %d mismatches',nrun,nbad));

if nfail>0
    error("Test_intop_ordered_intmat: %d check(s) failed.",nfail)
end
disp('Test_intop_ordered_intmat passed.')


%% Local functions
function out = ternary(c,a,b)
if c, out = a; else, out = b; end
end

function s = mat2str_cell(c)
s = strjoin(cellfun(@mat2str,c,'UniformOutput',false),' | ');
end

function ok = is_valid_row(r,d)
% Checks the validity of the row of an intmat (apart from the order), by
% assigning it to a functional and calling isvalid, which only reports the
% order after all other checks have passed
ok = false;
try
    Kt = intop({polynomial(1)},{repmat([0,d+1],1,d)},{'s'},[0,1],true);
    Kt.intmat{1} = r;
    [lv,em] = isvalid(Kt);
    ok = lv || contains(em,'not ordered');
catch
end
end

function I = zero_integral(row,P,pv)
% Integral over [0,1]^d of the polynomial P in the variables pv, with each
% variable j restricted to the interval defined by chi_{[x_lj,x_uj]}(x_j),
% taken to be a true indicator, by brute force over all orderings
d = numel(row)/2;
L = row(1:2:end);       U = row(2:2:end);
is_dlt = L==U;
tgt = 0:d+1;            tgt(find(is_dlt)+1) = L(is_dlt);
for k=1:d,  tgt = tgt(tgt+1);   end
% Substitute the multipliers
for j=find(is_dlt)
    t = tgt(j+1);
    if t==0,        val = polynomial(0);
    elseif t==d+1,  val = polynomial(1);
    else,           val = pv(t);
    end
    P = subs(P,pv(j),val);
end
fre = find(~is_dlt);        n = numel(fre);
Lr = tgt(L+1);              Ur = tgt(U+1);
I = 0;
if any(Lr(fre)==Ur(fre) & Lr(fre)>=1 & Lr(fre)<=d)
    return          % zero-width interval
end
perms_n = perms(1:n);
for r=1:size(perms_n,1)
    ord = fre(perms_n(r,:));            % x_ord(1) <= ... <= x_ord(n)
    pos = zeros(1,d);   pos(ord) = 1:n;
    good = true;
    for j=fre
        if Lr(j)>=1 && Lr(j)<=d,  good = good && pos(Lr(j))<=pos(j);  end
        if Ur(j)>=1 && Ur(j)<=d,  good = good && pos(j)<=pos(Ur(j));  end
    end
    if ~good,   continue,   end
    % Integrate over the cell a<=x_ord(1)<=...<=x_ord(n)<=b
    Q = P;
    for k=1:n
        up = polynomial(1);
        if k<n,     up = pv(ord(k+1));  end
        Q = int(Q,pv(ord(k)).varname{1},0,up);
    end
    I = I + double(Q);
end
end
