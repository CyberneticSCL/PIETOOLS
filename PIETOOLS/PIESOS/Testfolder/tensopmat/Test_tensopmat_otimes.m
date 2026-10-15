% Test_tensopmat_otimes.m -- tests @tensopmat/otimes.m.
%
% Structural checks (sizes, depmat row counts), plus two further checks:
% (1) A passed as a raw 'tensopvar' with B a genuine 'tensopmat' must
%     convert/check B;
% (2) the type(2)==1 column-dimension check errors on a mismatch, not on a
%     match.
% Numeric ground truth for type=[2,2] and [1,2] (otype(2) is fixed at 2,
% so otype(1) in {1,2} is the only axis of 'type' reachable). Oracle,
% built from the operands alone: for [2,2], C*Zfull = kron(A*ZA,B*ZB)
% row-block by row-block; for [1,2], C*Zfull = (A*ZA).*(B*ZB). Both
% distribute bilinearly over each operand's term/column sum. 20 fresh
% random reps.

% DJ (CC), 10/15/2026: Initial coding;

% NOTES
% - Each combined block concatenates into a 2-factor 'tensopvar' for
%   either 'type', so apply_op on the otimes result takes a NESTED 1x2
%   cell input per new column (A-sourced factor, B-sourced factor), unlike
%   other tests in this folder;
% - @tensopvar/otimes.m forces A.otype=B.otype=type for single-factor
%   operands, so the otype(1) sweep of the inputs is nearly a no-op here;
%   nterms/degs are what vary;
% - The combined block's otype(1) is forced to 'type', so the "_<k>"
%   renaming of the output variables by apply_op occurs whenever
%   type(1)==2; both sides are passed through strip_kron_suffix (see
%   Test_tensopmat_mtimes.m);
% - The numeric checks use disjoint pools and blk_dnsty=1 (every block
%   populated) to keep the expected-value bookkeeping simple; [1,1] and
%   [2,1] are only checked structurally, [1,1] only if the retry finds
%   matching dims.

tol = 1e-8;
nreps = 20;

for rep_idx=1:nreps

p = randi([1,3]);  d = randi([1,3]);
N = randi([1,3]);  M = randi([1,3]);
dom = [0,1];
if rand()<0.3
    dnsty = 0.3+0.4*rand();
else
    dnsty = [];
end

[Aop,~] = rand_tensopmat(p,d,N,M,dom,0.7,dnsty,struct(),'a','b','rand',3,3);
[Bop,~] = rand_tensopmat(p,d,N,M,dom,0.7,dnsty,struct(),'c','e','rand',3,3);

% Full Kronecker product (type=[2,2]): shape only, exact depmat values would
% duplicate the common_vars bookkeeping (numeric check below is independent)
C22 = otimes(Aop,Bop,[2,2]);
if ~isequal(size(C22),[size(Aop,1)*size(Bop,1),size(Aop,2)*size(Bop,2)])
    error("Test_tensopmat_otimes failed: size mismatch for type=[2,2].")
end
if ~isequal(size(C22.depmat.in,1),size(Aop.ops,2)*size(Bop.ops,2)) ...
        || ~isequal(size(C22.depmat.out,1),size(Aop.ops,1)*size(Bop.ops,1))
    error("Test_tensopmat_otimes failed: depmat row-count mismatch for type=[2,2].")
end

% Pointwise along both dimensions (type=[1,1]) needs matching row AND column
% dims: same block counts as Aop, retried (max 200) until the random
% per-block dims match
p11 = numel(Aop.dims.out);   d11 = numel(Aop.dims.in);
[Bop11,~] = rand_tensopmat(p11,d11,N,M,dom,0.7,dnsty,struct(),'x','y','rand',3,3);
ntries = 0;
while (~isequal(Bop11.dims.out,Aop.dims.out) || ~isequal(Bop11.dims.in,Aop.dims.in)) && ntries<200
    [Bop11,~] = rand_tensopmat(p11,d11,N,M,dom,0.7,dnsty,struct(),'x','y','rand',3,3);
    ntries = ntries+1;
end
if isequal(Bop11.dims.out,Aop.dims.out) && isequal(Bop11.dims.in,Aop.dims.in)
    C11 = otimes(Aop,Bop11,[1,1]);
    if ~isequal(size(C11),size(Aop))
        error("Test_tensopmat_otimes failed: size mismatch for type=[1,1].")
    end
end

% Regression (1): A a raw 'tensopvar', B a tensopmat
Araw = Aop.ops{find(~cellfun(@isempty,Aop.ops),1)};
Craw = otimes(Araw,Bop,[2,2]);
if ~isequal(size(Craw),[size(Araw,1)*size(Bop,1),size(Araw,2)*size(Bop,2)])
    error("Test_tensopmat_otimes failed: raw-tensopvar-times-tensopmat size mismatch (regression 1).")
end

% Regression (2): type(2)==1 column-dimension check; the matching case is
% the Bop11 product above, the mismatched case must error
if size(Aop,2)~=size(Bop,2)
    did_error = false;
    try
        otimes(Aop,Bop,[2,1]); %#ok<NASGU>
    catch
        did_error = true;
    end
    if ~did_error
        error("Test_tensopmat_otimes failed: mismatched column dims with type(2)==1 should have errored (regression 2).")
    end
end

% --- Numeric ground truth, type=[2,2] ---------
% disjoint pools ('ak'/'bk' vs 'ck'/'ek'); C*Zfull = kron(A*ZA,B*ZB)
% row-block by row-block, I=(ridx-1)*pB+ib
pAk = randi([1,3]);  dAk = randi([1,3]);  NAk = randi([1,3]);  MAk = randi([1,3]);
pBk = randi([1,3]);  dBk = randi([1,3]);  NBk = randi([1,3]);  MBk = randi([1,3]);
[Aop_k,infoAk] = rand_tensopmat(pAk,dAk,NAk,MAk,dom,1,dnsty,struct(),'ak','bk','rand',3,3);
[Bop_k,infoBk] = rand_tensopmat(pBk,dBk,NBk,MBk,dom,1,dnsty,struct(),'ck','ek','rand',3,3);

XAk = cell(1,dAk);
for j=1:dAk
    XAk{j} = rand_poly([infoAk.dims_in(j),1],infoAk.vars_in_pool(infoAk.sel_in{j}),1,3);
end
XBk = cell(1,dBk);
for j=1:dBk
    XBk{j} = rand_poly([infoBk.dims_in(j),1],infoBk.vars_in_pool(infoBk.sel_in{j}),1,3);
end
YAk = apply_op(Aop_k,XAk);
YBk = apply_op(Bop_k,XBk);

Ck = otimes(Aop_k,Bop_k,[2,2]);
Xfull_k = cell(1,dAk*dBk);
for cidx=1:dAk
    for jb=1:dBk
        Xfull_k{(cidx-1)*dBk+jb} = {XAk{cidx},XBk{jb}};
    end
end
Yk =apply_op(Ck,Xfull_k);

nnAk = [0;cumsum(Aop_k.dims.out)];
nnBk = [0;cumsum(Bop_k.dims.out)];
Yk_exp_blocks = cell(pAk*pBk,1);
for ridx=1:pAk
    YAridx = YAk(nnAk(ridx)+1:nnAk(ridx+1),:);
    for ib=1:pBk
        YBib = YBk(nnBk(ib)+1:nnBk(ib+1),:);
        Yk_exp_blocks{(ridx-1)*pBk+ib} = kron(YAridx,YBib);
    end
end
Yk_exp = vertcat(Yk_exp_blocks{:});

err_k = cleanpoly(strip_kron_suffix(Yk)-strip_kron_suffix(Yk_exp),1e-10);
if max(max(abs(full(err_k.C))))>tol
    error("Test_tensopmat_otimes failed: numeric mismatch for type=[2,2] (kron ground truth).")
end

% --- Numeric ground truth, type=[1,2] ---------
% type(1)==1 needs the same row-blocks and per-row dims.out: Bop_p's output
% side is taken from Aop_p via info_in. C*Zfull = (A*ZA).*(B*ZB)
pAp = randi([1,3]);  dAp = randi([1,3]);  NAp = randi([1,3]);  MAp = randi([1,3]);
dBp = randi([1,3]);  NBp = randi([1,3]);
[Aop_p,infoAp] = rand_tensopmat(pAp,dAp,NAp,MAp,dom,1,dnsty,struct(),'ap','bp','rand',3,3);
infoBp_out = struct('vars_out_pool',{infoAp.vars_out_pool},'dom_out_pool',infoAp.dom_out_pool,...
                     'sel_out',{infoAp.sel_out},'dims_out',infoAp.dims_out);
[Bop_p,infoBp] = rand_tensopmat([],dBp,NBp,[],dom,1,dnsty,infoBp_out,'cp','ep','rand',3,3);

XAp = cell(1,dAp);
for j=1:dAp
    XAp{j} = rand_poly([infoAp.dims_in(j),1],infoAp.vars_in_pool(infoAp.sel_in{j}),1,3);
end
XBp = cell(1,dBp);
for j=1:dBp
    XBp{j} = rand_poly([infoBp.dims_in(j),1],infoBp.vars_in_pool(infoBp.sel_in{j}),1,3);
end
YAp = apply_op(Aop_p,XAp);
YBp = apply_op(Bop_p,XBp);

Cp = otimes(Aop_p,Bop_p,[1,2]);
Xfull_p = cell(1,dAp*dBp);
for cidx=1:dAp
    for jb=1:dBp
        Xfull_p{(cidx-1)*dBp+jb} = {XAp{cidx},XBp{jb}};
    end
end
Yp =apply_op(Cp,Xfull_p);
Yp_exp = YAp.*YBp;

err_p = cleanpoly(strip_kron_suffix(Yp)-strip_kron_suffix(Yp_exp),1e-10);
if max(max(abs(full(err_p.C))))>tol
    error("Test_tensopmat_otimes failed: numeric mismatch for type=[1,2] (pointwise ground truth).")
end

end

disp('Test_tensopmat_otimes passed.')

function Pout = strip_kron_suffix(P)
% POUT = STRIP_KRON_SUFFIX(P) strips a trailing "_<int>" suffix
% unconditionally; see Test_tensopmat_mtimes.m (same helper, duplicated so
% each test file stays self-contained).
Pout = P;
vn = Pout.varname;
for k=1:numel(vn)
    tok = regexp(vn{k},'^(.*)_\d+$','tokens','once');
    if ~isempty(tok)
        vn{k} = tok{1};
    end
end
Pout.varname = vn;
end
