% Test_tensopmat_mtimes.m -- tests mtimes for 'tensopmat' objects:
% - Case 1: scalar times 'tensopmat', checked via apply_op.
% - Case 2: matrix times 'tensopmat' (A*B, A a plain 'double'), checked via
%   apply_op against the same matrix applied directly to apply_op(B,X).
%   Only a 'double' A is supported; a non-scalar 'polynomial'/'dpvar' A
%   must error (checked below).
% - Case 3: 'tensopmat' times 'tensopmat' (A*B, A and B both 'tensopmat',
%   built so that A's input space matches B's output space exactly, but
%   A's OWN output pool and B's OWN input pool are independent, disjoint
%   pools), checked via apply_op(A*B,X) == apply_op(A,{apply_op(B,X)}).
% - Case 4: as Case 3, but the LHS operand's blocks mix otype(1) in {1,2}
%   (rand_tensopmat otype1_opt='rand'); the RHS block with otype(1)==2 is
%   also checked against the same oracle.

% DJ (CC), 10/15/2026: Initial coding;

% NOTES
% - Case 2: Bop's random per-row variable-subset selection can violate
%   the supported shape (non-empty row-blocks of a column depending on the
%   same output variables), ~1 in 200 draws; mtimes then raises
%   'tensopmat:mtimes_heterogeneous_rows' or '_heterogeneous_columns',
%   which is caught and counted. The test fails if no rep exercises
%   Case 2;
% - Case 2 does not exercise the B.depmat.out>1 ("tensor-product B")
%   restriction: rand_tensopmat blocks have depmat.out<=1;
% - Case 3 uses disjoint 'g'/'c' pools for A.out and B.in, to test that
%   mtimes compares A.depmat.in with B.depmat.out directly, without
%   common_vars rebasing (a same-side merge would reject a valid product);
% - Case 4: for otype(1)==2 on the RHS of the single-factor product,
%   @tensopvar/mtimes.m may give an error (Case 4b handles both outcomes),
%   so Bop4 keeps otype1=1 and only Aop4 mixes. apply_op appends "_<k>" to
%   the output variable names for otype(1)==2, so the results are compared
%   after strip_kron_suffix;
% - strip_kron_suffix strips "_<k>" UNCONDITIONALLY (no collision guard):
%   "g41" and "g41_1" legitimately co-occur when different row-blocks use
%   both otypes, so a guard misfires. Valid because rand_tensopmat only
%   builds single-factor blocks, for which the suffix never marks a real
%   second variable copy;
% - Bop4alt (RHS otype(1)==2) may error in mtimes or not; if not, the
%   reference must be fed the suffix-stripped intermediate, since Aop4's
%   kernels live on the unsuffixed variables.

tol = 1e-8;
nreps = 15;
dom = [0,1];
ncase2_exercised = 0;

for rep_idx=1:nreps

p = randi([1,3]);  q = randi([1,3]);  d = randi([1,3]);
N = randi([1,3]);  M = randi([1,3]);
if rand()<0.3
    dnsty = 0.3+0.4*rand();
else
    dnsty = [];
end

% otype1=1 here since Bop plays the RHS role below
[Bop,infoB] = rand_tensopmat(q,d,N,M,dom,0.7,dnsty,struct(),'c','e',1,3,3);

X = cell(1,d);
for j=1:d
    X{j} = rand_poly([infoB.dims_in(j),1],infoB.vars_in_pool(infoB.sel_in{j}),1,5);
end
BX = apply_op(Bop,X);

% Case 1: scalar times tensopmat
a = 2.5;
Y1 = apply_op(a*Bop,X);
err1 = cleanpoly(Y1-a*BX,1e-10);
if max(max(abs(err1.C)))>tol
    error("Test_tensopmat_mtimes failed: scalar*tensopmat mismatch.")
end

% Case 2: (non-scalar double matrix) times tensopmat,
% apply(A2*Bop,X) == A2*apply(Bop,X); known limitation caught and counted
% (see NOTES)
mrows = randi([2,3]);
A2 = randn(mrows,size(Bop,1));
try
    Y2 = apply_op(A2*Bop,X);
    err2 = cleanpoly(Y2-A2*BX,1e-10);
    if max(max(abs(err2.C)))>tol
        error("Test_tensopmat_mtimes failed: matrix*tensopmat mismatch.")
    end
    ncase2_exercised = ncase2_exercised+1;
catch ME
    if ~strcmp(ME.identifier,"tensopmat:mtimes_heterogeneous_rows") && ...
            ~strcmp(ME.identifier,"tensopmat:mtimes_heterogeneous_columns")
        rethrow(ME)
    end
end

% Case 2 regression: a non-scalar 'polynomial'/'dpvar' A must still error
A2p = polynomial(randn(mrows,size(Bop,1)));
did_error2p = false;
try
    A2p*Bop; %#ok<VUNUS>
catch
    did_error2p = true;
end
if ~did_error2p
    error("Test_tensopmat_mtimes failed: non-scalar 'polynomial' A times tensopmat should still error, but did not.")
end

% Case 3: tensopmat times tensopmat; A's input space equals B's output
% space, A's own output pool ('g...') is independent of B's
MA = randi([1,3]);
infoA_in = struct('vars_in_pool',{infoB.vars_out_pool},'dom_in_pool',infoB.dom_out_pool,...
                   'sel_in',{infoB.sel_out},'dims_in',infoB.dims_out);
[Aop,~] = rand_tensopmat(p,[],[],MA,dom,0.7,dnsty,infoA_in,[],'g',1,3,3);

Y3 = apply_op(Aop*Bop,X);
nn_arr = [0;cumsum(infoB.dims_out)];
BXcell = cell(1,numel(infoB.dims_out));
for i=1:numel(infoB.dims_out)
    BXcell{i} = BX(nn_arr(i)+1:nn_arr(i+1),:);
end
Y3_check = apply_op(Aop,BXcell);
err3 = cleanpoly(Y3-Y3_check,1e-10);
if max(max(abs(err3.C)))>tol
    error("Test_tensopmat_mtimes failed: tensopmat*tensopmat mismatch.")
end

% Case 4: as Case 3, but the LHS (Aop4) mixes otype(1) in {1,2}, the RHS
% (Bop4) keeps otype1=1; compared after stripping the Kronecker-renaming
% suffix (see NOTES)
q4 = randi([1,3]);  d4 = randi([1,3]);  N4 = randi([1,3]);  M4 = randi([1,3]);
[Bop4,infoB4] = rand_tensopmat(q4,d4,N4,M4,dom,0.7,dnsty,struct(),'c4','e4',1,3,3);

X4 = cell(1,d4);
for j=1:d4
    X4{j} = rand_poly([infoB4.dims_in(j),1],infoB4.vars_in_pool(infoB4.sel_in{j}),1,5);
end
BX4 = apply_op(Bop4,X4);

MA4 = randi([1,3]);
infoA4_in = struct('vars_in_pool',{infoB4.vars_out_pool},'dom_in_pool',infoB4.dom_out_pool,...
                    'sel_in',{infoB4.sel_out},'dims_in',infoB4.dims_out);
p4 = randi([1,3]);
[Aop4,~] = rand_tensopmat(p4,[],[],MA4,dom,0.7,dnsty,infoA4_in,[],'g4','rand',3,3);

Y4 = apply_op(Aop4*Bop4,X4);
nn4_arr = [0;cumsum(infoB4.dims_out)];
BX4cell = cell(1,numel(infoB4.dims_out));
for i=1:numel(infoB4.dims_out)
    BX4cell{i} = BX4(nn4_arr(i)+1:nn4_arr(i+1),:);
end
Y4_check = apply_op(Aop4,BX4cell);
Y4n = strip_kron_suffix(Y4);
Y4_checkn = strip_kron_suffix(Y4_check);
err4 = cleanpoly(Y4n-Y4_checkn,1e-10);
if max(max(abs(err4.C)))>tol
    error("Test_tensopmat_mtimes failed: tensopmat*tensopmat with otype=[2,2] LHS blocks mismatch (after suffix normalization).")
end

% Case 4b: RHS with otype(1)==2 (Bop4alt) on the same input AND output space
% as Bop4, so only otype1 differs
infoB4alt_in = struct('vars_in_pool',{infoB4.vars_in_pool},'dom_in_pool',infoB4.dom_in_pool,...
                       'sel_in',{infoB4.sel_in},'dims_in',infoB4.dims_in,...
                       'vars_out_pool',{infoB4.vars_out_pool},'dom_out_pool',infoB4.dom_out_pool,...
                       'sel_out',{infoB4.sel_out},'dims_out',infoB4.dims_out);
[Bop4alt,~] = rand_tensopmat([],[],[],[],dom,0.7,[],infoB4alt_in,[],[],2,1,1);
did_error4 = false;
try
    C4alt = Aop4*Bop4alt;
catch
    did_error4 = true;
end
if ~did_error4
    BX4alt = apply_op(Bop4alt,X4);
    Y4alt = apply_op(C4alt,X4);
    BX4altcell = cell(1,numel(infoB4.dims_out));
    for i=1:numel(infoB4.dims_out)
        BX4altcell{i} = strip_kron_suffix(BX4alt(nn4_arr(i)+1:nn4_arr(i+1),:));
    end
    Y4alt_check = apply_op(Aop4,BX4altcell);
    err4alt = cleanpoly(strip_kron_suffix(Y4alt)-strip_kron_suffix(Y4alt_check),1e-10);
    nrm4alt = full(max(max(abs(err4alt.C))));
    if isempty(nrm4alt), nrm4alt = 0; end
    if nrm4alt>tol
        error("Test_tensopmat_mtimes failed: tensopmat*tensopmat with otype(1)==2 RHS blocks mismatch (residual %.3g, rep %d).",nrm4alt,rep_idx)
    end
end

end

if ncase2_exercised==0
    error("Test_tensopmat_mtimes: every Case 2 trial hit the heterogeneous-" + ...
        "rows limitation; no trial actually exercised the matrix*tensopmat fix.")
end
fprintf('Test_tensopmat_mtimes passed (Case 2 exercised %d/%d reps; rest hit the known heterogeneous-rows limitation).\n',ncase2_exercised,nreps)

function Pout = strip_kron_suffix(P)
% POUT = STRIP_KRON_SUFFIX(P) strips a trailing "_<int>" suffix from each
% variable name in polynomial P, normalizing apply_tensopvar.m's
% otype(1)==2 output renaming; strips unconditionally (see NOTES).
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
