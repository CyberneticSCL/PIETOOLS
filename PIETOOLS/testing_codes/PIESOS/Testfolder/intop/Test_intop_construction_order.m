% Test_intop_construction_order.m -- Tests that the 'intop' constructor
% sorts vars.in alphabetically regardless of the order in which the caller
% supplied it, correctly permuting dom.in/intmat/params to match (so the
% represented functional is unchanged), and that this makes two
% functionals over the same variable set, built in different given
% orders, combine via '+'. Also covers the N=0 and N=1 edge cases.
%
% The shuffled functional is declared twice: from the kernels,
% kernels(Aop), and from the coefficients in the basis of the shuffled
% variable order (columns relabeled in the test from the definition of the
% Kronecker basis).
%
% DJ (CC), 10/15/2026: Initial coding;

tol = 1e-8;

for trial = 1:10
    rng(6900+trial);

    N = randi([2,4]);
    dom = randi([-2,2],[N,1]) + [zeros(N,1),2*rand([N,1])];
    degs = randi([1,2],1,N);
    pdegs = randi([0,2],1,N);
    Aop = rand_intop([randi(2),randi(2)],N,dom,degs,pdegs);

    % rand_intop builds vars {'s1',...,'sN'}, already sorted for N<=4
    if ~isequal(Aop.vars.in,sort(Aop.vars.in))
        error("Test_intop_construction_order failed: vars.in not sorted after construction (trial=%d).",trial)
    end

    % Same functional with the variables (and dom.in/intmat/params) given in a
    % shuffled order: must land on the same representation as Aop
    perm = randperm(N);
    while isequal(perm,1:N)
        perm = randperm(N);
    end
    vars_shuffled = Aop.vars.in(perm);
    dom_shuffled = Aop.dom.in(perm,:);
    intmat_shuffled = Aop.intmat(perm);
    Akernels = kernels(Aop);
    params_shuffled = permute(reshape(Akernels,[size(Aop.params,1:N),1]),[perm,N+1]);

    Bop = intop(params_shuffled,intmat_shuffled,vars_shuffled,dom_shuffled);

    % From the coefficients in the basis of the shuffled variables: the monomial with
    % subscripts z (order of Aop.vars.in) is number sum((z(i)-1)*stride(i))+1 in the
    % Kronecker basis (first variable slowest)
    nm = cellfun(@(a) size(a,1),Aop.Z.in);
    nZ = prod(nm);
    sub = cell(1,N);
    [sub{:}] = ind2sub(fliplr(nm),(1:nZ)');
    sub = fliplr(sub);                  % subscripts of each variable, first variable first
    stride_new = [fliplr(cumprod(fliplr(nm(perm(2:end))))),1];
    znew = ones(nZ,1);
    for ip=1:N
        znew = znew + (sub{perm(ip)}-1)*stride_new(ip);
    end
    nd = numel(Aop.Zd);
    coeffs_shuffled = cell(size(Aop.params));
    for k=1:numel(Aop.params)
        [ii,jj,vv] = find(Aop.params{k});
        ii = ii(:);     jj = jj(:);     vv = vv(:);     % row vectors if there is one row
        jb = ceil(jj/nZ);   zz = jj-(jb-1)*nZ;
        coeffs_shuffled{k} = sparse(ii,(jb-1)*nZ+znew(zz),vv,size(Aop.params{k},1),size(Aop.params{k},2));
    end
    coeffs_shuffled = permute(reshape(coeffs_shuffled,[size(Aop.params,1:N),1]),[perm,N+1]);
    Bop_c = intop(coeffs_shuffled,intmat_shuffled,vars_shuffled,dom_shuffled,Aop.Z.in(perm),Aop.Zd);

    % Bop.params is not compared with isequal: isequal of matrix-valued polynomials is
    % elementwise and does not compose inside a cell-array isequal. The apply_op check
    % below is used instead.
    if ~isequal(Bop.vars.in,Aop.vars.in)
        error("Test_intop_construction_order failed: shuffled-order construction did not re-sort vars.in to match (trial=%d).",trial)
    elseif ~isequal(Bop.dom.in,Aop.dom.in)
        error("Test_intop_construction_order failed: dom.in not correctly permuted back (trial=%d).",trial)
    elseif ~isequal(Bop.intmat,Aop.intmat)
        error("Test_intop_construction_order failed: intmat not correctly permuted back (trial=%d).",trial)
    end

    % Functional check: both give the same value on the same input
    pvars_in = Aop.pvars.in;
    pvars_full = polynomial(zeros(0,1));
    for i=1:N
        pvars_full = [pvars_full; pvars_in{i}(:)];
    end
    v = rand_poly([size(Aop,2),1],pvars_full,2,20);
    valA = apply_op(Aop,v);
    valB = apply_op(Bop,v);
    err = cleanpoly(valA-valB,1e-10);
    if max(max(abs(full(err.C))))>tol
        error("Test_intop_construction_order failed: shuffled-order construction's functional does not match the original (trial=%d).",trial)
    end
    % Same for the declaration from the coefficients
    if ~isequal(Bop_c.vars.in,Aop.vars.in) || ~isequal(Bop_c.dom.in,Aop.dom.in) || ~isequal(Bop_c.intmat,Aop.intmat)
        error("Test_intop_construction_order failed: shuffled-order construction from coefficients did not restore the order (trial=%d).",trial)
    end
    valBc = apply_op(Bop_c,v);
    err = cleanpoly(valA-valBc,1e-10);
    if max(max(abs(full(err.C))))>tol
        error("Test_intop_construction_order failed: shuffled-order construction from coefficients does not match the original (trial=%d).",trial)
    end

    % Aop and Bop are internally identical, so Aop+Bop is 2*Aop
    Cop = Aop + Bop;
    valC = apply_op(Cop,v);
    err2 = cleanpoly(valC-2*valA,1e-10);
    if max(max(abs(full(err2.C))))>tol
        error("Test_intop_construction_order failed: Aop+Bop (built with different variable orders) did not give 2*Aop (trial=%d).",trial)
    end
end

%%%% Edge cases: N=0 and N=1 are unaffected (sorting 0 or 1 variables is a no-op)
Kop0 = intop(5);
if numel(Kop0.vars.in)~=0
    error("Test_intop_construction_order failed: N=0 case should have empty vars.in.")
end
Kop1 = rand_intop([1,1],1,[0,1],1,1);
if numel(Kop1.vars.in)~=1
    error("Test_intop_construction_order failed: N=1 case should have exactly 1 variable.")
end

disp('Test_intop_construction_order passed.')
