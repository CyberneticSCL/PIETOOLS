% Test_tensopmat_common_vars.m -- tests @tensopmat/common_vars.m directly
% at the field level. Oracle: the rebased vars/dom/depmat are checked
% against the ORIGINAL (pre-merge) operators by name lookup, independent of
% unique()'s sort order; catches a wrong-axis merge or a mis-mapped depmat
% column that an apply_op-based test could miss.
% A and B share their INPUT-side pool by name (partial overlap) but use
% disjoint OUTPUT-side pools ('aOut' vs 'bOut'), exercising both the
% overlapping and fully disjoint cases.

% DJ (CC), 10/15/2026: Initial coding;

pA = randi([1,3]);  dA = randi([1,3]);  NA = randi([2,4]);  MA = randi([2,3]);
pB = randi([1,3]);  dB = randi([1,3]);  NB = randi([2,4]);  MB = randi([2,3]);
dom = [0,1];

[Aop,~] = rand_tensopmat(pA,dA,NA,MA,dom,0.7,[],struct(),'aIn','aOut');
[Bop,~] = rand_tensopmat(pB,dB,NB,MB,dom,0.7,[],struct(),'aIn','bOut');

[A2,B2] = common_vars(Aop,Bop);

% Both operators must end up on the exact same basis
if ~isequal(A2.vars.in,B2.vars.in) || ~isequal(A2.dom.in,B2.dom.in)
    error("Test_tensopmat_common_vars failed: input vars/dom not unified.")
end
if ~isequal(A2.vars.out,B2.vars.out) || ~isequal(A2.dom.out,B2.dom.out)
    error("Test_tensopmat_common_vars failed: output vars/dom not unified.")
end
% New basis must be the union of the original var lists
if ~isequal(sort(A2.vars.in),sort(unique([Aop.vars.in,Bop.vars.in])))
    error("Test_tensopmat_common_vars failed: input basis is not the union of the originals.")
end
if ~isequal(sort(A2.vars.out),sort(unique([Aop.vars.out,Bop.vars.out])))
    error("Test_tensopmat_common_vars failed: output basis is not the union of the originals.")
end

% depmat column-by-column by name: zero for a variable the operator never
% depended on, the original value otherwise
check_depmat_side(A2.depmat.in, A2.vars.in, Aop.depmat.in, Aop.vars.in, 'A.in');
check_depmat_side(B2.depmat.in, B2.vars.in, Bop.depmat.in, Bop.vars.in, 'B.in');
check_depmat_side(A2.depmat.out, A2.vars.out, Aop.depmat.out, Aop.vars.out, 'A.out');
check_depmat_side(B2.depmat.out, B2.vars.out, Bop.depmat.out, Bop.vars.out, 'B.out');

% Block counts must be unchanged; only the variable axis grows
if size(A2.depmat.in,1)~=size(Aop.depmat.in,1) || size(B2.depmat.in,1)~=size(Bop.depmat.in,1)
    error("Test_tensopmat_common_vars failed: number of column-blocks changed by common_vars.")
end
if size(A2.depmat.out,1)~=size(Aop.depmat.out,1) || size(B2.depmat.out,1)~=size(Bop.depmat.out,1)
    error("Test_tensopmat_common_vars failed: number of row-blocks changed by common_vars.")
end

disp('Test_tensopmat_common_vars passed.')

function check_depmat_side(depmat_new,vars_new,depmat_old,vars_old,desc)
for k=1:numel(vars_new)
    v = vars_new{k};
    idx_old = find(strcmp(vars_old,v));
    if isempty(idx_old)
        if any(depmat_new(:,k)~=0)
            error("Test_tensopmat_common_vars failed: %s column for unseen variable '%s' is nonzero.",desc,v)
        end
    else
        if ~isequal(depmat_new(:,k),depmat_old(:,idx_old))
            error("Test_tensopmat_common_vars failed: %s column for variable '%s' does not match original.",desc,v)
        end
    end
end
end
