% Test_tensopmat_clean_op_dead_vars.m -- tests that @tensopmat/clean_op.m
% drops an input/output variable that no block structurally depends on (an
% entirely-zero column of depmat.in/depmat.out), leaving every used
% variable and the operator's action untouched; and that a 'tensopmat'
% with no dead variables is left unchanged. Swept over varying numbers of
% row/column blocks and spatial variables.
% Oracle: the original vars/dom/depmat and apply_op of the original.

% DJ (CC), 10/15/2026: Initial coding;

% NOTES
% - A fresh rand_tensopmat object never has a dead variable (initialize.m
%   seeds the pools so the contributing row/column has nonzero depmat), so
%   one is injected directly.

tol = 1e-8;

for trial = 1:10
    rng(6700+trial);

    p = randi(3)+1;   d = randi(3)+1;
    N = randi(3)+1;   M = randi(3)+1;

    Pop = rand_tensopmat(p,d,N,M,[0,1],0.7,0.8);

    % Input covering every column
    X = cell(1,d);
    for j=1:d
        i0 = find(~cellfun(@isempty,Pop.ops(:,j)),1);
        vars_ij = Pop.ops{i0,j}.vars.in;
        X{j} = rand_poly([Pop.dims.in(j),1],polynomial(vars_ij(:)),2,20);
    end
    Y0 = apply_op(Pop,X);

    %%%% Case 1: no dead variables -- vars/dom/depmat unchanged.
    [Pop_clean0,~] = clean_op(Pop);
    if ~isequal(Pop_clean0.vars,Pop.vars) || ~isequal(Pop_clean0.dom,Pop.dom) || ...
            ~isequal(Pop_clean0.depmat,Pop.depmat)
        error("Test_tensopmat_clean_op_dead_vars failed: vars/dom/depmat changed with no dead variables present (trial=%d).",trial)
    end

    %%%% Case 2: inject one dead output and one dead input variable (zero
    %%%% depmat column) -- clean_op should drop exactly these.
    Pop_dead = Pop;
    Pop_dead.vars.out = [Pop_dead.vars.out, {'dead_out'}];
    Pop_dead.dom.out = [Pop_dead.dom.out; [0,1]];
    Pop_dead.depmat.out = [Pop_dead.depmat.out, zeros(p,1)];
    Pop_dead.vars.in = [Pop_dead.vars.in, {'dead_in'}];
    Pop_dead.dom.in = [Pop_dead.dom.in; [0,1]];
    Pop_dead.depmat.in = [Pop_dead.depmat.in, zeros(d,1)];

    [Pop_clean,is_zero] = clean_op(Pop_dead);
    if is_zero
        error("Test_tensopmat_clean_op_dead_vars failed: operator incorrectly flagged entirely zero (trial=%d).",trial)
    end
    if any(strcmp(Pop_clean.vars.out,'dead_out')) || any(strcmp(Pop_clean.vars.in,'dead_in'))
        error("Test_tensopmat_clean_op_dead_vars failed: dead variable not removed (trial=%d).",trial)
    end
    if ~isequal(Pop_clean.vars,Pop.vars) || ~isequal(Pop_clean.dom,Pop.dom) || ~isequal(Pop_clean.depmat,Pop.depmat)
        error("Test_tensopmat_clean_op_dead_vars failed: cleaning the dead-variable object did not reduce exactly to the original vars/dom/depmat (trial=%d).",trial)
    end

    % Action must be unaffected (apply_op never reads vars/dom/depmat, so
    % this holds trivially)
    Y_dead = apply_op(Pop_dead,X);
    Y_clean = apply_op(Pop_clean,X);
    err1 = cleanpoly(Y_dead-Y0,1e-10);
    err2 = cleanpoly(Y_clean-Y0,1e-10);
    if max(max(abs(full(err1.C))))>tol || max(max(abs(full(err2.C))))>tol
        error("Test_tensopmat_clean_op_dead_vars failed: operator's action changed (trial=%d).",trial)
    end
end

disp('Test_tensopmat_clean_op_dead_vars passed.')
