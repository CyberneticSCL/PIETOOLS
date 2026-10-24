function Test_tensopmat_clean_op()
% Test_tensopmat_clean_op.m -- one local function for each case:
%   case_basic              that @tensopmat/clean_op.m cleans every
%                           (nonempty) block by dispatching to
%                           @tensopvar/clean_op.m, correctly determines
%                           whether the WHOLE block matrix is zero, and
%                           never
%   case_dead_vars          that @tensopmat/clean_op.m drops an
%                           input/output variable that no block
%                           structurally depends on (an entirely-zero
%                           column of depmat.in/depmat.out), leaving every
%                           used
%
% DJ (CC), 10/15/2026: Initial coding;
% DJ (CC), 10/23/2026: Merge Test_tensopmat_clean_op_dead_vars into this
%   test, as cases;

cases = {@case_basic, @case_dead_vars};
nfail = 0;
for i=1:numel(cases)
    tstart = tic;
    try
        cases{i}();
        fprintf('PASS  %s (%.1f s)\n',func2str(cases{i}),toc(tstart));
    catch err
        nfail = nfail+1;
        fprintf('FAIL  %s\n      %s\n',func2str(cases{i}),err.message);
        for k=1:min(3,numel(err.stack))
            fprintf('        in %s, line %d\n',err.stack(k).name,err.stack(k).line);
        end
    end
end
if nfail>0
    error('Test_tensopmat_clean_op: %d of %d cases failed.',nfail,numel(cases));
end
fprintf('Test_tensopmat_clean_op: all %d cases passed.\n',numel(cases));
end

%% --------------------------------------------------------------------------
function case_basic()
% case_basic -- tests that @tensopmat/clean_op.m cleans
% every (nonempty) block by dispatching to @tensopvar/clean_op.m,
% correctly determines whether the WHOLE block matrix is zero, and never
% empties a zeroed block (since Cop.ops{i,j}=[] is reserved for "block not
% provided", and initialize.m requires at least one nonempty block per row
% and per column). Swept over varying numbers of row/column blocks and
% spatial variables, and over single- vs multi-term blocks.
% Oracle: apply_op of the cleaned vs. the uncleaned operator on the same
% random input.

    tol = 1e-8;

    for trial = 1:10
        rng(6600+trial);

        p = randi(3)+1;   d = randi(3)+1;
        N = randi(3)+1;   M = randi(3)+1;
        nterms_max = 1 + mod(trial,2);   % alternate single-term / multi-term blocks

        Pop = rand_tensopmat(p,d,N,M,[0,1],0.7,0.8,struct(),'s','t',1,nterms_max,2);
        [p,d] = size(Pop.ops);
        populated = find(~cellfun(@isempty,Pop.ops));
        if numel(populated)<2
            continue   % need at least 2 populated blocks to test a non-trivial drop
        end

        % Input covering every column, from the first nonempty block in each
        X = cell(1,d);
        for j=1:d
            i0 = find(~cellfun(@isempty,Pop.ops(:,j)),1);
            vars_ij = Pop.ops{i0,j}.vars.in;
            X{j} = rand_poly([Pop.dims.in(j),1],polynomial(vars_ij(:)),2,20);
        end

        %%%% Case 1: zero out exactly one populated block -- not flagged zero,
        %%%% action unaffected.
        idx1 = populated(1);
        [i1,j1] = ind2sub([p,d],idx1);
        Pop_bad = Pop;
        Pop_bad.ops{i1,j1} = 0*Pop_bad.ops{i1,j1};
        [Pop_clean,is_zero1] = clean_op(Pop_bad);
        if is_zero1
            error("Test_tensopmat_clean_op failed: operator incorrectly flagged entirely zero after zeroing only one block (trial=%d).",trial)
        elseif isempty(Pop_clean.ops{i1,j1})
            error("Test_tensopmat_clean_op failed: a zeroed block was emptied out (trial=%d).",trial)
        end
        Y_bad = apply_op(Pop_bad,X);
        Y_clean = apply_op(Pop_clean,X);
        err = cleanpoly(Y_bad-Y_clean,1e-10);
        if max(max(abs(full(err.C))))>tol
            error("Test_tensopmat_clean_op failed: cleaned operator does not act identically (trial=%d).",trial)
        end

        %%%% Case 2: zero out EVERY populated block -- is_zero true, zero map.
        Pop_allzero = Pop;
        for k=1:numel(populated)
            [ik,jk] = ind2sub([p,d],populated(k));
            Pop_allzero.ops{ik,jk} = 0*Pop_allzero.ops{ik,jk};
        end
        [Pop_allzero_clean,is_zero2] = clean_op(Pop_allzero);
        if ~is_zero2
            error("Test_tensopmat_clean_op failed: entirely-zero operator not flagged is_zero (trial=%d).",trial)
        end
        for k=1:numel(populated)
            [ik,jk] = ind2sub([p,d],populated(k));
            if isempty(Pop_allzero_clean.ops{ik,jk})
                error("Test_tensopmat_clean_op failed: a zeroed block was emptied out in the all-zero case (trial=%d).",trial)
            end
        end
        Y_allzero = apply_op(Pop_allzero_clean,X);
        if max(max(abs(full(double(cleanpoly(Y_allzero,1e-10).C)))))>tol
            error("Test_tensopmat_clean_op failed: all-zero case did not act as the zero operator (trial=%d).",trial)
        end

        %%%% Case 3: nothing zero -- is_zero false, action unchanged.
        [Pop_clean3,is_zero3] = clean_op(Pop);
        if is_zero3
            error("Test_tensopmat_clean_op failed: nonzero operator incorrectly flagged is_zero (trial=%d).",trial)
        end
        Y0 = apply_op(Pop,X);
        Y3 = apply_op(Pop_clean3,X);
        err3 = cleanpoly(Y0-Y3,1e-10);
        if max(max(abs(full(err3.C))))>tol
            error("Test_tensopmat_clean_op failed: cleaning an already-clean operator changed its action (trial=%d).",trial)
        end

        %%%% Case 4: sub-tolerance noise in one block's sopvar coefficient is
        %%%% scrubbed, other blocks untouched.
        ztol = 1e-6;
        noise = ztol/1000;
        idx4 = populated(end);
        [i4,j4] = ind2sub([p,d],idx4);
        Pop_noisy = Pop;
        Pop_noisy.ops{i4,j4}.ops{1,1}.params{1}(1,1) = noise;
        [Pop_clean4,~] = clean_op(Pop_noisy,ztol);
        if Pop_clean4.ops{i4,j4}.ops{1,1}.params{1}(1,1)~=0
            error("Test_tensopmat_clean_op failed: sub-tolerance noise in a block was not zeroed (trial=%d).",trial)
        end
        for k=1:numel(populated)
            if populated(k)==idx4
                continue
            end
            [ik,jk] = ind2sub([p,d],populated(k));
            if ~isequal(Pop_clean4.ops{ik,jk},Pop.ops{ik,jk})
                error("Test_tensopmat_clean_op failed: cleaning noise in one block changed an untouched block (trial=%d).",trial)
            end
        end
    end

    disp('Test_tensopmat_clean_op passed.')
end

%% --------------------------------------------------------------------------
function case_dead_vars()
% case_dead_vars -- tests that @tensopmat/clean_op.m
% drops an input/output variable that no block structurally depends on (an
% entirely-zero column of depmat.in/depmat.out), leaving every used
% variable and the operator's action untouched; and that a 'tensopmat'
% with no dead variables is left unchanged. Swept over varying numbers of
% row/column blocks and spatial variables.
% Oracle: the original vars/dom/depmat and apply_op of the original.


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
            error("case_dead_vars failed: vars/dom/depmat changed with no dead variables present (trial=%d).",trial)
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
            error("case_dead_vars failed: operator incorrectly flagged entirely zero (trial=%d).",trial)
        end
        if any(strcmp(Pop_clean.vars.out,'dead_out')) || any(strcmp(Pop_clean.vars.in,'dead_in'))
            error("case_dead_vars failed: dead variable not removed (trial=%d).",trial)
        end
        if ~isequal(Pop_clean.vars,Pop.vars) || ~isequal(Pop_clean.dom,Pop.dom) || ~isequal(Pop_clean.depmat,Pop.depmat)
            error("case_dead_vars failed: cleaning the dead-variable object did not reduce exactly to the original vars/dom/depmat (trial=%d).",trial)
        end

        % Action must be unaffected (apply_op never reads vars/dom/depmat, so
        % this holds trivially)
        Y_dead = apply_op(Pop_dead,X);
        Y_clean = apply_op(Pop_clean,X);
        err1 = cleanpoly(Y_dead-Y0,1e-10);
        err2 = cleanpoly(Y_clean-Y0,1e-10);
        if max(max(abs(full(err1.C))))>tol || max(max(abs(full(err2.C))))>tol
            error("case_dead_vars failed: operator's action changed (trial=%d).",trial)
        end
    end

    disp('case_dead_vars passed.')
end
