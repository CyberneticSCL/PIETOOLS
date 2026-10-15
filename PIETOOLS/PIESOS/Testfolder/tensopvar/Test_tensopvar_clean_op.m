% Test_tensopvar_clean_op.m -- verify @tensopvar/clean_op.m correctly
% removes terms (rows of Cop.ops) that are zero up to tolerance, leaves
% nonzero terms and the represented operator untouched, and always retains
% at least one term so matdim/vars/dom/depmat stay well-defined.
%
% The oracle is the action on random inputs (apply_op) of the operator
% before and after cleaning.
%
% DJ (CC), 10/15/2026: Initial coding;

tol = 1e-8;

for trial = 1:8
    rng(6400+trial);

    dims = [randi(2),randi(2); randi(2),randi(2)];
    vars = struct('in',{{'s1','s2'}},'out',{{'s3','s4'}});
    dom = struct('in',repmat([0,1],2,1),'out',repmat([0,1],2,1));
    depmat = struct('in',logical(eye(2)),'out',logical(eye(2)));
    Bop = rand_tensopvar(dims,vars,dom,depmat,3,[2,2],1,0.5);

    % Case 1: zero out one term's first factor: that term is dropped, the action is unchanged
    Bop_bad = Bop;
    Bop_bad.ops{2,1} = 0*Bop_bad.ops{2,1};
    Bop_clean = clean_op(Bop_bad);
    if size(Bop_clean.ops,1)~=2
        error("Test_tensopvar_clean_op failed: expected 2 surviving terms, got %d (trial=%d).",size(Bop_clean.ops,1),trial)
    end
    if ~isequal(Bop_clean.matdim,Bop_bad.matdim) || ~isequal(Bop_clean.vars,Bop_bad.vars) || ...
            ~isequal(Bop_clean.dom,Bop_bad.dom) || ~isequal(Bop_clean.depmat,Bop_bad.depmat)
        error("Test_tensopvar_clean_op failed: matdim/vars/dom/depmat changed after cleaning (trial=%d).",trial)
    end
    pvars = polynomial(Bop.vars.in(:));
    v = {rand_poly([dims(1,2),1],pvars(1),3,10), rand_poly([dims(2,2),1],pvars(2),3,10)};
    Y_bad = apply_op(Bop_bad,v);
    Y_clean = apply_op(Bop_clean,v);
    err = cleanpoly(Y_bad-Y_clean,1e-10);
    if max(max(abs(err.C)))>tol
        error("Test_tensopvar_clean_op failed: cleaned operator does not act identically (trial=%d).",trial)
    end

    % Case 2: every term zero: exactly 1 term retained, acts as zero map
    Bop_allzero = Bop;
    for i=1:3
        Bop_allzero.ops{i,1} = 0*Bop_allzero.ops{i,1};
    end
    Bop_allzero_clean = clean_op(Bop_allzero);
    if size(Bop_allzero_clean.ops,1)~=1
        error("Test_tensopvar_clean_op failed: all-zero case should retain exactly 1 term, got %d (trial=%d).",size(Bop_allzero_clean.ops,1),trial)
    end
    if ~isequal(Bop_allzero_clean.matdim,Bop_allzero.matdim)
        error("Test_tensopvar_clean_op failed: matdim not preserved in all-zero case (trial=%d).",trial)
    end
    Y_allzero = apply_op(Bop_allzero_clean,v);
    if max(max(abs(double(cleanpoly(Y_allzero,1e-10).C))))>tol
        error("Test_tensopvar_clean_op failed: all-zero case did not act as the zero operator (trial=%d).",trial)
    end

    % Case 3: nothing zero: no term removed
    Bop_clean2 = clean_op(Bop);
    if size(Bop_clean2.ops,1)~=3
        error("Test_tensopvar_clean_op failed: no-op case incorrectly removed a term (trial=%d).",trial)
    end
end

disp('Test_tensopvar_clean_op passed.')
