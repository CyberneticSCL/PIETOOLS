function Test_polyopvar_clean_poly()
% Test_polyopvar_clean_poly.m -- one local function for each case:
%   case_basic              verify @polyopvar/clean_poly.m: (1)
%                           sub-tolerance coefficients are scrubbed without
%                           affecting other blocks or the polynomial's
%                           action;
%   case_canonical          verify that @polyopvar/clean_poly.m gives the
%                           same functional for equal functionals, with
%                           coefficients of type 'intop': the same
%                           monomials,
%   case_idempotent         verify that @polyopvar/clean_poly.m is
%                           idempotent on functionals of repeated states,
%                           clean_poly(clean_poly(f)) = clean_poly(f), and
%                           does not change
%
% DJ (CC), 10/15/2026: Initial coding;
% DJ (CC), 10/23/2026: Merge Test_polyopvar_clean_poly_canonical,
%   Test_polyopvar_clean_poly_idempotent into this test, as cases;

cases = {@case_basic, @case_canonical, @case_idempotent};
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
    error('Test_polyopvar_clean_poly: %d of %d cases failed.',nfail,numel(cases));
end
fprintf('Test_polyopvar_clean_poly: all %d cases passed.\n',numel(cases));
end

%% --------------------------------------------------------------------------
function case_basic()
% case_basic -- verify @polyopvar/clean_poly.m:
% (1) sub-tolerance coefficients are scrubbed without affecting other
%     blocks or the polynomial's action;
% (2) a monomial that does not contribute (all its blocks zero) is
%     dropped from Z/C, and a state variable left unused as a result is
%     dropped from varname/dims.in/depmat.in, with Z's indices remapped;
% (3) an entirely-zero polynomial is flagged is_zero, acts as the zero
%     map, and keeps exactly one (zeroed) monomial/state so its shape
%     stays well-defined;
% (4) an already-clean polynomial is left completely unchanged;
% (5)-(6) the same monomial-drop behavior holds for a degree-2
%     (self-product) monomial, both alone and alongside a surviving
%     degree-1 monomial of the same state;
% (7)-(8) the combine_terms(f) step: a duplicated monomial row is merged
%     (summed), doubling its contribution; and a duplicate whose copies are
%     negatives of each other sums to a NEW zero monomial that clean_poly's
%     own pruning then drops.
% Oracle: structure checked directly (Z, varname, dims, depmat) and values
% by apply_poly against the numeric polynomials.

    tol = 1e-8;

    for trial = 1:8
        rng(6800+trial);

        pvar s1 s2
        x = polyopvar('x',[s1,s2],[0,1]);
        y = polyopvar('y',[s1,s2],[0,1]);
        f = x+y;   % 2 monomials, degree 1 each, states 'x' and 'y'

        xv = rand_poly([1,1],{'s1','s2'},2,5);
        yv = rand_poly([1,1],{'s1','s2'},2,5);

        %%%% Case 1: sub-tolerance noise in a zero parameter slot is scrubbed;
        %%%% nothing else changes and the action is unaffected. The slot is
        %%%% exactly zero beforehand so the block stays nonzero (no monomial dropped).
        ztol = 1e-6;
        noise = ztol/1000;
        f_noisy = f;
        [irow,jcol] = find(~cellfun(@isempty,f_noisy.C.ops));
        pick = randi(numel(irow));
        blk = f_noisy.C.ops{irow(pick),jcol(pick)};
        sop = blk.ops{1,1};
        kzero = find(cellfun(@(c) nnz(c)==0, sop.params),1);
        sop.params{kzero}(1,1) = noise;
        blk.ops{1,1} = sop;
        f_noisy.C.ops{irow(pick),jcol(pick)} = blk;

        [f_clean1,is_zero1] = clean_poly(f_noisy,ztol);
        if is_zero1
            error("Test_polyopvar_clean_poly failed: noisy polynomial incorrectly flagged is_zero (trial=%d).",trial)
        elseif size(f_clean1.Z,1)~=size(f.Z,1) || numel(f_clean1.varname)~=numel(f.varname)
            error("Test_polyopvar_clean_poly failed: sub-tolerance noise caused a structural change (trial=%d).",trial)
        elseif f_clean1.C.ops{irow(pick),jcol(pick)}.ops{1,1}.params{kzero}(1,1)~=0
            error("Test_polyopvar_clean_poly failed: sub-tolerance noise not zeroed (trial=%d).",trial)
        end
        Y_noisy = apply_poly(f_noisy,{xv,yv});
        Y_clean1 = apply_poly(f_clean1,{xv,yv});
        err1 = cleanpoly(Y_noisy-Y_clean1,1e-10);
        if max(max(abs(full(err1.C))))>tol
            error("Test_polyopvar_clean_poly failed: cleaning sub-tolerance noise changed the polynomial's action (trial=%d).",trial)
        end

        %%%% Case 2: zero out the 'y' monomial: it and the unused state 'y' are dropped.
        f_half = f;
        y_idx = find(strcmp(f_half.varname,'y'));
        Zmat = full(f_half.Z);
        row_y = find(Zmat(:,1)==y_idx);
        for i=1:size(f_half.C.ops,1)
            if ~isempty(f_half.C.ops{i,row_y})
                f_half.C.ops{i,row_y} = 0*f_half.C.ops{i,row_y};
            end
        end

        [f_half_clean,is_zero2] = clean_poly(f_half);
        if is_zero2
            error("Test_polyopvar_clean_poly failed: half-zero polynomial incorrectly flagged is_zero (trial=%d).",trial)
        elseif size(f_half_clean.Z,1)~=1
            error("Test_polyopvar_clean_poly failed: zeroed monomial not dropped (trial=%d).",trial)
        elseif numel(f_half_clean.varname)~=1 || ~strcmp(f_half_clean.varname{1},'x')
            error("Test_polyopvar_clean_poly failed: unused state 'y' not dropped (trial=%d).",trial)
        end
        Y_half_clean = apply_poly(f_half_clean,{xv});
        err2 = cleanpoly(Y_half_clean-xv,1e-10);
        if max(max(abs(full(err2.C))))>tol
            error("Test_polyopvar_clean_poly failed: cleaned half-zero polynomial does not match the surviving 'x' term (trial=%d).",trial)
        end

        %%%% Case 3: entirely-zero polynomial: is_zero, zero map, one zeroed monomial/state kept.
        f_zero = f;
        for k=1:numel(f_zero.C.ops)
            if ~isempty(f_zero.C.ops{k})
                f_zero.C.ops{k} = 0*f_zero.C.ops{k};
            end
        end
        [f_zero_clean,is_zero3] = clean_poly(f_zero);
        if ~is_zero3
            error("Test_polyopvar_clean_poly failed: entirely-zero polynomial not flagged is_zero (trial=%d).",trial)
        elseif size(f_zero_clean.Z,1)~=1 || numel(f_zero_clean.varname)~=1
            error("Test_polyopvar_clean_poly failed: all-zero case should retain exactly 1 monomial/state (trial=%d).",trial)
        end
        xv_zero = {xv,yv};
        Y_zero = apply_poly(f_zero_clean,xv_zero(1:numel(f_zero_clean.varname)));
        if max(max(abs(full(double(cleanpoly(Y_zero,1e-10).C)))))>tol
            error("Test_polyopvar_clean_poly failed: entirely-zero polynomial did not act as the zero map (trial=%d).",trial)
        end

        %%%% Case 4: already clean: no structural change.
        [f_clean4,is_zero4] = clean_poly(f);
        if is_zero4
            error("Test_polyopvar_clean_poly failed: nonzero polynomial incorrectly flagged is_zero (trial=%d).",trial)
        elseif ~isequal(f_clean4.Z,f.Z) || ~isequal(f_clean4.varname,f.varname) || ...
                ~isequal(f_clean4.dims.in,f.dims.in) || ~isequal(f_clean4.depmat.in,f.depmat.in)
            error("Test_polyopvar_clean_poly failed: cleaning an already-clean polynomial changed its structure (trial=%d).",trial)
        end
        Y0 = apply_poly(f,{xv,yv});
        Y4 = apply_poly(f_clean4,{xv,yv});
        err4 = cleanpoly(Y0-Y4,1e-10);
        if max(max(abs(full(err4.C))))>tol
            error("Test_polyopvar_clean_poly failed: cleaning an already-clean polynomial changed its action (trial=%d).",trial)
        end

        %%%% Case 5: degree-2 monomial x.*x, entirely zeroed.
        g = x.*x;
        for i=1:size(g.C.ops,1)
            if ~isempty(g.C.ops{i,1})
                g.C.ops{i,1} = 0*g.C.ops{i,1};
            end
        end
        [g_clean,is_zero5] = clean_poly(g);
        if ~is_zero5
            error("Test_polyopvar_clean_poly failed: zeroed degree-2 polynomial not flagged is_zero (trial=%d).",trial)
        elseif size(g_clean.Z,1)~=1 || numel(g_clean.varname)~=1
            error("Test_polyopvar_clean_poly failed: zeroed degree-2 case should retain exactly 1 monomial/state (trial=%d).",trial)
        end
        Yg = apply_poly(g_clean,{xv});
        if max(max(abs(full(double(cleanpoly(Yg,1e-10).C)))))>tol
            error("Test_polyopvar_clean_poly failed: zeroed degree-2 polynomial did not act as the zero map (trial=%d).",trial)
        end

        %%%% Case 6: 'x' + 'x.*x', zero out the degree-2 monomial: it is dropped,
        %%%% 'x' (used by the degree-1 monomial) is kept.
        h = x + x.*x;
        Zmat_h = full(h.Z);
        row_deg2 = find(sum(Zmat_h>0,2)==2);
        for i=1:size(h.C.ops,1)
            if ~isempty(h.C.ops{i,row_deg2})
                h.C.ops{i,row_deg2} = 0*h.C.ops{i,row_deg2};
            end
        end
        [h_clean,is_zero6] = clean_poly(h);
        if is_zero6
            error("Test_polyopvar_clean_poly failed: half-zero mixed-degree polynomial incorrectly flagged is_zero (trial=%d).",trial)
        elseif size(h_clean.Z,1)~=1
            error("Test_polyopvar_clean_poly failed: zeroed degree-2 monomial not dropped (trial=%d).",trial)
        elseif numel(h_clean.varname)~=1 || ~strcmp(h_clean.varname{1},'x')
            error("Test_polyopvar_clean_poly failed: still-used state 'x' incorrectly dropped (trial=%d).",trial)
        end
        Yh = apply_poly(h_clean,{xv});
        errh = cleanpoly(Yh-xv,1e-10);
        if max(max(abs(full(errh.C))))>tol
            error("Test_polyopvar_clean_poly failed: cleaned mixed-degree polynomial does not match the surviving degree-1 term (trial=%d).",trial)
        end

        %%%% Case 7: duplicate the 'x' monomial by hand (copied, not negated):
        %%%% combine_terms sums the copies, doubling its contribution; nothing dropped.
        f_dup = f;
        x_idx = find(strcmp(f_dup.varname,'x'));
        Zmat7 = full(f_dup.Z);
        row_x7 = find(Zmat7(:,1)==x_idx);
        f_dup.Z = [f_dup.Z; f_dup.Z(row_x7,:)];
        f_dup.C.ops = [f_dup.C.ops, f_dup.C.ops(:,row_x7)];
        f_dup.C.dims.in = [f_dup.C.dims.in; f_dup.C.dims.in(row_x7)];
        f_dup.C.depmat.in = [f_dup.C.depmat.in; f_dup.C.depmat.in(row_x7,:)];
        f_dup.chkval = false;
        f_dup.C.chkval = false;

        [f_dup_clean,is_zero7] = clean_poly(f_dup);
        if is_zero7
            error("Test_polyopvar_clean_poly failed: duplicated-monomial polynomial incorrectly flagged is_zero (trial=%d).",trial)
        elseif size(f_dup_clean.Z,1)~=size(f.Z,1) || numel(f_dup_clean.varname)~=numel(f.varname)
            error("Test_polyopvar_clean_poly failed: duplicate monomial not merged back down to the original structure (trial=%d).",trial)
        end
        Y_dup = apply_poly(f_dup,{xv,yv});
        Y_dup_clean = apply_poly(f_dup_clean,{xv,yv});
        err7a = cleanpoly(Y_dup-(2*xv+yv),1e-10);
        if max(max(abs(full(err7a.C))))>tol
            error("Test_polyopvar_clean_poly failed: hand-duplicated setup did not read as 2*x+y (trial=%d).",trial)
        end
        err7b = cleanpoly(Y_dup_clean-Y_dup,1e-10);
        if max(max(abs(full(err7b.C))))>tol
            error("Test_polyopvar_clean_poly failed: merging the duplicate monomial changed the polynomial's action (trial=%d).",trial)
        end

        %%%% Case 8: duplicate the 'x' monomial with the second copy negated:
        %%%% combine_terms sums them to a new zero monomial, which clean_poly drops.
        f_cancel = f;
        Zmat8 = full(f_cancel.Z);
        row_x8 = find(Zmat8(:,1)==x_idx);
        f_cancel.Z = [f_cancel.Z; f_cancel.Z(row_x8,:)];
        f_cancel.C.ops = [f_cancel.C.ops, f_cancel.C.ops(:,row_x8)];
        for i=1:size(f_cancel.C.ops,1)
            if ~isempty(f_cancel.C.ops{i,end})
                f_cancel.C.ops{i,end} = -f_cancel.C.ops{i,end};
            end
        end
        f_cancel.C.dims.in = [f_cancel.C.dims.in; f_cancel.C.dims.in(row_x8)];
        f_cancel.C.depmat.in = [f_cancel.C.depmat.in; f_cancel.C.depmat.in(row_x8,:)];
        f_cancel.chkval = false;
        f_cancel.C.chkval = false;

        [f_cancel_clean,is_zero8] = clean_poly(f_cancel);
        if is_zero8
            error("Test_polyopvar_clean_poly failed: cancelling-duplicate polynomial incorrectly flagged is_zero (trial=%d).",trial)
        elseif size(f_cancel_clean.Z,1)~=1
            error("Test_polyopvar_clean_poly failed: cancelled monomial (x + x - x) not dropped (trial=%d).",trial)
        elseif numel(f_cancel_clean.varname)~=1 || ~strcmp(f_cancel_clean.varname{1},'y')
            error("Test_polyopvar_clean_poly failed: unused state 'x' not dropped after cancellation (trial=%d).",trial)
        end
        Y_cancel_clean = apply_poly(f_cancel_clean,{yv});
        err8 = cleanpoly(Y_cancel_clean-yv,1e-10);
        if max(max(abs(full(err8.C))))>tol
            error("Test_polyopvar_clean_poly failed: cleaned cancelling-duplicate polynomial does not match the surviving 'y' term (trial=%d).",trial)
        end
    end

    disp('Test_polyopvar_clean_poly passed.')
end

%% --------------------------------------------------------------------------
function case_canonical()
% case_canonical -- verify that
% @polyopvar/clean_poly.m gives the same functional for equal
% functionals, with coefficients of type 'intop': the same monomials,
% intmat, and kernels, and that the difference of equal functionals is
% recognized as zero (is_zero).
%
% Cases: (1) x o x with the full square [0,3,0,3] and kernel k, with
% kernel k(s2,s1), the chain [0,3,0,1] with k(s1,s2)+k(s2,s1), and the
% chain [0,2,0,3] with this kernel; (2) x o y for two states, where the
% square is the sum of the chains [0,3,0,1] and [0,2,0,3] with the same
% kernel, which cannot be relabeled; (3) a quartic functional in x with
% the full cube, and the same functional with the copies of x relabeled;
% (4) the same in two spatial variables, relabeling the copies in both;
% (5) a functional that is not equal is not zero; (6) x o x with a
% multiplier [0,3,1,1], whose kernel k(s1,s2) is the same functional as
% k(s1,s1), and as the same kernel on [2,2,0,3]; (7) a vector-valued
% state, int x'*M*x with M not symmetric: the kernels on the elements
% (1,2) and (2,1) are interchangeable, only their sum matters.
%

% NOTES:
% - Equal functionals must have identical clean_poly output because
%   piesos_eq equates the kernels cell by cell; same_poly compares
%   monomials, intmat, kernels.
% - Values are evaluated by eval_functional below, not apply_poly, which
%   does not substitute states in functionals of degree >1.

    tol = 1e-8;

    %% Case 1: x o x, one spatial variable
    pvar s s_dum1 s_dum2
    x = polyopvar('x1',s,[0,1]);
    k = s_dum1^2+3*s_dum1*s_dum2^3;
    ksw = subs(k,[s_dum1;s_dum2],[s_dum2;s_dum1]);
    mk = @(K,im) set_coeff(x*x,intop({K},{im},{'s'},[0,1]));
    forms = {mk(k,[0,3,0,3]), mk(ksw,[0,3,0,3]), mk(k+ksw,[0,3,0,1]), mk(k+ksw,[0,2,0,3])};
    xv = rand_poly([1,1],{'s'},2,4);
    vals = zeros(1,4);
    for j=1:4
        vals(j) = double(cleanpoly(eval_functional(forms{j},{xv}),1e-12));
    end
    if max(vals)-min(vals)>tol*max(1,max(abs(vals)))
        error("case_canonical failed: Case 1, the four forms are not equal functionals.")
    end
    for j=2:4
        if ~same_poly(clean_poly(forms{1}),clean_poly(forms{j}))
            error("case_canonical failed: Case 1, form %d gives another clean_poly than form 1.",j)
        end
        [~,is_zero] = clean_poly(forms{1}-forms{j});
        if ~is_zero
            error("case_canonical failed: Case 1, form 1 - form %d is not zero.",j)
        end
    end

    %% Case 2: x o y, two states: the square is the sum of two chains with the same kernel
    y = polyopvar('y1',s,[0,1]);
    fsq = set_coeff(x*y,intop({k},{[0,3,0,3]},{'s'},[0,1]));
    fch = set_coeff(x*y,intop({k;k},{[0,3,0,1;0,2,0,3]},{'s'},[0,1]));
    fo = set_coeff(x*y,intop({k},{[0,3,0,1]},{'s'},[0,1]));
    [~,is_zero] = clean_poly(fsq-fch);
    if ~is_zero
        error("case_canonical failed: Case 2, the square is not the sum of the chains.")
    end
    [~,is_zero] = clean_poly(fsq-fo);
    if is_zero
        error("case_canonical failed: Case 2, the square is not one chain for two states.")
    end
    if ~same_poly(clean_poly(fsq),clean_poly(fch))
        error("case_canonical failed: Case 2, other form after clean_poly.")
    end

    %% Case 3: quartic functional in x: the cube, with the copies of x relabeled
    pvar s_dum3 s_dum4
    vv = [s_dum1;s_dum2;s_dum3;s_dum4];
    rng(1);
    k4 = rand_poly([1,1],{'s_dum1','s_dum2','s_dum3','s_dum4'},2,6);
    perm = [3,1,4,2];
    k4p = subs(k4,vv,vv(perm));
    mk4 = @(K) set_coeff(x*x*x*x,intop({K},{[0,5,0,5,0,5,0,5]},{'s'},[0,1]));
    f4 = mk4(k4);       g4 = mk4(k4p);
    [~,is_zero] = clean_poly(f4-g4);
    if ~is_zero
        error("case_canonical failed: Case 3, relabeled quartic functional is not equal.")
    end
    if ~same_poly(clean_poly(f4),clean_poly(g4))
        error("case_canonical failed: Case 3, relabeled quartic functional gives another clean_poly.")
    end
    xv = rand_poly([1,1],{'s'},1,3);
    err = cleanpoly(eval_functional(f4,{xv})-eval_functional(clean_poly(f4),{xv}),1e-10);
    if ~isdouble(err) || abs(double(err))>tol*max(1,abs(double(cleanpoly(eval_functional(f4,{xv}),1e-10))))
        error("case_canonical failed: Case 3, clean_poly changed the value.")
    end

    %% Case 4: two spatial variables, copies relabeled in both
    pvar s1 s2
    x2 = polyopvar('x2',[s1,s2],[0,1]);
    dn = {'s1_dum1','s1_dum2','s2_dum1','s2_dum2'};
    w = polynomial(dn(:));
    kk2 = rand_poly([1,1],dn,2,6);
    kk2p = subs(kk2,w,w([2,1,4,3]));
    mk2 = @(K,I1,I2) set_coeff(x2*x2,intop({K},{I1,I2},{'s1','s2'},[0,1;0,1]));
    f2 = mk2(kk2,[0,3,0,3],[0,3,0,3]);
    g2 = mk2(kk2p,[0,3,0,3],[0,3,0,3]);
    [~,is_zero] = clean_poly(f2-g2);
    if ~is_zero
        error("case_canonical failed: Case 4, relabeled functional is not equal.")
    end
    if ~same_poly(clean_poly(f2),clean_poly(g2))
        error("case_canonical failed: Case 4, relabeled functional gives another clean_poly.")
    end
    % [2,3,0,3] along s1 is only one of the two chains of the square: another functional
    h2 = mk2(kk2p,[2,3,0,3],[0,3,0,3]);
    [~,is_zero] = clean_poly(f2-h2);
    if is_zero
        error("case_canonical failed: Case 4, different functionals are equal.")
    end

    %% Case 5: a functional that is not equal is not zero
    [~,is_zero] = clean_poly(forms{1}-mk(k,[0,3,0,1]));
    if is_zero
        error("case_canonical failed: Case 5, a different functional was recognized as zero.")
    end

    %% Case 6: multiplier: the kernel k(s1,s2) on [0,3,1,1] (s2=s1) is the functional with kernel k(s1,s1)
    kd = subs(k,s_dum2,s_dum1);
    f6 = mk(k,[0,3,1,1]);
    g6 = mk(kd,[0,3,1,1]);
    h6 = mk(k,[2,2,0,3]);                   % s1=s2, with s2 the integrated variable
    [~,is_zero] = clean_poly(f6-g6);
    if ~is_zero
        error("case_canonical failed: Case 6, k(s1,s2) and k(s1,s1) on [0,3,1,1] are not equal.")
    end
    [~,is_zero] = clean_poly(f6-h6);
    if ~is_zero
        error("case_canonical failed: Case 6, [0,3,1,1] and [2,2,0,3] are not equal.")
    end
    if ~same_poly(clean_poly(f6),clean_poly(g6)) || ~same_poly(clean_poly(f6),clean_poly(h6))
        error("case_canonical failed: Case 6, other form after clean_poly.")
    end
    [~,is_zero] = clean_poly(f6-mk(kd+1,[0,3,1,1]));
    if is_zero
        error("case_canonical failed: Case 6, a different functional was recognized as zero.")
    end

    %% Case 7: vector-valued state with two elements, int x'*M*x: the weights M(1,2) and M(2,1) are interchangeable
    xv2 = polyopvar('y1',2,s,[0,1]);
    M = [s+1, 2*s^2; -3*s, 5];                          % not symmetric
    Msym = [s+1, (2*s^2-3*s)/2; (2*s^2-3*s)/2, 5];      % the same functional
    Mdif = [s+1, 2*s^2; 2*s^2, 5];                      % another functional
    f7 = innerprod(xv2,M,xv2);
    g7 = innerprod(xv2,Msym,xv2);
    h7 = innerprod(xv2,Mdif,xv2);
    [~,is_zero] = clean_poly(f7-g7);
    if ~is_zero
        error("case_canonical failed: Case 7, asymmetric and symmetric weights of x'*M*x are not equal.")
    end
    if ~same_poly(clean_poly(f7),clean_poly(g7))
        error("case_canonical failed: Case 7, other form after clean_poly.")
    end
    [~,is_zero] = clean_poly(f7-h7);
    if is_zero
        error("case_canonical failed: Case 7, a different functional was recognized as zero.")
    end
    disp('case_canonical passed.')
end

%% --------------------------------------------------------------------------
function case_idempotent()
% case_idempotent -- verify that
% @polyopvar/clean_poly.m is idempotent on functionals of repeated
% states, clean_poly(clean_poly(f)) = clean_poly(f), and does not change
% the functional.
%
% For each of the functionals f below, with coefficients of type 'intop':
% (1) clean_poly is applied four times; the results of the second, third
%     and fourth call must be identical to that of the first: the same
%     monomials Z, the same intmat, monomials Z.in, and parameters of
%     every functional (so also the same kernels() and intmat);
% (2) the cleaned functionals must have the same value as f at random
%     polynomial states (apply_op of the coefficients on copies of the
%     states; apply_poly does not substitute states in functionals of
%     degree >1).
% Cases: x o x with the tied full square [0,3,0,3] and with the two
% orderings (chains) of the variables, a quartic functional in x (4
% factors, with ties), x o x in two spatial variables (x(s1,s2)), a
% quartic functional in two spatial variables, and a functional of the
% monomial x o y o x o y of two states.
%

% NOTES:
% - The test checks that clean_poly does not swap the dummy variables
%   tied in the order of integration (full square [0,3,0,3]) on every
%   call (@intop/combine_terms.m); piesos_eq and piesos_ineq compare
%   kernels without symmetrizing.

    tol = 1e-8;

    %% Case 1: x o x in one spatial variable, tied and ordered variables
    pvar s s_dum1 s_dum2
    x = polyopvar('x1',s,[0,1]);
    xv = rand_poly([1,1],{'s'},2,5);
    intmat_arr = {[0,3,0,3],[0,3,0,1],[0,2,0,3]};
    for j=1:numel(intmat_arr)
        f = set_coeff(x*x,intop({s_dum1*(1+s_dum2^2)},intmat_arr(j),{'s'},[0,1]));
        check_idempotent(f,{xv},tol,sprintf('Case 1, intmat %s',mat2str(intmat_arr{j})))
    end

    %% Case 2: quartic functional in x, with ties, (x o x o x o x), 3 terms
    rng(1);
    dumvars = {'s_dum1','s_dum2','s_dum3','s_dum4'};
    intmat = [0,5,0,5,0,5,0,5; 0,5,0,5,0,1,0,2; 0,2,0,3,0,4,0,5];
    kern = cell(3,1);
    for k=1:3
        kern{k} = rand_poly([1,1],dumvars,1,6);
    end
    f = set_coeff(x*x*x*x,intop(kern,{intmat},{'s'},[0,1]));
    xv = rand_poly([1,1],{'s'},1,3);
    check_idempotent(f,{xv},tol,'Case 2')

    %% Case 3: x o x in two spatial variables, ordered differently along s1 and s2
    pvar s1 s2
    x2 = polyopvar('x2',[s1,s2],[0,1]);
    dumvars = {'s1_dum1','s1_dum2','s2_dum1','s2_dum2'};
    xv2 = rand_poly([1,1],{'s1','s2'},2,5);
    % intmat along s1, along s2, and the number of kernels
    cases = {{[0,2,0,3;0,3,0,1],[0,3,0,1]}, ...                % two orders along s1, twisted along s2
             {[0,2,0,3;0,3,0,1;0,3,0,3],[0,3,0,1]}, ...        % plus a tied term along s1
             {[0,3,0,3],[0,3,0,1;0,2,0,3]}, ...                % tied along s1, orders along s2
             {[0,3,0,3],[0,3,0,3]}};                           % tied along both
    for j=1:numel(cases)
        nt = [size(cases{j}{1},1),size(cases{j}{2},1)];
        kern = cell(nt);
        for k=1:numel(kern)
            kern{k} = rand_poly([1,1],dumvars,2,6);
        end
        f = set_coeff(x2*x2,intop(kern,cases{j},{'s1','s2'},[0,1;0,1]));
        check_idempotent(f,{xv2},tol,sprintf('Case 3.%d',j))
    end

    %% Case 4: quartic functional in two spatial variables
    dumvars = {'s1_dum1','s1_dum2','s1_dum3','s1_dum4','s2_dum1','s2_dum2','s2_dum3','s2_dum4'};
    I1 = [0,5,0,5,0,5,0,5; 0,5,0,5,0,1,0,2];
    I2 = [0,5,0,5,0,1,0,5; 0,2,0,3,0,4,0,5];
    kern = cell(2,2);
    for k=1:numel(kern)
        kern{k} = rand_poly([1,1],dumvars,1,4);
    end
    f = set_coeff(x2*x2*x2*x2,intop(kern,{I1,I2},{'s1','s2'},[0,1;0,1]));
    xv2 = rand_poly([1,1],{'s1','s2'},1,3);
    check_idempotent(f,{xv2},tol,'Case 4')

    %% Case 5: x o y o x o y, two states, in one spatial variable
    y = polyopvar('y1',s,[0,1]);
    dumvars = {'s_dum1','s_dum2','s_dum3','s_dum4'};
    intmat = [0,5,0,5,0,5,0,5; 0,2,0,3,0,4,0,5; 0,5,0,5,0,1,0,2];
    kern = cell(3,1);
    for k=1:3
        kern{k} = rand_poly([1,1],dumvars,1,6);
    end
    f = x*y*x*y;
    f = set_coeff(f,intop(kern,{intmat},{'s'},[0,1]));
    xv = rand_poly([1,1],{'s'},1,3);
    yv = rand_poly([1,1],{'s'},1,3);
    check_idempotent(f,{xv,yv},tol,'Case 5')

    disp('case_idempotent passed.')
end

%% ---- helpers --------------------------------------------------------------

function V = set_coeff(V,K)
% Set the coefficient of the monomial of V to the functional K (as in piesos_sosvar)
V.C = tensopmat({K});
V.dims.out = V.C.dims.out;      V.vars.out = V.C.vars.out;
V.dom.out = V.C.dom.out;        V.depmat.out = V.C.depmat.out;
V.chkval = true;
end

function Y = eval_functional(f,xvals)
% Value of f = sum_j K_j (x o ... o x) at xvals: apply K_j to the product of the states, with
% in each factor the spatial variable v renamed to v_k, k = number of factors so far depending on v
Y = 0;
svars = f.vars.in;
for j=1:size(f.Z,1)
    cnt = zeros(1,numel(svars));
    xval = 1;
    for st = nonzeros(f.Z(j,:))'
        dep = f.depmat.in(st,:);
        cnt(dep) = cnt(dep)+1;
        newnames = arrayfun(@(i) [svars{i},'_',num2str(cnt(i))],find(dep),'UniformOutput',false);
        xval = kron(xval,subs(xvals{st},polynomial(svars(dep)'),polynomial(newnames(:))));
    end
    Y = Y + apply_op(f.C.ops{1,j},xval);
end
end

function tf = same_poly(A,B)
% Same monomials, and for each block the same intmat and kernels
tf = isequal(full(A.Z),full(B.Z)) && isequal(size(A.C.ops),size(B.C.ops));
if ~tf
    return
end
for j=1:numel(A.C.ops)
    a = A.C.ops{j};     b = B.C.ops{j};
    if isempty(a) && isempty(b)
        continue
    end
    if ~isa(a,'intop') || ~isa(b,'intop') || ~isequal(a.intmat,b.intmat) || ~isequal(size(a.params),size(b.params))
        tf = false;
        return
    end
    ka = kernels(a);    kb = kernels(b);
    for r=1:numel(ka)
        e = cleanpoly(ka{r}-kb{r},1e-10);
        if ~isdouble(e) || any(any(abs(double(e))>1e-9))
            tf = false;
            return
        end
    end
end
end

function check_idempotent(f,xvals,tol,label)
% clean_poly applied four times: calls 2-4 leave the result of call 1 unchanged,
% and the functional keeps its value
W = cell(1,4);
W{1} = clean_poly(f);
for k=2:4
    W{k} = clean_poly(W{k-1});
    if ~same_poly_idempotent(W{k},W{1})
        error("case_idempotent failed: %s, call %d of clean_poly changed the functional.",label,k)
    end
end
Y0 = eval_functional_idempotent(f,xvals);
Y1 = eval_functional_idempotent(W{1},xvals);
err = cleanpoly(Y0-Y1,1e-10);
nrmref = norm(double(cleanpoly(Y0,1e-10)),'fro');
if ~isdouble(err) || norm(double(err),'fro')>tol*max(1,nrmref)
    error("case_idempotent failed: %s, clean_poly changed the value of the functional.",label)
end
end

function Y = eval_functional_idempotent(f,xvals)
% Value of f = sum_j K_j (x o ... o x), K_j = f.C.ops{1,j} of type 'intop', at the states xvals,
% by applying K_j to the product of the states, with in each factor the spatial variable v renamed to v_k,
% k = number of factors so far depending on v (as in Test_intop_combine_terms_ND_random)
Y = 0;
svars = f.vars.in;
for j=1:size(f.Z,1)
    cnt = zeros(1,numel(svars));
    xval = 1;
    for st = nonzeros(f.Z(j,:))'
        dep = f.depmat.in(st,:);
        cnt(dep) = cnt(dep)+1;
        newnames = arrayfun(@(i) [svars{i},'_',num2str(cnt(i))],find(dep),'UniformOutput',false);
        xval = kron(xval,subs(xvals{st},polynomial(svars(dep)'),polynomial(newnames(:))));
    end
    Y = Y + apply_op(f.C.ops{1,j},xval);
end
end

function tf = same_poly_idempotent(A,B)
% Same monomials, and for each block the same intmat, monomials, parameters, and kernels
tf = isequal(full(A.Z),full(B.Z)) && isequal(size(A.C.ops),size(B.C.ops));
if ~tf
    return
end
for j=1:numel(A.C.ops)
    a = A.C.ops{j};     b = B.C.ops{j};
    if isempty(a) && isempty(b)
        continue
    end
    if ~isa(a,'intop') || ~isa(b,'intop') || ~same_op(a,b)
        tf = false;
        return
    end
    ka = kernels(a);    kb = kernels(b);
    for k=1:numel(ka)
        err = cleanpoly(ka{k}-kb{k},1e-12);
        if ~isdouble(err) || any(any(double(err)))
            tf = false;
            return
        end
    end
end
end

function tf = same_op(A,B)
tf = isequal(A.intmat,B.intmat) && isequal(A.Z.in,B.Z.in) && isequal(A.Zd,B.Zd) ...
    && isequal(size(A.params),size(B.params)) && all(cellfun(@isequal,A.params(:),B.params(:)));
end
