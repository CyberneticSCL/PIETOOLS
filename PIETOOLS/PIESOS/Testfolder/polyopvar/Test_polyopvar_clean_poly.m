% Test_polyopvar_clean_poly.m -- verify @polyopvar/clean_poly.m:
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

% DJ (CC), 10/15/2026: Initial coding;

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
