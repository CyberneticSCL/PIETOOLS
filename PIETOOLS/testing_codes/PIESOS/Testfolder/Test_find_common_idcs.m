% Test_find_common_idcs.m -- verify misc/find_common_idcs.m: both the
% general ('unique'-based) path and the 'nrange' path, against a
% brute-force reference, for shape as well as content.
% Part 1: idcs1/idcs2 are columns also when one input has exactly 1 element
% (find(S1*S2) is row-shaped then). Part 2: random inputs, membership and
% completeness against brute force.
%
% DJ (CC), 10/15/2026: Initial coding;

% ===== PART 1: shape contract for singleton inputs =====
shapes = {[3,1],[1,3]};   % column, row
for i1 = 1:numel(shapes)
    arr1 = reshape([2,5,7],shapes{i1});
    arr2 = 5;

    [idcs1,idcs2] = find_common_idcs(arr1,arr2);
    if ~iscolumn(idcs1) || ~iscolumn(idcs2)
        error("Test_find_common_idcs failed: general path returned non-column idcs for singleton arr2 (shape arr1=%s).",mat2str(shapes{i1}))
    elseif numel(idcs1)~=1 || arr1(idcs1)~=arr2(idcs2)
        error("Test_find_common_idcs failed: general path found wrong match for singleton arr2.")
    end

    [idcs1,idcs2] = find_common_idcs(arr2,arr1);
    if ~iscolumn(idcs1) || ~iscolumn(idcs2)
        error("Test_find_common_idcs failed: general path returned non-column idcs for singleton arr1 (shape arr2=%s).",mat2str(shapes{i1}))
    elseif numel(idcs1)~=1 || arr2(idcs1)~=arr1(idcs2)
        error("Test_find_common_idcs failed: general path found wrong match for singleton arr1.")
    end

    % nrange path: values must be positive integers in [1,nrange]
    nrange = 7;
    [idcs1,idcs2] = find_common_idcs(arr1,arr2,nrange);
    if ~iscolumn(idcs1) || ~iscolumn(idcs2)
        error("Test_find_common_idcs failed: nrange path returned non-column idcs for singleton arr2.")
    elseif numel(idcs1)~=1 || arr1(idcs1)~=arr2(idcs2)
        error("Test_find_common_idcs failed: nrange path found wrong match for singleton arr2.")
    end

    [idcs1,idcs2] = find_common_idcs(arr2,arr1,nrange);
    if ~iscolumn(idcs1) || ~iscolumn(idcs2)
        error("Test_find_common_idcs failed: nrange path returned non-column idcs for singleton arr1.")
    elseif numel(idcs1)~=1 || arr2(idcs1)~=arr1(idcs2)
        error("Test_find_common_idcs failed: nrange path found wrong match for singleton arr1.")
    end
end
% No match: must return empty (0x1) columns, not error
arr1 = [2;4;6];   arr2 = 9;
[idcs1,idcs2] = find_common_idcs(arr1,arr2);
if ~iscolumn(idcs1) || ~iscolumn(idcs2) || ~isempty(idcs1)
    error("Test_find_common_idcs failed: no-match singleton case did not return empty columns (general path).")
end
[idcs1,idcs2] = find_common_idcs(arr1,arr2,9);
if ~iscolumn(idcs1) || ~iscolumn(idcs2) || ~isempty(idcs1)
    error("Test_find_common_idcs failed: no-match singleton case did not return empty columns (nrange path).")
end
disp('Part 1 (singleton-input shape contract) passed.')

% ===== PART 2: random-sized correctness, both paths, vs brute force =====
tol_trials = 30;
for trial = 1:tol_trials
    rng(160000+trial)
    m = randi(10);   n = randi(10);
    nrange = randi(5);                      % small range forces repeats
    arr1 = randi(nrange,[m,1]);
    arr2 = randi(nrange,[n,1]);
    % Either orientation must be accepted
    if rand()<0.5, arr1 = arr1'; end
    if rand()<0.5, arr2 = arr2'; end

    % Brute-force reference: every (i,j) with arr1(i)==arr2(j). Index only
    % through plain column vectors: A(mask) with a non-vector A (II/JJ) and
    % a differently shaped mask silently returns wrong shape/content.
    [II,JJ] = ndgrid(1:m,1:n);
    ii = II(:);   jj = JJ(:);
    a1 = arr1(:); a2 = arr2(:);
    match = a1(ii)==a2(jj);
    ref_pairs = sortrows([ii(match),jj(match)]);

    [idcs1,idcs2] = find_common_idcs(arr1,arr2);
    if ~iscolumn(idcs1) || ~iscolumn(idcs2)
        error("Test_find_common_idcs failed (trial %d): general path returned non-column idcs.",trial)
    elseif numel(idcs1)~=numel(idcs2) || any(a1(idcs1(:))~=a2(idcs2(:)))
        error("Test_find_common_idcs failed (trial %d): general path returned a non-matching pair.",trial)
    end
    got_pairs = sortrows([idcs1(:),idcs2(:)]);
    if ~isequal(got_pairs,ref_pairs)
        error("Test_find_common_idcs failed (trial %d): general path is not complete (missed or extra pairs vs brute force).",trial)
    end

    [idcs1,idcs2] = find_common_idcs(arr1,arr2,nrange);
    if ~iscolumn(idcs1) || ~iscolumn(idcs2)
        error("Test_find_common_idcs failed (trial %d): nrange path returned non-column idcs.",trial)
    end
    got_pairs = sortrows([idcs1(:),idcs2(:)]);
    if ~isequal(got_pairs,ref_pairs)
        error("Test_find_common_idcs failed (trial %d): nrange path is not complete (missed or extra pairs vs brute force).",trial)
    end
end
fprintf('Part 2 (general random-sized correctness) passed: %d/%d trials checked against brute force.\n',tol_trials,tol_trials);
