% Test_remove_badints_fuzz.m -- Fuzz-tests the cycle-detection logic of
% remove_badints.m against the independent brute-force ground truth in
% is_valid_intmat_bruteforce.m. remove_badints.m and the internal
% is_consistent() of isvalid.m use the same "safe-set closure" pattern, so
% comparing them to each other would not catch a bug shared by both; this
% test checks against a structurally different (permutation-enumeration)
% method instead.
%
% Rows containing a delta/multiplier sub-term (A(v)==B(v) for some v) are
% skipped: remove_badints eliminates and renumbers those variables before
% returning, so the surviving row need not match the input row verbatim,
% and comparing by row-membership would misreport agreement as mismatch.
%
% DJ (CC), 10/15/2026: Initial coding;
%
% NOTES
% - The delta-elimination path of remove_badints is not exercised: rows
%   with a delta/multiplier sub-term are skipped.

ntests = 500;
max_di = 4;
max_nterms = 15;

ntested = 0;
nmismatch = 0;

for t=1:ntests
    di = randi([1,max_di]);
    nterms = randi([1,max_nterms]);
    intmat_i = zeros(nterms,2*di);
    intmat_i(:,1:2:end) = randi([0,di],[nterms,di]);
    intmat_i(:,2:2:end) = randi([1,di+1],[nterms,di]);

    % Skip rows with a delta/multiplier sub-term (see header note)
    has_delta = any(intmat_i(:,1:2:end)==intmat_i(:,2:2:end),2);
    intmat_i = intmat_i(~has_delta,:);
    if isempty(intmat_i)
        continue
    end
    nterms = size(intmat_i,1);

    % Ground truth
    truth = is_valid_intmat_bruteforce(intmat_i);

    % remove_badints verdict: row kept <=> valid (rows are never modified
    % on this delta-free path, so exact row membership is a fair check)
    intmat_after = remove_badints({intmat_i});
    kept_rows = false(nterms,1);
    for j=1:nterms
        kept_rows(j) = ismember(intmat_i(j,:),intmat_after{1},'rows');
    end

    ntested = ntested + nterms;
    mismatch = kept_rows~=truth;
    if any(mismatch)
        nmismatch = nmismatch + nnz(mismatch);
        bad_rows = find(mismatch);
        for r=bad_rows(:)'
            fprintf('Mismatch: di=%d, row=%s, remove_badints kept=%d, bruteforce valid=%d\n',...
                di,mat2str(intmat_i(r,:)),kept_rows(r),truth(r));
        end
    end
end

fprintf('\n%d individual rows compared; %d mismatches found between remove_badints and brute-force ground truth.\n',ntested,nmismatch);
if nmismatch==0
    fprintf('remove_badints agrees with the independent brute-force check on all tested rows.\n');
end
