function is_valid = is_valid_intmat_bruteforce(intmat_i)
% IS_VALID = IS_VALID_INTMAT_BRUTEFORCE(INTMAT_I) checks, by brute-force
% enumeration of all integration orders, whether each row of INTMAT_I (in
% the 'intop' intmat format for a single spatial direction) specifies a
% consistent (non-cyclic) nested integral, i.e. whether some order exists
% in which the di variables can be integrated (inner to outer) such that
% the bounds of every variable are the domain boundary (0 or di+1) or a
% variable already integrated at that point. Used as ground truth to
% fuzz-test the cycle detection of remove_badints.m and isvalid.m.
%
% INPUTS
% - intmat_i:   qi x 2*di array, as for the field intmat{i} of an 'intop'
%               object. Rows with a delta/multiplier sub-term (some
%               variable v with A(v)==B(v)) are not supported;
%
% OUTPUTS
% - is_valid:   qi x 1 logical array, true if row j admits at least one
%               valid integration order;
%
% NOTES
% - All di! integration orders are enumerated per row, so only use this
%   function for small di.
% - Rows with a delta/multiplier sub-term are handled differently by
%   remove_badints.m (variable elimination and renumbering), so the fuzz
%   test skips them.
% - A valid order exists iff the dependency graph "v needs A(v),B(v)
%   resolved first" is topologically sortable, which is what
%   is_consistent() in isvalid.m tests (safe-set closure, by BFS layers).
%   Non-degeneracy of the chain of inequalities (e.g. a<=x<=y<=b) is not
%   required.

% DJ (CC), 10/15/2026: Initial coding;

[nterms,d2] = size(intmat_i);
di = d2/2;
is_valid = false(nterms,1);
if di==0
    is_valid = true(nterms,1);
    return
end
perms_di = perms(1:di);        % all di! candidate integration orders

for j=1:nterms
    Aj = intmat_i(j,1:2:end);
    Bj = intmat_i(j,2:2:end);
    for p=1:size(perms_di,1)
        order = perms_di(p,:);     % order(1) integrated first (innermost)
        resolved = [0,di+1];
        ok = true;
        for k=1:di
            v = order(k);
            if ~ismember(Aj(v),resolved) || ~ismember(Bj(v),resolved)
                ok = false;
                break
            end
            resolved = [resolved,v];
        end
        if ok
            is_valid(j) = true;
            break
        end
    end
end

end
