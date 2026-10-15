function intmat_i = rand_valid_intmat(d,nterms)
% INTMAT_I = RAND_VALID_INTMAT(D,NTERMS) generates NTERMS independently
% random, individually valid nested integration structures for a single
% spatial direction of degree D, in the 'intop' intmat row format.
%
% INPUTS
% - d:          scalar integer, degree (number of dummy variables) along
%               this direction;
% - nterms:     scalar integer, number of (independently random) terms to
%               generate;
%
% OUTPUTS
% - intmat_i:   nterms x 2*d array, suitable as one element of the
%               'intmat' cell of an 'intop' object (one row per term),
%               always passing the per-row checks of @intop/isvalid.m;
%
% NOTES
% - Each term is a "chain": a uniformly random permutation of 1:d giving
%   the order innermost-to-outermost, with the innermost variable bounded
%   by the next one out and the outermost bounded by the full domain
%   [0,d+1]. A chain is a single unambiguous nested order, so no
%   normalization by split_ints.m (as in rand_intop.m) is needed.

% DJ (CC), 10/15/2026: Initial coding;

intmat_i = zeros(nterms,2*d);
for t=1:nterms
    perm = randperm(d);   % perm(1)=innermost ... perm(d)=outermost
    row = zeros(1,2*d);
    row(2*(perm(d)-1)+1:2*perm(d)) = [0,d+1];        % outermost: full domain
    for k=1:d-1
        row(2*(perm(k)-1)+1:2*perm(k)) = [0,perm(k+1)];  % bounded by next one out
    end
    intmat_i(t,:) = row;
end

end
