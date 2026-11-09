function intmat_i = rand_ties_intmat(d,nterms,pdelta)
% INTMAT_I = RAND_TIES_INTMAT(D,NTERMS,PDELTA) generates NTERMS
% independently random, individually valid nested integration structures
% for a single spatial direction of degree D, in the 'intop' intmat row
% format, in which several variables may be TIED in the order of
% integration.
%
% INPUTS
% - d:          scalar integer, degree (number of dummy variables) along
%               this direction;
% - nterms:     scalar integer, number of (independently random) terms to
%               generate;
% - pdelta:     (optional) scalar in [0,1], probability that a variable
%               that is not integrated over the full domain is a
%               multiplier, [q,q], instead of an integral. Defaults to 0;
%
% OUTPUTS
% - intmat_i:   nterms x 2*d array, suitable as one element of the
%               'intmat' cell of an 'intop' object (one row per term),
%               always passing the per-row checks of @intop/isvalid.m;
%
% NOTES
% - Each variable is assigned a random level (at most 3 levels). The
%   variables at level 1 are integrated over the full domain [0,d+1], so
%   that e.g. d=2 with both variables at level 1 gives the full square
%   [0,3,0,3]. A variable at level l>1 is bounded by a random variable of
%   level l-1: [0,q] or [q,d+1], or [q,q] for a multiplier.
% - Variables at the same level are tied: their order of integration is
%   not determined by the term, see int_order in 'intop'. This differs
%   from rand_valid_intmat.m, which only generates chains (no ties).

% DJ (CC), 10/15/2026: Initial coding;

if nargin<3
    pdelta = 0;
end

intmat_i = zeros(nterms,2*d);
for t=1:nterms
    [~,~,lvl] = unique(randi(min(d,3),1,d));      % dense levels 1,...,L
    lvl = lvl(:)';
    row = zeros(1,2*d);
    for j=1:d
        if lvl(j)==1
            row(2*j-1:2*j) = [0,d+1];
        else
            cand = find(lvl==lvl(j)-1);
            q = cand(randi(numel(cand)));
            if rand<pdelta
                row(2*j-1:2*j) = [q,q];
            elseif rand<0.5
                row(2*j-1:2*j) = [0,q];
            else
                row(2*j-1:2*j) = [q,d+1];
            end
        end
    end
    intmat_i(t,:) = row;
end

end
