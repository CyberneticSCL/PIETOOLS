function ords = enumerate_valid_orders(Llims,Ulims,expnd_vars,var_nums,n)
% ORDS = ENUMERATE_VALID_ORDERS(LLIMS,ULIMS,EXPND_VARS,VAR_NUMS,N) returns
% all possible orderings of the variables VARS(EXPND_VARS) that are 
% consistent with the domains of integration specified by LLIMS and ULIMS.
%
% INPUTS
% - Llims,Ulims:  1 x d arrays, the indices of the lower and upper limits
%                 of integration of each of the d variables (0 and d+1
%                 for the boundaries a and b of the domain);
% - expnd_vars:   1 x d logical, the variables whose orderings are
%                 enumerated;
% - var_nums:     1 x n array, find(expnd_vars)
% - n:            numel(var_nums);
%
% OUTPUTS
% - ords:         nprms x n array as described above; nprms==0 is possible
%                 (an unsatisfiable/contradictory constraint set).
%
% NOTES
%
% DETAILS
% We consider the variables specified by expnd_vars, and integrate them all
% over the full domain, [a,b], represented by [0,di+1].
% We split this full domain into n! subdomains, considering all possible
% orderings of the n variables.
% We keep track of the ordering by assigning a magnitude to each variable,
% 'ord', so that ord(i)>=ord(j) implies si>=sj. For example, if n=3,
%   ord(i,:) = [3,1,2] implies s2 <= s3 <= s1 in term i

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PIETOOLS - enumerate_valid_orders
%
% Copyright (C) 2026 PIETOOLS Team
%
% This program is free software; you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation; either version 2 of the License, or
% (at your option) any later version.
%
% This program is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
%
% You should have received a copy of the GNU General Public License
% along with this program; if not, write to the Free Software
% Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% If you modify this code, document all changes carefully and include date
% authorship, and a brief description of modifications
%
% For support, contact M. Peet, Arizona State University at mpeet@asu.edu
% DJ, 10/22/2026: Initial coding (moved from 'mtimes_tensopvar');

% For each expand variable (by its position 1..n in var_nums), find the
% position of its Llims/Ulims reference WITHIN var_nums, when that
% reference is itself an expand variable (0 otherwise, meaning "no
% constraint" -- the reference is a fixed/non-expand variable or a domain
% endpoint).
d = numel(expnd_vars);
pred = zeros(1,n);      succ = zeros(1,n);
for c = 1:n
    lv = Llims(var_nums(c));
    if lv>=1 && lv<=d && expnd_vars(lv)
        pred(c) = find(var_nums==lv,1);
    end
    uv = Ulims(var_nums(c));
    if uv>=1 && uv<=d && expnd_vars(uv)
        succ(c) = find(var_nums==uv,1);
    end
end

% Build a "must come before" adjacency and each node's indegree, for
% Kahn-style backtracking: a node is "ready" to be placed next once all its
% prerequisites are already placed. pred(c)=p and succ(p)=c can assert the
% IDENTICAL constraint "p before c" (discovered from two different sides:
% c's own Llims, and separately p's own Ulims) -- collect each node's
% prerequisite set and deduplicate with unique() before computing indegree,
% or such a pair double-counts c's indegree and c can never become "ready"
% (confirmed by direct trace: d=2 chain Llims=[0,1] Ulims=[2,3] gives
% pred(2)=1 AND succ(1)=2, the same "1 before 2" constraint from both
% sides).
parents = cell(1,n);
for c = 1:n
    if pred(c)~=0
        parents{c}(end+1) = pred(c);
    end
    if succ(c)~=0
        s = succ(c);
        parents{s}(end+1) = c;
    end
end
children = cell(1,n);
indeg0 = zeros(1,n);
for c = 1:n
    % unique([]) returns a 0x1 (column), not 1x0 -- 'for p = (0x1 empty)'
    % iterates ONCE with p empty (only a 1x0 row iterates zero times), so
    % force row shape or an empty prerequisite set spuriously adds an edge
    % to an empty index (confirmed: this exact issue caused every closure
    % with an unconstrained variable to find zero valid orderings, i.e.
    % every real case -- every variable has at least one free bound).
    parents{c} = reshape(unique(parents{c}),1,[]);
    indeg0(c) = numel(parents{c});
    for p = parents{c}
        children{p}(end+1) = c;
    end
end

rows = cell(0,1);
seq = zeros(1,n);
placed = false(1,n);
indeg = indeg0;
backtrack(1);

if isempty(rows)
    ords = zeros(0,n);
else
    ords = vertcat(rows{:});
end

    function backtrack(pos)
        % Nested function: shares seq/placed/indeg/rows with the parent.
        % Places a variable at position POS of the ordering, for each
        % variable whose predecessors are placed, and recurses; a
        % complete ordering is appended to the rows.
        %
        % INPUT
        % - pos:    integer, the position to fill (1..n+1);
        %
        % OUTPUT
        % - none: the orderings are collected in the parent workspace.
        if pos>n
            row = zeros(1,n);
            row(seq) = 1:n;
            rows{end+1,1} = row;
            return
        end
        avail = find(~placed & indeg==0);
        for v = avail
            seq(pos) = v;
            placed(v) = true;
            ch = children{v};
            indeg(ch) = indeg(ch)-1;
            backtrack(pos+1);
            indeg(ch) = indeg(ch)+1;
            placed(v) = false;
        end
    end

end
