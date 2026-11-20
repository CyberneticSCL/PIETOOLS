function tf = eq_copvar(A,B,tol)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% TF = EQ_COPVAR(A,B,TOL) tests whether two 'copvar' or 'cdopvar'
% containers represent the same operator.
%
% INPUTS
% - A, B:   'copvar' or 'cdopvar' objects, single 'sopvar'/'sdopvar' blocks
%           (taken as 1 x 1 containers), or the double 0. Any other number
%           is an error, as it is for the blocks;
% - tol:    (optional) absolute tolerance on coefficients of the blockwise
%           comparison. Default 1e-14. [] means the default;
%
% OUTPUTS
% - tf:     logical scalar. true iff A and B have the same ORDERED sequence
%           of component spaces on each side (each component a space by
%           name and domain, in the row / column order of the container) and
%           equal kernels componentwise, to within tol; the two block
%           partitions may differ, and a block row or column with zero
%           components has no effect. For decision containers, "equal" is
%           equality as affine functions of the decision variables,
%           matched by name.
%
% NOTES
% A == 0 tests that every block is zero, [] blocks trivially.
%
% Spaces are compared by variable NAME, row by row and column by column,
% with the domains of the variables used; so two containers over merged and
% unmerged registries compare equal when they are the same operator.
% Operators between different spaces are unequal, not an error, as for the
% blocks.
%
% One operator can be stored in different block PARTITIONS, e.g.
% R^1 x R^1 as two block rows or as one R^2 row. When the grids or the
% dimensions differ, both containers are cut to the common refinement of
% their row partitions and of their column partitions (component ranges
% between the union of the breakpoints cumsum(dim)), the space of each
% segment is compared, and the cut pair goes through the blockwise
% comparison below. A different ordered sequence of component spaces gives
% false: it is a different operator, on a different ordered product space.
% The tolerance means the same on both paths, since the block 'eq' tests
% the largest coefficient of the difference, and the largest over the
% sub-blocks is the largest over the block.
% A block row or column with zero components takes the refinement path
% too, where it holds no segment, so its space cannot decide the verdict:
% [R^0; L2^2] and [L2^0; L2^2] both equal [L2^2].
%
% Blocks are compared pairwise: [] against [] is equal; [] against a block
% is the zero test of that block; two blocks of one class use the 'eq' of
% that class, which aligns monomial and decision bases; a 'sopvar' against
% an 'sdopvar' is compared as (a-b)==0.
%
% Cost: one block 'eq' per populated pair; each forms the block difference.
% Two decision containers on DIFFERENT lists reconcile them once per block
% pair, O(q log q) each, not once for the container; the zero test A == 0
% and containers on one list never reconcile.
% Different partitions add one sub-block slice per cut cell where the
% segment is a proper part of a block: for 'sdopvar', a gather of the kept
% columns of each B, O(nnz) of those columns, so linear in q.
% Every A.name read in this plain function goes through the overloaded
% indexing of the container, about 13 us a read (measured): each container
% is read once, into a struct, independent of q.
%
% See also EQ, PLUS, MINUS, COPVAR, CDOPVAR.
%
% For support, contact M. Peet, Arizona State University at mpeet@asu.edu

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PIETOOLS - eq_copvar
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
% Initial coding MMP, 09/25/2026

% Default tolerance
if nargin<3 || isempty(tol)
    tol = 1e-14;
end

% % % The double 0, on either side.
Az = is_numeric_zero(A);        Bz = is_numeric_zero(B);
if (isnumeric(A) && ~Az) || (isnumeric(B) && ~Bz)
    error('eq:badInput',['A container can be compared with another '...
        'container, a sopvar/sdopvar block, or 0; not with a nonzero number.'])
end
if Az && Bz,    tf = true;                  return,     end
if Bz,          tf = is_zero(as_cont(A),tol);   return,     end
if Az,          tf = is_zero(as_cont(B),tol);   return,     end

A = as_cont(A);     B = as_cont(B);
% Each A.name read in this plain function costs about 13 us through the
% overloaded indexing: read the fields once. 'metadata' is a class method,
% so its own reads are builtin.
mA = metadata(A);   mA.C = A.C;     mB = metadata(B);   mB.C = B.C;

% % % Same grid, same spaces by name, same dimensions, same domains.
tf = false;
% A different grid or different dimensions is a different PARTITION, not
% by itself a different operator: decide on the common refinement. Equal
% grid and dimensions give identical partitions, so a row or column space
% mismatch below is a different component sequence and returns false.
% A zero-count row or column goes there as well: it holds no component, so
% its space must not decide (else [R^0;L2^2] ~= [L2^0;L2^2], though each
% equals [L2^2]).
dims = [mA.dim_out(:);mA.dim_in(:);mB.dim_out(:);mB.dim_in(:)];
if ~isequal(size(mA.C),size(mB.C)) || ~isequal(mA.dim_out(:),mB.dim_out(:)) ...
        || ~isequal(mA.dim_in(:),mB.dim_in(:)) || ~all(dims)
    % Refine both sides, check the component space of every segment, then
    % compare the cut pair blockwise by calling eq_copvar again. The cut pair has one grid, equal dimensions, no
    % zero-count row or column (each segment holds at least one component)
    % and (checked here) equal spaces, so that call takes the same-partition
    % path and cannot come back here. A side with no components in either
    % operand gives a 0 x N (or M x 0) cut grid, with no kernel entries: the
    % component spaces of the other side decide.
    tf = false;
    [ok,rA,rB] = refine(mA.dim_out,mB.dim_out);     if ~ok,     return,     end
    [ok,cA,cB] = refine(mA.dim_in,mB.dim_in);       if ~ok,     return,     end
    % The space of every row segment, then of every column segment
    for t = 1:numel(rA.blk)
        if ~same_space(mA,mA.space_out(rA.blk(t),:),mB,mB.space_out(rB.blk(t),:)),  return,  end
    end
    for u = 1:numel(cA.blk)
        if ~same_space(mA,mA.space_in(cA.blk(u),:),mB,mB.space_in(cB.blk(u),:)),    return,  end
    end
    tf = eq_copvar(cut(A,mA,rA,cA),cut(B,mB,rB,cB),tol);
    return
end
% Same partitions: the same space in every row and every column
for i = 1:size(mA.C,1)
    if ~same_space(mA,mA.space_out(i,:),mB,mB.space_out(i,:)), return, end
end
for j = 1:size(mA.C,2)
    if ~same_space(mA,mA.space_in(j,:),mB,mB.space_in(j,:)),   return, end
end

% % % Blockwise.
for ii = 1:numel(mA.C)
    a = mA.C{ii};    b = mB.C{ii};
    % [] against [] is equal, [] against a block is its zero test, two blocks
    % of one class are compared directly, of two classes by their difference
    if isempty(a) && isempty(b)
        continue
    elseif isempty(a)
        ok = eq(b,0,tol);
    elseif isempty(b)
        ok = eq(a,0,tol);
    elseif strcmp(class(a),class(b))
        ok = eq(a,b,tol);
    else
        ok = eq(a-b,0,tol);
    end
    if ~ok,     return,     end
end
tf = true;

end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function tf = is_numeric_zero(X)
% TF = IS_NUMERIC_ZERO(X) is true for the numeric scalar 0.
% INPUT   X: any value.  OUTPUT  tf: logical.
tf = isnumeric(X) && isscalar(X) && X==0;
end


function P = as_cont(X)
% P = AS_CONT(X): a single block is a 1 x 1 container; containers pass
% through; any other class is an error.
% INPUT   X: 'copvar', 'cdopvar', 'sopvar' or 'sdopvar' object.
% OUTPUT  P: 'copvar' or 'cdopvar' object.
if isa(X,'copvar') || isa(X,'cdopvar')
    P = X;
elseif isa(X,'sdopvar')
    P = cdopvar({X});
elseif isa(X,'sopvar')
    P = copvar({X});
else
    error('eq:badInput',['Cannot compare a container with a ''%s''. '...
        'Operands must be copvar/cdopvar, sopvar/sdopvar, or 0.'],class(X))
end
end


function tf = is_zero(P,tol)
% TF = IS_ZERO(P,TOL) is true when every populated block of the container P
% is zero within TOL.
% INPUTS  P: 'copvar' or 'cdopvar' object; tol: absolute tolerance.
% OUTPUT  tf: logical.
tf = true;
C = P.C;            % one read through the overloaded indexing, as above
for ii = 1:numel(C)
    if ~isempty(C{ii}) && ~eq(C{ii},0,tol)
        tf = false;
        return
    end
end
end


function tf = same_space(A,mA,B,mB)
% TF = SAME_SPACE(A,MA,B,MB): the same variable names, each on the same
% domain in both containers.
% INPUTS  A, B: containers, or structs with their fields vars and dom;
%         mA, mB: logical masks of one space over A.vars and B.vars.
% OUTPUT  tf: logical.
vA = A.vars(mA);    vB = B.vars(mB);
tf = isequal(sort(vA(:))',sort(vB(:))');
if ~tf,     return,     end
[~,iA] = ismember(vA,A.vars);
[~,iB] = ismember(vA,B.vars);
tf = isequal(A.dom(iA,:),B.dom(iB,:));
end


function [ok,a,b] = refine(dA,dB)
% Common refinement of two partitions dA, dB (component counts per block)
% of one component range; ok = false when the totals differ. Segment t is
% components br(t)+1..br(t+1), br the union of the breakpoints; a.blk(t) is
% the block of A holding it and a.loc{t} its indices inside that block,
% likewise b. A zero-count block holds no segment and drops out.
% INPUTS  dA, dB: component counts per block of the two partitions.
% OUTPUT  ok: logical; a, b: structs with fields blk (segment -> block) and
%         loc (segment -> indices inside the block).
cA = [0;cumsum(dA(:))];     cB = [0;cumsum(dB(:))];
ok = cA(end)==cB(end);
a = struct('blk',zeros(0,1),'loc',{cell(0,1)});     b = a;
if ~ok,     return,     end
br = unique([cA;cB]);
for t = 1:numel(br)-1
    k = br(t)+1:br(t+1);
    i = sum(cA(2:end)<k(1))+1;      a.blk(t,1) = i;     a.loc{t,1} = k-cA(i);
    j = sum(cB(2:end)<k(1))+1;      b.blk(t,1) = j;     b.loc{t,1} = k-cB(j);
end
end


function R = cut(P,m,r,c)
% R = CUT(P,M,R,C) is P on the refined grid: block (t,u) is rows r.loc{t},
% columns c.loc{u} of block (r.blk(t),c.blk(u)) of P, sliced by block
% indexing (which keeps the monomial bases and, for sdopvar, the decision
% list) only where the segment is a proper part of the block; [] stays [].
% Registry, domains and Zd are those of P, so the two-argument constructor
% may trust them. m: the fields of P read once. Slicing keeps a [] /
% scalar-0 parameter cell as a zero cell.
% INPUTS  P: container; m: its fields; r, c: the refinements of 'refine'.
% OUTPUT  R: container of the class of P on the refined grid.
C = cell(numel(r.blk),numel(c.blk));
for t = 1:numel(r.blk)
    for u = 1:numel(c.blk)
        X = m.C{r.blk(t),c.blk(u)};
        if isempty(X),  continue,   end
        if numel(r.loc{t})<m.dim_out(r.blk(t)) || numel(c.loc{u})<m.dim_in(c.blk(u))
            X = X(r.loc{t},c.loc{u});
        end
        C{t,u} = X;
    end
end
% Metadata of the refined grid: one space and count per segment
meta = metadata(P);
meta.space_out = m.space_out(r.blk,:);      meta.space_in = m.space_in(c.blk,:);
meta.dim_out = cellfun(@numel,r.loc);       meta.dim_in = cellfun(@numel,c.loc);
R = feval(class(P),C,meta);
end
