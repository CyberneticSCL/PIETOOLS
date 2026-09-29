function tf = eq_copvar(A,B,tol)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% TF = EQ_COPVAR(A,B,TOL) tests whether two containers represent the same
% operator, for the 'eq' methods of 'copvar' and 'cdopvar'.
%
% INPUTS
% - A, B:   'copvar' or 'cdopvar' objects, single 'sopvar'/'sdopvar' blocks
%           (taken as 1 x 1 containers), or the double 0. Any other number
%           is an error, as it is for the blocks;
% - tol:    (optional) absolute tolerance on coefficients, passed to the
%           block 'eq' methods. Default 1e-14, their default. [] means the
%           default;
%
% OUTPUTS
% - tf:     logical scalar. true iff A and B have the same ORDERED sequence % MMP, 09/29/2026
%           of component spaces on each side (each component a space by     % MMP, 09/29/2026
%           name and domain, in the container's row / column order) and     % MMP, 09/29/2026
%           equal kernels componentwise, to within tol; the two block       % MMP, 09/29/2026
%           partitions may differ, and a block row or column with zero      % MMP, 09/29/2026
%           components has no effect. For decision containers, "equal" is   % MMP, 09/29/2026
%           equality as affine functions of the decision variables,         % MMP, 09/29/2026
%           matched by name.                                                % MMP, 09/29/2026
%
% NOTES
% A == 0 is the zero test the executives use ('~(Twop==0)'): every block is
% zero, [] blocks trivially.
%
% Spaces are compared by variable NAME, row by row and column by column,
% with the domains of the variables used; so two containers over merged and
% unmerged registries compare equal when they are the same operator.
% Operators between different spaces are unequal, not an error, following
% the blocks.
%
% One operator can be stored in different block PARTITIONS: blkdiag(Rc,Pc)  % MMP, 09/29/2026
% keeps R^1 x R^1 as two block rows where opvar2copvar(blkdiag(Rop,Pop))    % MMP, 09/29/2026
% has one R^2 row. When the grids or the dimensions differ, both            % MMP, 09/29/2026
% containers are cut to the common refinement of their row partitions and   % MMP, 09/29/2026
% of their column partitions (component ranges between the union of the     % MMP, 09/29/2026
% breakpoints cumsum(dim)), each segment's space is compared, and the cut   % MMP, 09/29/2026
% pair goes through the blockwise comparison below ('eq_refined'). A        % MMP, 09/29/2026
% different ordered sequence of component spaces gives false: it is a       % MMP, 09/29/2026
% different operator, on a different ordered product space. The tolerance   % MMP, 09/29/2026
% means the same on both paths, since the block 'eq' tests the largest      % MMP, 09/29/2026
% coefficient of the difference, and the largest over the sub-blocks is     % MMP, 09/29/2026
% the largest over the block.                                               % MMP, 09/29/2026
% A block row or column with zero components goes to 'eq_refined' too,      % MMP, 09/29/2026
% where it holds no segment: its space then cannot decide the verdict. On   % MMP, 09/29/2026
% the same-partition path it did, so eq was not transitive: [R^0; L2^2]     % MMP, 09/29/2026
% and [L2^0; L2^2] were unequal while each equalled [L2^2] (measured).      % MMP, 09/29/2026
%
% Blocks are compared pairwise: [] against [] is equal; [] against a block
% is that block's zero test; two blocks of one class use that class's own
% 'eq', which aligns monomial and decision bases; a 'sopvar' against an
% 'sdopvar' is compared as (a-b)==0, since @sdopvar/minus accepts a fixed
% operand on either side. Shared by both containers for the same reason as
% 'derive_copvar_meta': it is logic that would silently diverge.
%
% Cost: one block 'eq' per populated pair; each forms the block difference.
% Two decision containers on DIFFERENT lists reconcile them once per block  % MMP, 09/25/2026
% pair ('sync_basis' in @sdopvar/eq), O(q log q) each, not once for the     % MMP, 09/25/2026
% container. Left so: eq is a test-side routine, and the zero test A == 0   % MMP, 09/25/2026
% and containers on one list never reconcile.                               % MMP, 09/25/2026
% Different partitions add one sub-block slice per cut cell where the       % MMP, 09/29/2026
% segment is a proper part of a block: for 'sdopvar', a gather of the kept  % MMP, 09/29/2026
% columns of each B, O(nnz) of those columns, so linear in q.               % MMP, 09/29/2026
% This is a plain function, so every A.name read goes through the           % MMP, 09/29/2026
% overloaded @copvar/@cdopvar subsref, about 13 us a read (measured): each  % MMP, 09/29/2026
% container is read once, into a struct, independent of q.                  % MMP, 09/29/2026
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
% MMP, 09/25/2026: Documented the per-block decision-list reconciliation
%                  cost of comparing containers on different lists. No code
%                  change.
% MMP, 09/29/2026: A == B for one operator stored in two block partitions
%                  (e.g. blkdiag(Rc,Pc), grid 3x3, against
%                  opvar2copvar(blkdiag(Rop,Pop)), grid 2x2) returned false,
%                  since a grid or dimension mismatch returned at once. Those
%                  two returns now go to 'eq_refined', which compares the
%                  ordered component sequences and the kernels on the common
%                  refinement of the partitions. The logic of the
%                  same-partition path, the space checks and the '== 0'
%                  path is unchanged. A block row or column with zero
%                  components also goes to 'eq_refined', so it has no
%                  effect on the verdict (eq was not transitive otherwise).
%                  Each container is read once into a struct ('mA', 'mB';
%                  the C of the '== 0' test), since each property read in
%                  this plain function costs about 13 us through the
%                  overloaded subsref; verdicts are unchanged. The OUTPUTS
%                  entry for tf was:
% - tf:     logical scalar. true iff A and B map between the same spaces    % MMP, 09/29/2026 (was)
%           with the same dimensions and every pair of blocks is equal to   % MMP, 09/29/2026 (was)
%           within tol. For decision containers, "equal" is equality as     % MMP, 09/29/2026 (was)
%           affine functions of the decision variables, matched by name.    % MMP, 09/29/2026 (was)

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
% Each A.name read in this plain function costs about 13 us through the     % MMP, 09/29/2026
% overloaded @copvar/@cdopvar subsref: read the fields once. 'metadata' is  % MMP, 09/29/2026
% a class method, so its own reads are builtin. Values are unchanged.       % MMP, 09/29/2026
mA = metadata(A);   mA.C = A.C;     mB = metadata(B);   mB.C = B.C;         % MMP, 09/29/2026

% % % Same grid, same spaces by name, same dimensions, same domains.
tf = false;
% A different grid or different dimensions is a different PARTITION, not    % MMP, 09/29/2026
% by itself a different operator: decide on the common refinement. Equal    % MMP, 09/29/2026
% grid and dimensions give identical partitions, so a row or column space   % MMP, 09/29/2026
% mismatch below is a different component sequence and returns false.       % MMP, 09/29/2026
% A zero-count row or column goes there as well: it holds no component, so  % MMP, 09/29/2026
% its space must not decide (else [R^0;L2^2] ~= [L2^0;L2^2], though each    % MMP, 09/29/2026
% equals [L2^2]).                                                           % MMP, 09/29/2026
% if ~isequal(size(A.C),size(B.C)),               return,     end           % MMP, 09/29/2026 (was)
% if ~isequal(A.dim_out(:),B.dim_out(:)) || ~isequal(A.dim_in(:),B.dim_in(:)) % MMP, 09/29/2026 (was)
%   return                                                                  % MMP, 09/29/2026 (was)
dims = [mA.dim_out(:);mA.dim_in(:);mB.dim_out(:);mB.dim_in(:)];             % MMP, 09/29/2026
if ~isequal(size(mA.C),size(mB.C)) || ~isequal(mA.dim_out(:),mB.dim_out(:)) ...
        || ~isequal(mA.dim_in(:),mB.dim_in(:)) || ~all(dims)                % MMP, 09/29/2026
    tf = eq_refined(A,B,mA,mB,tol);                                         % MMP, 09/29/2026
    return                                                                  % MMP, 09/29/2026
end
% for i = 1:size(A.C,1)                                                     % MMP, 09/29/2026 (was)
%     if ~same_space(A,A.space_out(i,:),B,B.space_out(i,:)),  return,  end  % MMP, 09/29/2026 (was)
for i = 1:size(mA.C,1)                                                      % MMP, 09/29/2026
    if ~same_space(mA,mA.space_out(i,:),mB,mB.space_out(i,:)), return, end  % MMP, 09/29/2026
end
% for j = 1:size(A.C,2)                                                     % MMP, 09/29/2026 (was)
%     if ~same_space(A,A.space_in(j,:),B,B.space_in(j,:)),    return,  end  % MMP, 09/29/2026 (was)
for j = 1:size(mA.C,2)                                                      % MMP, 09/29/2026
    if ~same_space(mA,mA.space_in(j,:),mB,mB.space_in(j,:)),   return, end  % MMP, 09/29/2026
end

% % % Blockwise.
% for ii = 1:numel(A.C)                                                     % MMP, 09/29/2026 (was)
%     a = A.C{ii};    b = B.C{ii};                                          % MMP, 09/29/2026 (was)
for ii = 1:numel(mA.C)                                                      % MMP, 09/29/2026
    a = mA.C{ii};    b = mB.C{ii};                                          % MMP, 09/29/2026
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
tf = isnumeric(X) && isscalar(X) && X==0;
end


function P = as_cont(X)
% A single block is a 1 x 1 container; containers pass through.
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
tf = true;
% for ii = 1:numel(P.C)                                                     % MMP, 09/29/2026 (was)
%     if ~isempty(P.C{ii}) && ~eq(P.C{ii},0,tol)                            % MMP, 09/29/2026 (was)
C = P.C;            % one read through the overloaded subsref, as above     % MMP, 09/29/2026
for ii = 1:numel(C)                                                         % MMP, 09/29/2026
    if ~isempty(C{ii}) && ~eq(C{ii},0,tol)                                  % MMP, 09/29/2026
        tf = false;
        return
    end
end
end


function tf = same_space(A,mA,B,mB)
% Same variable names, and each on the same domain in both containers.
% A, B: containers, or structs with their fields vars and dom.              % MMP, 09/29/2026
vA = A.vars(mA);    vB = B.vars(mB);
tf = isequal(sort(vA(:))',sort(vB(:))');
if ~tf,     return,     end
[~,iA] = ismember(vA,A.vars);
[~,iB] = ismember(vA,B.vars);
tf = isequal(A.dom(iA,:),B.dom(iB,:));
end


% % % BEGIN new code, MMP, 09/29/2026 - the three functions below, to the END
% % % marker, compare containers whose block partitions differ (NOTES).

function tf = eq_refined(A,B,mA,mB,tol)
% A == B when the grids or dimensions differ, or a row or column has zero
% components. mA, mB: the fields of A, B read once (main function). Refine
% both sides, check the component space of every segment, then compare the
% cut pair blockwise by calling eq_copvar again. The cut pair has one grid,
% equal dimensions, no zero-count row or column (each segment holds at
% least one component) and (checked here) equal spaces, so that call takes
% the same-partition path and cannot come back here. A side with no
% components in either operand gives a 0 x N (or M x 0) cut grid, with no
% kernel entries: the component spaces of the other side decide.
tf = false;
[ok,rA,rB] = refine(mA.dim_out,mB.dim_out);     if ~ok,     return,     end
[ok,cA,cB] = refine(mA.dim_in,mB.dim_in);       if ~ok,     return,     end
for t = 1:numel(rA.blk)
    if ~same_space(mA,mA.space_out(rA.blk(t),:),mB,mB.space_out(rB.blk(t),:)),  return,  end
end
for u = 1:numel(cA.blk)
    if ~same_space(mA,mA.space_in(cA.blk(u),:),mB,mB.space_in(cB.blk(u),:)),    return,  end
end
tf = eq_copvar(cut(A,mA,rA,cA),cut(B,mB,rB,cB),tol);
end


function [ok,a,b] = refine(dA,dB)
% Common refinement of two partitions dA, dB (component counts per block)
% of one component range; ok = false when the totals differ. Segment t is
% components br(t)+1..br(t+1), br the union of the breakpoints; a.blk(t) is
% the block of A holding it and a.loc{t} its indices inside that block,
% likewise b. A zero-count block holds no segment and drops out.
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
% P on the refined grid: block (t,u) is rows r.loc{t}, columns c.loc{u} of
% block (r.blk(t),c.blk(u)) of P, sliced by @sopvar/@sdopvar subsref (which
% keeps the monomial bases and, for sdopvar, the decision list) only where
% the segment is a proper part of the block; [] stays []. Registry, domains
% and Zd are P's own, so the two-argument constructor may trust them. m: the
% fields of P read once. Slicing needs fully sized block parameters: a
% block with the [] / scalar-0 parameter shorthand, which the constructor
% accepts, errors in @sopvar/subsref (pre-existing).
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
meta = metadata(P);
meta.space_out = m.space_out(r.blk,:);      meta.space_in = m.space_in(c.blk,:);
meta.dim_out = cellfun(@numel,r.loc);       meta.dim_in = cellfun(@numel,c.loc);
R = feval(class(P),C,meta);
end

% % % END new code, MMP, 09/29/2026
