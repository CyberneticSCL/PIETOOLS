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
% - tf:     logical scalar. true iff A and B map between the same spaces
%           with the same dimensions and every pair of blocks is equal to
%           within tol. For decision containers, "equal" is equality as
%           affine functions of the decision variables, matched by name.
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

% % % Same grid, same spaces by name, same dimensions, same domains.
tf = false;
if ~isequal(size(A.C),size(B.C)),               return,     end
if ~isequal(A.dim_out(:),B.dim_out(:)) || ~isequal(A.dim_in(:),B.dim_in(:))
    return
end
for i = 1:size(A.C,1)
    if ~same_space(A,A.space_out(i,:),B,B.space_out(i,:)),  return,  end
end
for j = 1:size(A.C,2)
    if ~same_space(A,A.space_in(j,:),B,B.space_in(j,:)),    return,  end
end

% % % Blockwise.
for ii = 1:numel(A.C)
    a = A.C{ii};    b = B.C{ii};
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
for ii = 1:numel(P.C)
    if ~isempty(P.C{ii}) && ~eq(P.C{ii},0,tol)
        tf = false;
        return
    end
end
end


function tf = same_space(A,mA,B,mB)
% Same variable names, and each on the same domain in both containers.
vA = A.vars(mA);    vB = B.vars(mB);
tf = isequal(sort(vA(:))',sort(vB(:))');
if ~tf,     return,     end
[~,iA] = ismember(vA,A.vars);
[~,iB] = ismember(vA,B.vars);
tf = isequal(A.dom(iA,:),B.dom(iB,:));
end
