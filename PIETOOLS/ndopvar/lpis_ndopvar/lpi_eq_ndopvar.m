function prog = lpi_eq_ndopvar(prog,P,opts)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PROG = LPI_EQ_NDOPVAR(PROG,P) takes an LPI optimization program structure
% 'prog' and an 'ndopvar' decision variable P, and adds equality constraints
% enforcing P==0. It is the 'ndopvar' counterpart of 'lpi_eq'.
%
% An 'ndopvar' object represents the operator
%
%   (P(r) x)(s) = sum_{j in {0,1,2}^N} int_{[a,b]} I_j(s,t) R_j(s,t;r) x(t) dt
%
%   R_j(s,t;r) = (Im o Zd(s))^T (Ik o [1;r])^T C_j (In o Zd(t)),
%
% for k = m*prod(deg+1). Since the monomials Zd(s) o Zd(t) are linearly
% independent and the 3^N indicator functions I_j act on disjoint regions,
% the operator vanishes if and only if every entry of every coefficient
% matrix C_j vanishes. Each entry is affine in the decision variables r, so
% P==0 amounts to one linear equality constraint per entry.
%
% Those constraints are handed to 'soseq' by wrapping each C_j in a 'dpvar'
% object. This needs no conversion work: the rows of C_j are grouped as
% (Ik o [1;r]), i.e. one block of q+1 consecutive rows per slot with the
% constant term first, which is exactly the row layout of a 'dpvar'
% coefficient matrix. Folding the monomial indices into the matrix
% dimensions gives a variable-free 'dpvar' whose entries are precisely the
% affine expressions to be set to zero, so no dummy variable substitution is
% needed and the constraint does not depend on 'prog.vartable'.
%
% INPUT
% - prog:   'struct' specifying an LPI program structure to modify (see also
%           'lpiprogram'). All decision variables of P must already appear in
%           'prog.decvartable';
% - P:      'ndopvar' object representing a PI operator decision variable, of
%           which to enforce P==0. Coefficient matrices P.C{i} that are empty
%           or identically zero are taken to be zero and generate no
%           constraints;
% - opts:   (optional) set opts = 'symmetric' if P is known to be
%           self-adjoint, to enforce only one coefficient matrix per adjoint
%           pair to be zero (the other then follows by self-adjointness).
%           Taking the adjoint maps a lower integral to an upper integral in
%           every direction, so the parameters pair up under the index map
%           0->0, 1->2, 2->1, and it suffices to constrain one member of each
%           pair;
%
% OUTPUT
% - prog:   Same LPI program structure as the input, but with the added
%           constraints enforcing P==0.
%
% NOTES
% To enforce P==Q for two 'ndopvar' objects, call LPI_EQ_NDOPVAR(PROG,P-Q).
% The 'minus' method aligns the monomial degrees and the decision variable
% lists of the two operators, so no separate basis handling is needed here.
%
% This routine does not accept 'nopvar' objects: those carry no decision
% variables, so P==0 is a fact to be checked with 'eq', not a constraint to
% be imposed.
%
% See also LPI_EQ, LPI_EQ_2D, SOSEQ, NDOPVAR.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PIETOOLS - lpi_eq_ndopvar
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
% MMP, 08/28/2026: Initial coding


% % % Check the inputs
if isa(P,'polynomial') || isa(P,'double')
    error('Enforcing equality constraints on fixed values or polynomials is not supported.')
elseif isa(P,'nopvar')
    error("Input of type 'nopvar' carries no decision variables; use 'eq' to test whether it is zero.")
elseif ~isa(P,'ndopvar')
    error("Input must be of type 'ndopvar'.")
end
if ~isa(prog,'struct') || ~isfield(prog,'decvartable')
    error("First input must be an LPI program structure; see 'lpiprogram'.")
end

% % % Extract the dimensions of the operator
N = numel(P.deg);
q = numel(P.dvarname);
dvars = cellstr(string(P.dvarname(:)));
nZ = prod(P.deg+1);

% An operator all of whose parameters are empty or zero is already zero, so
% it imposes no constraints. This is checked before the dimensions are read,
% since 'dim' is inferred from the coefficient matrices and is NaN when none
% of them are specified.
is_zero = true;
for ii=1:numel(P.C)
    if ~isempty(P.C{ii}) && nnz(P.C{ii})>0
        is_zero = false;
        break
    end
end
if is_zero
    return
end

dims = P.dim;
if any(isnan(dims))
    error("Dimensions of the operator could not be determined from its coefficient matrices.")
end
m = dims(1);    n = dims(2);

% Every decision variable of the operator must be known to the program,
% otherwise the constraint would silently refer to a variable that does not
% exist in the optimization problem.
if q>0
    known = cellstr(string(prog.decvartable(:)));
    is_new = ~ismember(dvars,known);
    if any(is_new)
        error("Decision variable '"+string(dvars{find(is_new,1)})+"' does not appear "...
              +"in the program; declare it with 'lpidecvar' before imposing constraints.")
    end
end

% % % Determine which coefficient matrices to constrain
sz_C = [size(P.C),1];
n_C = numel(P.C);
if nargin>=3 && ~isempty(opts)
    if ~(ischar(opts) || isstring(opts)) || ~strcmpi(opts,'symmetric')
        error("Third input is not recognized; the only supported option is 'symmetric'.")
    end
    use_symmetry = true;
else
    use_symmetry = false;
end

% % % Impose the constraints, one coefficient matrix at a time
for ii=1:n_C
    Cii = P.C{ii};
    % An empty or zero parameter is already zero, so it needs no constraint.
    if isempty(Cii) || nnz(Cii)==0
        continue
    end

    % Determine the index of element ii along each dimension of the cell C,
    % and hence which directions carry an integral.
    idcs = cell(1,N);
    [idcs{:}] = ind2sub(sz_C,ii);
    idcs = cell2mat(idcs);
    is_int = logical(idcs-1);

    % Under self-adjointness a lower integral pairs with an upper integral,
    % so only one member of each pair has to be constrained.
    if use_symmetry
        adj_idcs = idcs;
        adj_idcs(idcs==2) = 3;
        adj_idcs(idcs==3) = 2;
        adj_cell = num2cell(adj_idcs);
        adj_ii = sub2ind(sz_C,adj_cell{:});
        if adj_ii<ii
            continue
        end
    end

    % A multiplier direction identifies the dummy variable with the primary
    % one, so it contributes no monomials in the dummy variable.
    nZ_t = prod(P.deg(is_int)+1);
    if ~isequal(size(Cii),[m*nZ*(q+1),n*nZ_t])
        error("Coefficient matrix P.C{"+num2str(ii)+"} has size "+mat2str(size(Cii))...
              +", expected "+mat2str([m*nZ*(q+1),n*nZ_t])+".")
    end

    % Wrap the coefficient matrix as a variable-free 'dpvar' whose entries
    % are the affine expressions to be set to zero, and impose the equality.
    Dii = dpvar(sparse(Cii),zeros(1,0),{},dvars,[m*nZ,n*nZ_t]);
    prog = soseq(prog,Dii);
end

end
