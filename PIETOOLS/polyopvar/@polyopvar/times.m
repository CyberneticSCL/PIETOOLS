function C = times(A,B)
% C = times(A,B) returns the 'polyopvar' object C representing the
% elementwise product of two polyopvar objects.
%
% INPUTS
% - A:     m x n 'polyopvar' object; 
% - B:     m x n 'polyopvar' object;
%
% OUTPUTS
% - C:     'polyopvar' object representing A.*B.
%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PIETOOLS - times
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
% DJ, 03/05/2026: Initial coding
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Make sure the sizes are appropriate for elementwise multiplication
sz_A = size(A);     is_dim1_A = sz_A==1;
sz_B = size(B);     is_dim1_B = sz_B==1;
if any(~is_dim1_A & ~is_dim1_B & sz_A~=sz_B)
    error("Matrix dimensions must match for elementwise multiplication.")
end

% Account for case where one object is not a distributed polynomial
if ~isa(A,'polyopvar')
    % Multiplication of constant with polynomial
    C = B;
    C.C = A.*B.C;
    return
elseif ~isa(B,'polyopvar')
    % Multiplication of polynomial with constant
    C = A;
    C.C = A.C.*B;
    return
end

% Express A and B in terms of the same variables
[A,B] = common_vars(A,B);
C = A;

% Declare the full list of monomials appearing in the product of A and B
degmat_full = repmat(A.degmat, size(B.degmat,1), 1) + repelem(B.degmat,size(A.degmat,1), 1);

% Declare the coefficients acting on the full vector of monomials
params_full = otimes(A.C,B.C);
% Reorder factors to account for order of state variables in monomial
d1 = size(A.degmat,1);
d2 = size(B.degmat,1);
nvars = numel(A.varname);
for k = 1:size(degmat_full,1)
    % For monomials of degree 1, no reordering needs to be done
    if sum(degmat_full(k,:))<=1
        continue
    end
    [k1,k2] = ind2sub([d1,d2],k);
    deg1 = A.degmat(k1,:);
    deg2 = B.degmat(k2,:);
    Ck = params_full.ops{k};

    % Determine in which order the states appear in the old monomials
    state_idcs = [repelem(1:nvars,deg1),repelem(1:nvars,deg2)];
    % Set the order of the states in the combined monomial
    [~,var_order] = sort(state_idcs);
    if isa(Ck,'cell')
        % For cell of factors, we need only reorder the factors
        Ck = Ck(var_order);
    elseif isa(Ck,'intop')
        % For functional operator, we need to update the order of the dummy
        % variables used for integration
        Ck.pvarname = Ck.pvarname(var_order);
        % This also requires updating the order of the integrals
        [~,old_order] = sort(var_order);
        zero_idcs = all(Ck.omat==0,2);
        omat1 = Ck.omat(~zero_idcs,:);
        omat1 = old_order(omat1);
        Ck.omat(~zero_idcs,:) = omat1;
        % Combine terms
    end
    params_full.ops{k} = Ck;
end

% Finally, combine terms involving the same monomial
[Pmat,degmat_new] = uniquerows_integerTable(degmat_full);   % deg_full = Pmat*deg_new
C.degmat = degmat_new;
C.C = tensopvar();
for j=1:size(degmat_new,1)
    % Establish which parameters act on the jth monomial
    param_idcs = find(Pmat(:,j));
    % Add the different operators acting on this monomial
    Ctmp = params_full.ops{param_idcs(1)};
    for k=2:numel(param_idcs)
        if isa(Ctmp,'cell')
            % For terms involvin multiple factors, we store the sum using
            % different rows in a cell
            Ctmp = [Ctmp; params_full.ops{param_idcs(k)}];
        else
            Ctmp = Ctmp + params_full.ops{param_idcs(k)};
        end
    end
    % Set the operator acting on the jth monomial
    C.C.ops{j} = Ctmp;
end

end