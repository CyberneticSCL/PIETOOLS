function C = polyopvar_times_v2(A,B)                         
% C = polyopvar_times_v2(A,B) returns the 'polyopvar' object C
% representing the elementwise product of two distributed polynomials.
%
% This file follows PIETOOLS/polyopvar/@polyopvar/times.m wherever possible.
% Lines marked [V2-ADDED] or [V2-REMOVED] identify the
% differences wrt times.m.
%
% It is used for p1*g in LocalStability.m, which is the product of
% functional distributed polynomials described by Lemma 5 in the Automatica
% paper. The substantive change is the treatment of zero entries in an
% integration-order row (see the marked block below).
%
% In Lemma 5, tensor products of F-PI operators generate integration-order
% matrices containing zero entries. A zero means that the corresponding
% dummy variable is not present in that term; it is not an index that
% should be passed to old_order. The reordering therefore has to be
% applied only to positive entries of omat.
%
% INPUTS
% - A:     m x n 'polyopvar' object.
% - B:     m x n 'polyopvar' object.
%
% OUTPUTS
% - C:     'polyopvar' object representing A.*B.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PIETOOLS - polyopvar_times_v2                         % [V2-CHANGED]
%
% Copyright (C) 2026 PIETOOLS Team
%
% This program is free software; you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation; either version 2 of the License, or
% (at your option) any later version.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% CR, 09/28/2026: LocalStability copy of polyopvar/times.m.
% CR, 09/28/2026: Preserve zero entries in integration-order matrices while
%                   reordering the positive entries required by Lemma 5.

% Make sure the sizes are appropriate for elementwise multiplication.
sz_A = size(A);     is_dim1_A = sz_A==1;
sz_B = size(B);     is_dim1_B = sz_B==1;
if any(~is_dim1_A & ~is_dim1_B & sz_A~=sz_B)
    error("Matrix dimensions must match for elementwise multiplication.")
end

% Account for the case where one object is not a distributed polynomial.
if ~isa(A,'polyopvar')
    C = B;
    C.C = A.*B.C;
    return
elseif ~isa(B,'polyopvar')
    C = A;
    C.C = A.C.*B;
    return
end

% Express A and B in terms of the same variables.
[A,B] = common_vars(A,B);
C = A;

% Declare the full list of monomials in the product.
degmat_full = repelem(A.degmat, size(B.degmat,1), 1) + ...
              repmat(B.degmat,size(A.degmat,1), 1);

% Declare the coefficient operators acting on the full monomial list.
params_full = otimes(A.C,B.C);

% Reorder factors to account for the order of state variables in each
% combined monomial.
d1 = size(A.degmat,1);
d2 = size(B.degmat,1);
nvars = numel(A.varname);
for k = 1:size(degmat_full,1)
    if sum(degmat_full(k,:))<=1
        continue
    end

    [k2,k1] = ind2sub([d2,d1],k);
    deg1 = A.degmat(k1,:);
    deg2 = B.degmat(k2,:);
    Ck = params_full.ops{k};

    % Determine the order in which the state factors occur before sorting.
    state_idcs = [repelem(1:nvars,deg1),repelem(1:nvars,deg2)];
    [~,var_order] = sort(state_idcs);

    if isa(Ck,'tensopvar')
        % For tensor-valued coefficients, only the factor order changes.
        Ck.ops = Ck.ops(:,var_order);
    elseif isa(Ck,'intop')
        % For a functional coefficient, update the names of the dummy
        % variables and the order encoded in its integration-order matrix.
        Ck.pvarname = Ck.pvarname(var_order);
        [~,old_order] = sort(var_order);

        zero_idcs = all(Ck.omat==0,2);
        omat1 = Ck.omat(~zero_idcs,:);

        % [V2-REMOVED] Original times.m applies this line to all entries,
        % including zero placeholders:
        % omat1 = old_order(omat1);
        %
        % [V2-ADDED] Lemma 5 produces rows with zero placeholders. Keep
        % those zeros unchanged and reorder only valid positive indices.
        is_present = omat1>0;
        omat1(is_present) = old_order(omat1(is_present));
        Ck.omat(~zero_idcs,:) = omat1;
    end
    params_full.ops{k} = Ck;
end

% Finally, combine terms involving the same monomial.
[Pmat,degmat_new] = uniquerows_integerTable(degmat_full);
C.degmat = degmat_new;
C.C = params_full;
C.C.ops = cell(size(params_full.ops,1),size(degmat_new,1));
C.C.depmat2 = zeros(size(degmat_new,1),nvars);
for i=1:size(params_full.ops,1)
    for j=1:size(degmat_new,1)
        param_idcs = find(Pmat(:,j));
        C.C.depmat2(j,:) = params_full.depmat2(param_idcs(1),:);
        Ctmp = 0;
        for k=1:numel(param_idcs)
            Ctmp = Ctmp + params_full.ops{i,param_idcs(k)};
        end
        C.C.ops{i,j} = Ctmp;
    end
end

end
