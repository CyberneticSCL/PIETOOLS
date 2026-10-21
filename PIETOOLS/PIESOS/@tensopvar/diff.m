function dRop = diff(Rop,vars,d)
% DROP = DIFF(ROP,VARS,D) takes the composition of the tensor-PI operator
% ROP with a Dth order differential operator in variables VARS
%   DROP = d^D/(d vars^D) o ROP
% VARS and D can be row vectors to simultaneously differentiate with
% respect to multiple variables. D can also be a column vector or matrix to
% take multiple derivatives, concatenated vertically.
%
% INPUTS
% - Rop:    m x n 'tensopvar' object representing a tensor-PI operator;
% - vars:   1 x r 'cellstr' or 'polynomial' (pvar) array specifying the
%           names of the spatial variables with respect to which to take
%           the derivative;
% - d:      (optional) p x r array of nonnegative integers specifying the
%           order of the derivatives to be taken with respect to each of
%           the spatial variables. Defaults to 1;
%
% OUTPUTS
% - dRop:   m*p x n 'tensopvar' representing the composition of the
%           requested differential operator with the tensor-PI operator
%           Rop;
%
% NOTES
% - For derivatives of higher order, or in multiple variables, the general
%   Leibniz rule is used: d^D(f1*...*fJ) is the sum over all splits a of
%   the derivatives over the factors (a(j,k) the order of factor j in the
%   kth variable), with the multinomial coefficient, for each variable,
%       D!/(a(1,k)!...a(J,k)!)
%   times the product of the derivatives d^a(j,:) fj. Terms in which a
%   factor that does not depend on the variable is differentiated are
%   zero, and are dropped;
% - The factors are differentiated with 'diff' for 'sopvar', so a factor
%   with a multiplier part in a differentiated variable throws an error;
% - A nested factor (a 'tensopvar' object as factor of Rop) is
%   differentiated recursively, by the same Leibniz rule applied to its own
%   factors. Its derivative is again a nested factor, which can have more
%   terms than the factor itself;
% - Type-2 tensor products of the outputs of multiple factors are not
%   supported, as the outputs depend on different copies of the variable.
%   This holds at every level of nesting;
%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PIETOOLS - diff
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
% DJ, 10/20/2026: Initial coding;

% Check that the operator is properly specified
if ~isa(Rop,'tensopvar')
    error("Tensor-PI operator must be specified as object of type 'tensopvar'.")
end

% Check that the variables and orders of differentiation are properly
% specified
if nargin<=2
    [vars,d] = parse_diff_inputs(vars);
else
    [vars,d] = parse_diff_inputs(vars,d);
end

% Differentiate, applying the Leibniz rule recursively to nested factors
dRop = diff_tensopvar(Rop,vars,d);

end
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %



%% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %%
function [dRop,is_zero] = diff_tensopvar(Rop,vars,d)
% [DROP,IS_ZERO] = DIFF_TENSOPVAR(ROP,VARS,D) returns the derivative of the
% tensor-PI operator ROP, for VARS a 1 x r cellstr of unique names and D a
% p x r array of nonnegative integers (as returned by parse_diff_inputs).
% IS_ZERO is true if the derivative is known to be zero, and DROP then is
% the zero operator.
if ~Rop.chkval
    Rop = initialize(Rop);
end
% Check that the derivative is supported, for this operator and (through the
% recursion) for all nested operators
if Rop.otype(1)==2 && size(Rop.ops,2)~=1
    error("Differentiation of type-2 tensor products of PI operators is currently not supported.")
end
if ~all(cellfun(@(P) isa(P,'sopvar') || isa(P,'tensopvar'),Rop.ops(:)))
    error("Differentiation of tensor-PI operators with factors other than 'sopvar' or 'tensopvar' objects is not supported.")
end
% Return quickly if no derivative is taken: no variables, or all orders zero
is_zero = false;
p = size(d,1);
if all(d(:)==0)
    dRop = repmat(Rop,p,1);
    return
end

% % We compute the derivative using general Leibniz rule:
% %   d^n/dx^n [f(x)g(x)] = 
% %   sum_{k=0}^{n} (n choose k) (d^{n-k} f/dx^{n-k}) (d^{k} g/dx^{k})
% First, establish the terms that appear in dRop, for each determining what
% derivative is taken of each of the factors in Rop. Collect the orders in
% "dmat", with the rows corresponding to the total number of terms, and
% each column corresponding to a factor/variable pair
[ntrms,nfctrs] = size(Rop.ops);
r = numel(vars);
% Collect the splits of each of the p derivative orders listed in d
dmat = zeros(0,nfctrs*r);
for rr=1:p
    % For each term in dRop (element of splits), determine for each of the
    % factors what derivative is taken with respect to each of the
    % variables 
    splits = {zeros(nfctrs,r)};
    for k=1:r
        % Determine all derivative orders of the factors that sum up to
        % d(rr,k): the splits of d(rr,k) into nfctrs nonnegative integers,
        % from the positions of nfctrs-1 separators among d(rr,k)+nfctrs-1
        % slots (stars and bars)
        if nfctrs==1
            C = d(rr,k);
        else
            bars = nchoosek(1:d(rr,k)+nfctrs-1,nfctrs-1);
            nbars = size(bars,1);
            C = [bars,(d(rr,k)+nfctrs)*ones(nbars,1)] - [zeros(nbars,1),bars] - 1;
        end
        % Add the derivative orders to splits
        new_splits = cell(1,numel(splits)*size(C,1));
        idx = 0;
        for s=1:numel(splits)
            for c=1:size(C,1)
                idx = idx+1;
                new_splits{idx} = splits{s};
                new_splits{idx}(:,k) = C(c,:)';
            end
        end
        splits = new_splits;
    end
    % Concatenate the derivative orders for each row of d
    dmat = [dmat; cell2mat(cellfun(@(a) a(:)',splits(:),'UniformOutput',false))];
end
% Splits shared by different rows of d are only computed once
dmat = unique(dmat,'rows');

% Compute the derivatives listed in dmat
newops = cell(0,nfctrs);
for i=1:ntrms
    for u=1:size(dmat,1)
        d_u = reshape(dmat(u,:),nfctrs,r);
        % Differentiate each factor, for each row of d to which the split 
        % applies
        is_valid = false(p,1);      % set to false if term is zero
        coeff = zeros(p,1);         % coefficient with which to multiply each term
        % blocks{rr,j} is the derivative of factor j for row rr of d. Rows
        % to which the split does not apply keep the original factor, and
        % are zeroed through their coefficient 0
        blocks = repmat(Rop.ops(i,:),p,1);
        for rr=1:p
            % Only take derivative that corresponds with d(rr,:)
            if any(sum(d_u,1)~=d(rr,:))
                continue
            end
            % Take the derivative of order d_u of the factors in Rop
            is_zero_trm = false;
            blk = cell(1,nfctrs);
            for j=1:nfctrs
                % vars and d_u(j,:) are normalized by parse_diff_inputs, and
                % the factors are known to be 'sopvar' or 'tensopvar':
                % skip further checks for a 'sopvar', recurse for nested
                if isa(Rop.ops{i,j},'tensopvar')
                    [blk{j},is_zero_j] = diff_tensopvar(Rop.ops{i,j},vars,d_u(j,:));
                else
                    [blk{j},is_zero_j] = diff(Rop.ops{i,j},vars,d_u(j,:),true);
                end
                if is_zero_j
                    % Factor does not depend on a variable differentiated
                    is_zero_trm = true;
                    break
                end
            end
            if ~is_zero_trm
                % Establish the factor with which to multiply this term: the
                % multinomial coefficient D!/(a(1,k)!...a(J,k)!) of the split
                is_valid(rr) = true;
                coeff(rr) = prod(factorial(d(rr,:)))/prod(factorial(d_u(:)));
                blocks(rr,:) = blk;
            end
        end
        if ~any(is_valid)
            continue
        end
        % Stack the rows in each factor, multiplying the first factor with
        % the coefficient
        row_ops = cell(1,nfctrs);
        for j=1:nfctrs
            if j==1
                for rr=1:p
                    if coeff(rr)~=1
                        blocks{rr,1} = coeff(rr)*blocks{rr,1};
                    end
                end
            end
            if p==1
                row_ops{j} = blocks{1,j};
            elseif isa(Rop.ops{i,j},'tensopvar')
                % Nested factors cannot be concatenated with vertcat (it
                % returns a 'tensopmat'): stack them term by term instead
                row_ops{j} = stack_nested(blocks(:,j),is_valid);
            else
                row_ops{j} = vertcat(blocks{:,j});
            end
        end
        newops = [newops; row_ops];
    end
end
if isempty(newops)
    % All terms vanish: return the zero operator of the appropriate size
    is_zero = true;
    newops = cell(1,nfctrs);
    for j=1:nfctrs
        newops{j} = repmat(Rop.ops{1,j},p,1);
    end
    newops{1} = 0*newops{1};
end
% Build the derivative operator
dRop = Rop;
dRop.ops = newops;
dRop = initialize(dRop);

end
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %



%% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %%
function C = stack_nested(B,is_valid)
% C = STACK_NESTED(B,IS_VALID) returns the vertical concatenation of the
% nested operators B{rr}, rr=1,...,p, with zero rows for rr with
% IS_VALID(rr)=false (vertcat of 'tensopvar' objects returns a 'tensopmat').
% The operators B{rr} can have different terms, so the stack is the sum of
% the operators [0;..;B{rr};..;0], each term of which is the stack of
%   [0;F1] .* [F2;F2] .* ...
% for a pointwise product, so that the term is zero outside of block row rr.
% Differentiation does not change the dimensions and variables of nested
% factors, so all B{rr} have the same number of columns.
p = numel(B);
rr_valid = find(is_valid(:))';
C = B{rr_valid(1)};
ops = cell(0,size(C.ops,2));
% Select the terms of each of the valid rows
for rr=rr_valid
    S = B{rr}.ops;
    otype1 = B{rr}.otype(1);
    for i=1:size(S,1)
        for j=1:size(S,2)
            F = S{i,j};
            if j==1
                % First factor: zero in all rows except row rr
                if isa(F,'tensopvar')
                    is_rr = false(p,1);     is_rr(rr) = true;
                    Brr = cell(p,1);        Brr{rr} = F;
                    F = stack_nested(Brr,is_rr);
                else
                    blks = repmat({0*F},p,1);
                    blks{rr} = F;
                    F = vertcat(blks{:});
                end
            elseif otype1==1
                % Other factors of a pointwise product: repeated in all rows
                F = repmat(F,p,1);
            end
            S{i,j} = F;
        end
    end
    ops = [ops; S];
end
C.ops = ops;
C = initialize(C);
end
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
