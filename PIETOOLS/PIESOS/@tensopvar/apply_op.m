function Y = apply_op(Cop,X)
% Y = APPLY_OP(COP,X) applies the tensor-PI operator COP to the
% distributed monomial X to return the function Y.
%
% INPUTS
% - Cop:    m x n 'tensopvar' object mapping
%               Cop: L2^{n1}[Omega1] x ... x L2^{nd}[Omegad]
%                       --> L2^{m1}[Omega1] x ... x L2^{mp}[Omegap];
% - X:      1 x d cell of 'polynomial' objects (or a single 'polynomial'
%           if d=1), representing a distributed monomial
%            X := x1 x... xd in L2^{n1}[Omega1] x ... x L2^{nd}[Omegad];
%           if Cop has nested factors (a 'tensopvar' acting on some of the
%           x_j), x_j is a cell of the inputs of these nested factors, or X
%           is a cell with one input for each of the factors of Cop,
%           counting those of nested operators (num_factors(Cop)), in
%           order;
%
% OUTPUTS
% - Y:      m x r 'polynomial' array, representing the output of Cop upon
%           acting on X;
%
% NOTES
% - The polynomials in X must be defined in terms of the same variables as
%   kernels defining Cop.
% - We must have Cop.otype(2)==2, meaning a proper tensor product is taken
%   of the inputs. If Cop.otype(1)==1, then m = m1 =...= mp, and the output
%   represents the elementwise product of each of the operators applied to
%   each of the factors. If Cop.otype(1)==2, then m = m1*...*mp.
% - Each operator in Cop.ops may be a 'tensopvar' or 'intop' object, to
%   which its own 'apply_op' is applied, or a 'sopvar' object, to which
%   'apply_sopvar' is applied.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PIETOOLS - apply_op (tensopvar)
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
% DJ, 10/15/2026: Initial coding;


% Check that the tensor-PI operator is properly specified
if ~isa(Cop,'tensopvar')
    error("First argument must be of type 'tensopvar'.")
end
if ~Cop.chkval
    Cop = initialize(Cop);
end
if Cop.otype(2)==1
    error("Pointwise products of the input monomial factors are not supported.")
end
ops = Cop.ops;
[ntrms,d] = size(ops);

% Check that the distributed monomial is properly specified
is_flat = false;        % true if X has an input for each factor, also those of nested operators
if ~isa(X,'cell')
    if d==1
        X = {X};
    else
        error("Distributed monomial must be specified as 1xd cell.")
    end
elseif numel(X)~=d
    % Nested operators: X may have one input for each of the num_factors(Cop) factors
    if numel(X)~=num_factors(Cop)
        error("Distributed monomial must be specified as 1xd cell, or with an input for each factor, including those of nested operators.")
    end
    is_flat = true;
end

% Split the inputs over the nested factors, as the 1xd cell with the inputs of the nested factors as cells
Xs = X;
if is_flat
    Xs = split_inputs(X,ops(1,:));
end

% Check that the dimensions match
dims = Cop.dims;
for j=1:d
    xj = Xs{j};
    if isa(xj,'double') || isa(xj,'polynomial') || isa(xj,'dpvar')
        if size(xj,1)~=dims(j,2)
            error("Column dimension of PI operators in tensor-PI must match row dimension of the distributed monomial factors.")
        end
    elseif ~(isa(X,'cell') && isa(ops{1,j},'tensopvar'))
        error("Factors in distributed monomial must be specified as objects of type 'polynomial'.")
    end
end


% Apply each of the factors in each term in the tensor-PI operator to the
% associated factor in the distributed monomial
Y = 0;
vars_out = Cop.vars.out;
depmat_out = logical(Cop.depmat.out);
degs = cumsum(depmat_out,1);
for i=1:ntrms
    Yi = 1;
    % The nested factors may differ between terms, so split the inputs for each term
    Xi = Xs;
    if is_flat && i>1
        Xi = split_inputs(X,ops(i,:));
    end
    for j=1:d
        % 'sopvar' has no apply_op; the others dispatch by class
        if isa(ops{i,j},'sopvar')
            Yij = apply_sopvar(ops{i,j},Xi{j});
        else
            Yij = apply_op(ops{i,j},Xi{j});
        end
        if Cop.otype(1)==1
            % Pointwise product
            Yi = Yi.*Yij;
        else
            % Tensor product
            % --> need to replace variables
            vars_ij = vars_out(depmat_out(j,:));
            degs_ij = degs(j,depmat_out(j,:));
            %vars_new = vars_ij;
            for k=1:numel(vars_ij)
                is_varij = ismember(Yij.varname,vars_ij{k});
                if any(is_varij)
                    Yij.varname(is_varij) = {[vars_ij{k},'_',num2str(degs_ij(k))]};
                    %vars_new{k} = [vars_ij{k},'_',num2str(degs_ij(k))];
                end
            end
            Yi = kron(Yi,Yij);
        end
    end
    Y = Y + Yi;
end

end
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %



%% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %%
function Xs = split_inputs(X,fctrs)
% XS = SPLIT_INPUTS(X,FCTRS) splits the cell X of inputs, one for each factor of the operator, including
% those of nested operators, over the 1 x d cell FCTRS of the factors of one term: XS is the 1 x d cell
% with the input of each factor, and, for a nested factor, the 1 x n cell of the n=num_factors inputs
% that belong to it, which its own apply_op splits further if it has nested factors.
d = numel(fctrs);
is_nested = cellfun('isclass',fctrs,'tensopvar');
n = ones(1,d);
for j=find(is_nested)
    n(j) = num_factors(fctrs{j});
end
if sum(n)~=numel(X)
    error("The number of inputs must match the number of factors in each term, including those of nested operators.")
end
last = cumsum(n);
Xs = cell(1,d);
for j=1:d
    if is_nested(j)
        Xs{j} = X(last(j)-n(j)+1:last(j));      % inputs of the nested factors
    else
        Xs{j} = X{last(j)};
    end
end
end
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %