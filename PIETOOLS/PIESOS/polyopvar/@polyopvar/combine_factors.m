function f = combine_factors(f)
% F = COMBINE_FACTORS(F) combines consecutive factors of the operators
% acting on the monomials of F that are multipliers from R, i.e. that are
% associated with consecutive zeros in F.Z, into a single factor, e.g.
%   F.Z(i,:) = [1,0,0,2]  -->  F.Z(i,:) = [1,0,2].
%
% INPUTS
% - f:      'polyopvar' object representing a distributed polynomial, see
%           'polyopvar' for the factors that are multipliers from R;
%
% OUTPUTS
% - f:      'polyopvar' object representing the same distributed
%           polynomial, but with each run of consecutive factors that are
%           multipliers from R in the operators acting on a monomial
%           replaced by a single factor. The multiplier of this factor is
%           the elementwise product of the multipliers in the run, mapping
%           to the function space of all variables on which any factor in
%           the run depends (the maximum of the dependency arrays of the
%           factors). The zeros in the monomial that correspond to the run
%           are replaced by a single zero;
%
% NOTES
% - Except as in the next note, only 'tensopvar' operators with
%   otype(1)=1 are supported, i.e. for which the outputs of the factors
%   are multiplied pointwise (the product of the multipliers is then the
%   multiplier of the product of the factors). Nested 'tensopvar'
%   factors, and multipliers depending on decision variables, are not
%   supported either. A monomial acted on by such an operator (in any of
%   the rows of F.C) is left as it is, no error is thrown, and the other
%   monomials are still reduced.
% - A monomial of which all factors act on the constant state (no state in
%   F.Z(i,:), all zeros) is not reduced to a single factor, but to a
%   single block of type 'polynomial' in each of the rows of F.C, a known
%   function (see 'tensopmat'), by converting the operators acting on it
%   (tensopvar2polynomial), which also supports a Kronecker product of
%   the factors in different variables. If an operator cannot be
%   converted, the monomial is reduced as above, or left as it is.
% - Monomials with no consecutive multipliers from R, and monomials acted
%   on by functionals ('intop'), are left as they are.
%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PIETOOLS - combine_factors
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

if ~f.chkval
    f = initialize(f);
end

% Determine which factors in each monomial are multipliers from R: the zeros in positions of
% Z up to the number of factors acting on the monomial, see `polyopvar`
[~,nops] = num_factors(f);
nZ = size(f.Z,1);
W = max([size(f.Z,2);nops]);
Zp = full(f.Z);     Zp(:,end+1:W) = 0;
isR = repmat(1:W,nZ,1)<=nops & Zp==0;
cols = find(any(isR(:,1:end-1) & isR(:,2:end),2))';
if isempty(cols)
    return
end

% Combine the consecutive multipliers in the operators acting on each such monomial
ops = f.C.ops;
nops_new = nops;
Zrows = cell(nZ,1);
for i=cols
    n = nops(i);
    % If all factors act on the constant state (no state in Z(i,:)), the operators are multipliers from R, which
    % are represented by a single block of type polynomial, a known function (see tensopmat, tensopvar2polynomial)
    if ~any(Zp(i,:))
        [new_col,is_done] = convert_column(ops(:,i));
        if is_done
            ops(:,i) = new_col;
            nops_new(i) = 1;                % the block takes up one factor, see num_factors
            Zrows{i} = 0;
            continue
        end
    end
    % Factor j starts a new factor, unless it and the factor before it are multipliers from R
    is_new = true(1,n);
    is_new(2:n) = ~(isR(i,2:n) & isR(i,1:n-1));
    grp = cumsum(is_new);                   % new factor that factor j is part of
    ng = grp(end);
    % Combine the factors in all operators acting on the monomial, or in none of them if this is not supported
    % for one, as Z is the same for all: then the monomial is left as it is
    [new_col,is_done] = combine_column(ops(:,i),grp);
    if ~is_done
        continue
    end
    ops(:,i) = new_col;
    nops_new(i) = ng;
    Zrows{i} = Zp(i,is_new);                % zero for each new factor that is a multiplier from R
end
if all(cellfun(@isempty,Zrows))
    return                                  % nothing could be combined
end

% Remove the zeros of Z that were combined, shifting the other elements to the left
lastnz = zeros(nZ,1);
for i=1:nZ
    if ~isempty(Zrows{i})
        Zp(i,:) = 0;
        Zp(i,1:numel(Zrows{i})) = Zrows{i};
    end
    lastnz(i) = max([0,find(Zp(i,:),1,'last')]);
end
Zp = Zp(:,1:max([0;nops_new;lastnz]));

% Set the new polynomial, which initialize checks
Cop = f.C;
Cop.ops = ops;
f.C = Cop;
f.Z = sparse(Zp);
f.chkval = false;
f = initialize(f);

end


function [col,is_done] = convert_column(col)
% [COL,IS_DONE] = CONVERT_COLUMN(COL) converts each of the operators in the cell COL, the blocks acting
% on a monomial of which all factors act on the constant state, to a polynomial (tensopvar2polynomial).
% Blocks that are empty, or already of type polynomial/double/dpvar, are kept. IS_DONE is false, and COL
% empty, if one block cannot be converted (as Z is the same for all blocks): then none is.
is_done = true;
for r=1:numel(col)
    op = col{r};
    if isempty(op) || isa(op,'double') || isa(op,'polynomial') || isa(op,'dpvar')
        continue
    end
    [col{r},is_done] = tensopvar2polynomial(op);
    if ~is_done
        col = [];
        return
    end
end
end


function [col,is_done] = combine_column(col,grp)
% [COL,IS_DONE] = COMBINE_COLUMN(COL,GRP) combines the factors of the tensor-PI operators in the
% cell COL, the blocks acting on a single monomial, such that factor j is combined with the
% other factors with the same value GRP(j), which are consecutive multipliers from R. IS_DONE is
% false, and COL empty, if this is not supported for one of the operators: with otype(1)=2, with
% nested factors, or with multipliers that depend on decision variables.
is_done = true;
ng = grp(end);
for r=1:numel(col)
    op = col{r};
    if ~isa(op,'tensopvar')
        continue
    elseif op.otype(1)~=1
        is_done = false;                    % the outputs of the factors are not multiplied pointwise
        break
    end
    fctrs = op.ops;                         % K x n array of factors
    if size(fctrs,2)~=numel(grp)
        is_done = false;                    % nested factors
        break
    end
    new_fctrs = cell(size(fctrs,1),ng);
    for g=1:ng
        js = find(grp==g);
        for k=1:size(fctrs,1)
            if isscalar(js)
                new_fctrs{k,g} = fctrs{k,js};
            else
                [new_fctrs{k,g},is_done] = multiply_multipliers(fctrs(k,js));
                if ~is_done
                    break
                end
            end
        end
        if ~is_done
            break
        end
    end
    if ~is_done
        break
    end
    op.ops = new_fctrs;
    col{r} = initialize(op);
end
if ~is_done
    col = [];
end
end


function [Pop,is_done] = multiply_multipliers(fctrs)
% [POP,IS_DONE] = MULTIPLY_MULTIPLIERS(FCTRS) takes a 1 x L cell of multipliers from R to L2,
% 'sopvar' objects, and returns the multiplier from R to L2 of the elementwise product of the
% functions that they map the dummy input 1 to, in the variables on which any of them depends.
% IS_DONE is false if a factor is no multiplier from R, or depends on decision variables.
is_done = true;     Pop = [];
vars_out = cell(1,0);       dom_out = zeros(0,2);
p = 1;
for l=1:numel(fctrs)
    fl = fctrs{l};
    if ~isa(fl,'sopvar') || size(fl,2)~=1 || ~isempty(fl.vars.in)
        is_done = false;
        return
    end
    pl = apply_sopvar(fl,1);
    if ~isa(pl,'polynomial')
        is_done = false;
        return
    end
    p = p.*pl;
    vars_out = [vars_out,reshape(fl.vars.out,1,[])];
    dom_out = [dom_out; fl.dom.out];
end
[vars_out,idcs] = unique(vars_out);
vars_out = reshape(vars_out,1,[]);
vars = struct('in',{cell(1,0)},'out',{vars_out});
dom = struct('in',zeros(0,2),'out',dom_out(idcs,:));
Pop = mat2sopvar(p,vars,dom);
end
