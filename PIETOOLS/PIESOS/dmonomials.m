function Z = dmonomials(vartab,degs)
% Z = DMONOMIALS(VARTAB,DEGS) constructs a vector of distributed monomials
% of degrees DEGS in the variables VARTAB.
%
% INPUTS
% - vartab: 'polyopvar' object, specifying the variables in which to
%           construct the monomials. Each row of blocks of vartab.C.ops
%           (each entry of vartab.dims.out) is one (vector-valued) variable
%           of which monomial combinations are taken;
% - degs:   n x 1 array specifying the cumulative degrees of the monomials
%           to be included in the vector, or n x q array specifying the
%           degree of each monomial in each of the q variables;
%
% OUTPUTS
% - Z:      m x 1 'polyopvar' object representing a vector of all
%           monomials of the cumulative degrees specified by DEGS (or of
%           the degrees specified row-wise by DEGS);
%
% NOTES
% - The monomials are tensor products of the variables: for degree 2 in x,
%   x o x, a function of two copies of each variable that x depends on,
%   (x o x)(s1,s2) = kron(x(s1),x(s2)).
% - A monomial has the same degree in all entries of a variable. For a
%   variable y of size n, the monomial of degree d is the single vector
%   y o ... o y of size n^d.
%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PIETOOLS - dmonomials
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
% DJ, 10/23/2026: Add support for polynomial expressions in vartab;

if ~isa(vartab,'polyopvar')
    error("Variables to appear in distributed monomials must be specified as object of type 'polyopvar'.")
end

% Check that the degrees make sense
if any(any(degs<0)) || any(any(round(degs)~=degs))
    error("Monomial degrees must be specified as nonnegative integers");
end

% The variables are the rows of blocks of vartab.C.ops 
% (entries of dims.out)
if ~vartab.chkval
    vartab = initialize(vartab);        % dims.out must be up to date
end
dims_var = reshape(vartab.dims.out,1,[]);

% Build the matrix specifying the degrees of the monomials in each variable
nvars = numel(dims_var);
if size(degs,2)==nvars
    % Assume monomial degrees are already specified for all variables
    degmat = degs;
elseif size(degs,2)==1 || size(degs,1)==1
    % Assume cumulative monomial degrees are specified
    degmat = zeros(0,nvars);
    for ii=1:numel(degs)
        d = degs(ii);
        degmat_ii = sparse(1,nvars);
        for jj = 1:nvars
            nZ = size(degmat_ii,1);
            degmat_ii = sprepmat(degmat_ii,d+1,1);
            for kk = 0:d
                degmat_ii(nZ*kk+(1:nZ),jj) = kk;
            end
            degmat_ii = degmat_ii(sum(degmat_ii,2)<=d,:);   % Discard invalid monomials
        end
        degmat = [degmat;degmat_ii(sum(degmat_ii,2)==d,:)];
    end
else
    error("Degrees must be specified as n x 1 array or n x q array for q variables.")
end

% Convert degmat to state variable indices
dmax = max(sum(degmat,2));
nZ = size(degmat,1);
degmat_cell = mat2cell(degmat,ones(nZ,1),nvars);
Zvals = cellfun(@(a) repelem(1:nvars,a),degmat_cell,'UniformOutput',false);
Zvals = cellfun(@(a) [a,zeros(1,dmax-numel(a))],Zvals,'UniformOutput',false);
Zvals = sparse(cell2mat(Zvals));

% Extract the rows of each variable that appears in a monomial, once: the
% monomials are products of these (n1 x 1, n2 x 1,... polyopvar objects).
% vartab{j,:} keeps only the nonempty terms of block row j: vartab(rows,:) keeps all
% monomials of vartab, which leaves unused monomials in the products
vars_i = cell(1,nvars);
for j=find(any(degmat,1))
    vars_i{j} = vartab{j,:};
end

% Declare the distributed monomial basis
Z = [];
for i=1:size(degmat,1)
    % Take the product of the states specified by Zvals
    if ~any(Zvals(i,:))
        % Constant monomial
        Z = [Z;constant2polyopvar(1)];
        continue
    end
    Zvals_i = Zvals(i,Zvals(i,:)>0);
    Z_i = vars_i{Zvals_i(1)};
    for j=Zvals_i(2:end)
        % Take tensor product of the states to compute the monomials
        Z_i = otimes(Z_i,vars_i{j},2);
    end
    Z = [Z; Z_i];
end

end
