function E = build_exponent_grid(caps,jointcap,subsetcaps)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% E = BUILD_EXPONENT_GRID(CAPS,JOINTCAP,SUBSETCAPS) returns the matrix of
% exponents of all monomials in numel(CAPS) variables whose degree in
% variable i is at most CAPS(i), whose total degree is at most JOINTCAP, and
% whose degree summed over each variable subset is at most the corresponding
% entry of SUBSETCAPS.
%
% INPUT
% - caps:       'double', 1 x n, nonnegative integers: the largest degree
%               of each variable; caps(i) = 0 excludes variable i;
% - jointcap:   (optional) 'double', scalar, the largest total degree
%               (default sum(caps): no further restriction);
% - subsetcaps: (optional) 'double', 2^n nonnegative integers, one bound
%               per subset of the variables: the bound for the subset whose
%               members are the variables b sits at linear index
%               1+sum(2.^(b-1)), so a singleton i is at 1+2^(i-1) and the
%               full set at 2^n; entry 1, the empty subset, is not used.
%               The bounds apply in addition to CAPS and JOINTCAP (default
%               []: none).
%
% OUTPUT
% - E:          'double', T x n, nonnegative integers, rows in ascending
%               order (sortrows): row t the degree of each variable in
%               monomial t; zeros(1,0), the constant monomial, for n = 0.
%
% NOTES
% CAPS or SUBSETCAPS with negative or non-integer entries, and SUBSETCAPS
% with other than 2^n entries, are errors.
%
% For support, contact M. Peet, Arizona State University at mpeet@asu.edu

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PIETOOLS - build_exponent_grid
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
% Initial coding MMP, 08/22/2026

% Caps as a row of nonnegative integers
caps = reshape(caps,1,[]);
n = numel(caps);

if any(caps<0) || any(caps~=round(caps))
    error("Monomial degrees must be nonnegative integers.")
end

% Default total cap: sum(caps)
if nargin<2 || isempty(jointcap)
    jointcap = sum(caps);
end

% No variables: the only monomial is the constant 1.
if n==0
    E = zeros(1,0);
    return
end

% Tensor grid of all admissible per-variable degrees.
grids = cell(1,n);
for i=1:n
    grids{i} = (0:caps(i))';
end
subs = cell(1,n);
[subs{:}] = ndgrid(grids{:});

% One row per grid point, column i the degree of variable i
E = zeros(numel(subs{1}),n);
for i=1:n
    E(:,i) = subs{i}(:);
end

% Impose the joint degree bound.
E = E(sum(E,2)<=jointcap,:);

% Impose one bound per variable subset, when given. Singletons and the full
% set duplicate 'caps' and 'jointcap'; applying them again is harmless and
% keeps the caller free to pass a cap array verbatim.
if nargin>=3 && ~isempty(subsetcaps)
    if numel(subsetcaps)~=2^n
        error("Subset degrees should be given as a 2^n array for n "...
              +"variables.")
    end
    sc = subsetcaps(:);
    if any(sc<0) || any(sc~=round(sc))
        error("Subset degrees must be nonnegative integers.")
    end
    % Entry i bounds the summed degree of the variables of the set bits of
    % i-1
    for i = 2:2^n
        b = bitget(i-1,1:n)>0;
        E = E(sum(E(:,b),2)<=sc(i),:);
    end
end

% Rows in ascending order
E = sortrows(E);

end
