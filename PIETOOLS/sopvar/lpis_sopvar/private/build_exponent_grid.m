function E = build_exponent_grid(caps,jointcap)
% E = BUILD_EXPONENT_GRID(CAPS,JOINTCAP) returns the matrix of exponents of
% all monomials in numel(CAPS) variables whose degree in variable i is at
% most CAPS(i) and whose total degree is at most JOINTCAP.
%
% INPUTS
% - caps:       1 x n array of nonnegative integers, the per-variable
%               maximal degrees. caps(i)=0 excludes variable i;
% - jointcap:   (optional) scalar bound on the total degree. Defaults to
%               sum(caps), which imposes no additional restriction;
%
% OUTPUTS
% - E:          T x n array of nonnegative integers, sorted by rows, with
%               row t giving the degrees of each variable in monomial t
%
% MP, 08/22/2026: Initial coding

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

caps = reshape(caps,1,[]);
n = numel(caps);

if any(caps<0) || any(caps~=round(caps))
    error("Monomial degrees must be nonnegative integers.")
end

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

E = zeros(numel(subs{1}),n);
for i=1:n
    E(:,i) = subs{i}(:);
end

% Impose the joint degree bound.
E = E(sum(E,2)<=jointcap,:);

E = sortrows(E);

end
