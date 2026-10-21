function E = build_exponent_grid(caps,jointcap,subsetcaps)
% E = BUILD_EXPONENT_GRID(CAPS,JOINTCAP,SUBSETCAPS) returns the matrix of
% exponents of all monomials in numel(CAPS) variables whose degree in
% variable i is at most CAPS(i), whose total degree is at most JOINTCAP, and
% whose degree summed over each variable subset is at most the corresponding
% entry of SUBSETCAPS.
%
% INPUTS
% - caps:       1 x n array of nonnegative integers, the per-variable
%               maximal degrees. caps(i)=0 excludes variable i;
% - jointcap:   (optional) scalar bound on the total degree. Defaults to
%               sum(caps), which imposes no additional restriction;
% - subsetcaps: (optional) 2^n array of nonnegative integers, one bound per
%               subset of the variables. The bound for the subset whose
%               members are the bits set in b sits at linear index
%               1+sum(2.^(b-1)), so a singleton i is at 1+2^(i-1) and the
%               full set at 2^n; entry 1, the empty subset, is unused. This
%               is the convention 'poslpivar_2d/build_monoms' uses, so a
%               cap array can be carried across unchanged apart from the
%               variable ordering. Bounds are ADDITIONAL: 'caps' and
%               'jointcap' still apply, and where they overlap the tighter
%               one governs. Defaults to empty, which restores the
%               per-variable-plus-total behaviour exactly;
%
% OUTPUTS
% - E:          T x n array of nonnegative integers, sorted by rows, with
%               row t giving the degrees of each variable in monomial t
%
% MP, 08/22/2026: Initial coding
% MMP, 09/21/2026: Accept one cap per variable SUBSET. Per-variable plus
%                  one total cap cannot express the bases the 2D settings
%                  files describe -- 'build_monoms' prunes with 2^nvars
%                  subset caps -- so converting such a file could only
%                  produce a covering superset, measured at 1.05x to 1.80x
%                  the basis. With subset caps the conversion is exact.

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

% Impose one bound per variable subset, when given. Singletons and the full  % MMP, 09/21/2026
% set duplicate 'caps' and 'jointcap'; applying them again is harmless and   % MMP, 09/21/2026
% keeps the caller free to pass a cap array verbatim.                        % MMP, 09/21/2026
if nargin>=3 && ~isempty(subsetcaps)                                        % MMP, 09/21/2026
    if numel(subsetcaps)~=2^n                                               % MMP, 09/21/2026
        error("Subset degrees should be given as a 2^n array for n "...
              +"variables.")                                                % MMP, 09/21/2026
    end                                                                     % MMP, 09/21/2026
    sc = subsetcaps(:);                                                     % MMP, 09/21/2026
    if any(sc<0) || any(sc~=round(sc))                                      % MMP, 09/21/2026
        error("Subset degrees must be nonnegative integers.")               % MMP, 09/21/2026
    end                                                                     % MMP, 09/21/2026
    for i = 2:2^n                                                           % MMP, 09/21/2026
        b = bitget(i-1,1:n)>0;                                              % MMP, 09/21/2026
        E = E(sum(E(:,b),2)<=sc(i),:);                                      % MMP, 09/21/2026
    end                                                                     % MMP, 09/21/2026
end                                                                         % MMP, 09/21/2026

E = sortrows(E);

end
