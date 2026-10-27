function Pop = fill_empty_C(Pop)
% POP = FILL_EMPTY_C(POP) returns the 'nopvar' or 'ndopvar' object POP with
% every empty coefficient matrix POP.C{i1,...,iN} replaced by an all-zero
% sparse matrix of the appropriate dimensions.
%
% INPUTS
% - Pop:    'nopvar' or 'ndopvar' object, of which some of the coefficient
%           matrices Pop.C{i1,...,iN} may be empty, in which case the
%           associated parameter is assumed to be zero;
%
% OUTPUTS
% - Pop:    'nopvar' or 'ndopvar' object defining the same PI operator, but
%           with all coefficient matrices explicitly specified as
%               m*prod(deg+1)*(q+1) x n*prod(deg(is_int)+1)
%           sparse arrays;
%
% NOTES
% The dimensions [m,n] of the operator are inferred from the coefficient
% matrices that are specified. If none are specified, or if the specified
% ones are mutually inconsistent, an error is thrown.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PIETOOLS - fill_empty_C
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
% MMP, 08/23/2026: Initial coding

if ~isa(Pop,'nopvar') && ~isa(Pop,'ndopvar')
    error("Input must be a 'nopvar' or 'ndopvar' object.")
end

% Nothing to do if all parameters are already specified
is_mt = cellfun(@isempty,Pop.C);
if ~any(is_mt(:))
    return
end

deg = Pop.deg;
N = numel(deg);
nZ = prod(deg+1);
if isa(Pop,'ndopvar')
    q = numel(Pop.dvarname);
else
    q = 0;
end

% The dimensions are inferred from the parameters that are specified
dim = Pop.dim;
if any(isnan(dim))
    error("Dimensions of the operator cannot be inferred from the specified parameters.")
end

sz_C = [size(Pop.C),1];
idcs = cell(1,N);
for ii=1:numel(Pop.C)
    if ~is_mt(ii)
        continue
    end
    % Establish along which directions element ii represents an integral,
    % as only those directions contribute a monomial basis in the dummy
    % variables
    [idcs{:}] = ind2sub(sz_C,ii);
    is_int = logical(cell2mat(idcs)-1);
    nZ_t = prod(deg(is_int)+1);
    Pop.C{ii} = sparse(dim(1)*nZ*(q+1),dim(2)*nZ_t);
end

end
