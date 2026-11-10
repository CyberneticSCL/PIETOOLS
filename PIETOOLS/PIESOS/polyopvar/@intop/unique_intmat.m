function Kop = unique_intmat(Kop)
% KOP = UNIQUE_INTMAT(KOP) adds the parameters of the functional KOP that
% are associated with the same interval of integration, so that the rows of
% KOP.intmat{i} are distinct for each i.
%
% INPUTS
% - Kop:    'intop' object representing a functional;
%
% OUTPUTS
% - Kop:    'intop' object representing the same functional as the input,
%           with the terms associated with duplicate integrals combined;
%
% NOTES
% - The parameters are sparse matrices in the same basis, so combining
%   terms amounts to adding sparse matrices.
%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PIETOOLS - unique_intmat
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

% Make sure the operator is properly specified
if ~Kop.chkval && ~isvalid(Kop)
    error("The functional is not properly declared.")
end

% Extract the parameters and intervals of integration
params = Kop.params;
intmat = Kop.intmat;
N = numel(Kop.intmat);
sz_K = size(params,1:N);
full_idcs = cellfun(@(a) 1:a,num2cell(sz_K),'UniformOutput',false);
for i=1:N
    % Determine unique intervals of integration
    intmat_i = intmat{i};
    [P,intmat_i] = uniquerows_integerTable(intmat_i);
    if size(intmat_i,1)==size(intmat{i},1)
        continue
    end
    % Combine terms associated with the same integral
    sz_K(i) = size(intmat_i,1);
    params_new = cell([sz_K,1]);
    old_idcs = full_idcs;
    new_idcs = full_idcs;
    for j=1:size(P,1)
        idx_new = find(P(j,:),1,'first');
        old_idcs{i} = j;
        new_idcs{i} = idx_new;
        params_j_old = params(old_idcs{:});
        params_j_new = params_new(new_idcs{:});
        params_j_new = cellfun(@add_or_set,params_j_new,params_j_old,'UniformOutput',false);
        params_new(new_idcs{:}) = params_j_new;
    end
    % Update the parameters and intmat
    params = params_new;
    intmat{i} = intmat_i;
    full_idcs{i} = 1:sz_K(i);
end
% Set the new intervals and parameters
Kop.params = params;
Kop.intmat = intmat;
Kop.chkval = true;

end


function c = add_or_set(a,b)
% Add b to a, or set c=b if a has not been set
if isempty(a)
    c = b;
else
    c = a + b;
end
end
