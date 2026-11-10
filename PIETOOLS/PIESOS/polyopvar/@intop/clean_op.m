function [Kop,is_zero_op] = clean_op(Kop,ztol)
% [KOP,IS_ZERO_OP] = CLEAN_OP(KOP,ZTOL) sets the coefficients of the
% 'intop' object KOP that are below the tolerance ZTOL to zero, and
% discards the terms in the functional that are zero.
%
% INPUTS
% - Kop:    'intop' object representing a PI functional;
% - ztol:   (optional) scalar value specifying the tolerance below which
%           coefficients are assumed to be zero. Defaults to 1e-12;
%
% OUTPUTS
% - Kop:    'intop' object representing the same functional as the input up
%           to the specified tolerance, with coefficients in Kop.params
%           of magnitudes smaller than ztol set to 0. If all terms
%           involving a particular integral in the functional (particular
%           row of Kop.intmat{i} for some i) are zero, this integral (row
%           in Kop.intmat{i}) is removed -- except the last remaining row
%           along a given direction, which is always kept (with its zero
%           coefficients) so Kop.matdim stays correct;
% - is_zero_op: true if every entry of Kop.params is zero up to ztol, i.e.
%           Kop represents the zero functional; false otherwise;
%
% NOTES
% - The coefficients are cleaned as 'P.*(abs(P)>ztol)', and tested for
%   being zero as 'nnz(P)==0', so that only the nonzeros of the sparse
%   matrix P are accessed.
%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PIETOOLS - clean_op
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

% Check that the functional is properly specified
if ~Kop.chkval && ~isvalid(Kop)
    error("The functional is not properly declared.")
end

% Set default tolerance
if nargin==1
    ztol = 1e-12;
end

% Set terms with zero parameters to 0
Kparams = Kop.params;
is_zero = false([size(Kparams),1]);
for k=1:numel(Kparams)
    param_k = Kparams{k};
    param_k = param_k.*(abs(param_k)>ztol);
    is_zero(k) = nnz(param_k)==0;
    Kparams{k} = param_k;
end
is_zero_op = all(is_zero(:));

% Check if any of the integrals can be removed
if ~any(is_zero(:))
    % Nothing to remove, but still return the parameters with entries below ztol zeroed
    Kop.params = Kparams;
    return
end
intmat = Kop.intmat;
full_idcs = [cellfun(@(a) 1:size(a,1),intmat,'UniformOutput',false),{1}];
for i=1:numel(intmat)
    intmat_i = intmat{i};
    tmp_idcs = full_idcs;
    j = 1;
    while j<=size(intmat_i,1)
        % Check if all parameters in term j along direction i are zero
        tmp_idcs{i} = j;
        is_zero_ij = is_zero(tmp_idcs{:});
        if all(is_zero_ij(:)) && size(intmat_i,1)>1
            % Remove term j along the ith direction
            rtn_idcs = full_idcs;
            rtn_idcs{i} = [(1:j-1),j+1:size(intmat_i,1)];
            Kparams = Kparams(rtn_idcs{:});
            is_zero = is_zero(rtn_idcs{:});
            intmat_i(j,:) = [];
            full_idcs{i} = 1:size(intmat_i,1);
        else
            j = j+1;
        end
    end
    Kop.intmat{i} = intmat_i;
end
Kop.params = Kparams;
Kop.chkval = true;

end
