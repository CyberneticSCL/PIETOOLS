function [vars,d] = parse_diff_inputs(vars,d)
% [VARS,D] = PARSE_DIFF_INPUTS(VARS,D) checks and normalizes the variables
% and orders of differentiation specified for the 'diff' functions of the
% sopvar and tensor-PI classes.
%
% INPUTS
% - vars:   1 x r 'cellstr', char, or 1 x r 'polynomial' (pvar) array
%           specifying the names of the spatial variables to differentiate
%           with respect to. Names may be repeated, in which case the
%           orders of differentiation in the same variable are added;
% - d:      (optional) p x r array of nonnegative integers specifying the
%           order of the derivative in each of the variables, in each of p
%           derivatives. A scalar is used for all variables. Defaults to 1;
%
% OUTPUTS
% - vars:   1 x r' cellstr of the unique names, in sorted order, of the
%           variables to differentiate in at least one of the p
%           derivatives;
% - d:      p x r' array of the corresponding orders of differentiation;
%
% NOTES
% - r'=0, with d a p x 0 array, if all orders are zero;
% - Normalized inputs are returned unchanged if passed in again.
%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PIETOOLS - parse_diff_inputs
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

% Check that the variables are properly specified
if isa(vars,'polynomial')
    if ~ispvar(vars)
        error("Differentiation with respect to polynomial combinations of variables is not supported.")
    end
    vars = pvar2varname(vars);
elseif ischar(vars) || isstring(vars)
    vars = {char(vars)};
elseif ~iscellstr(vars)
    error("Variable names should be specified as 1 x r 'cellstr' or 'polynomial' array.")
end
if size(vars,1)~=1
    error("Variable names must be specified as row array.")
end
% Get unique variable names
[vars,~,idcs2] = unique(vars);
vars = vars(:)';

% Check that the order of differentiation is properly specified
if nargin<=1
    % Assume first-order derivative if no order is specified
    d = ones(1,numel(idcs2));
elseif isscalar(d)
    % Use the same order of differentiation in all variables
    d = d*ones(1,numel(idcs2));
elseif size(d,2)~=numel(idcs2)
    error("An order of differentiation must be specified for each of the spatial variables.")
end
if any(d(:)<0 | d(:)~=round(d(:)))
    error("Orders of differentiation must be nonnegative integers.")
end
% Add orders of derivatives in the same variable
d = full(d*sparse(1:numel(idcs2),idcs2,1));
% Remove derivatives of order zero
rm_cols = all(d==0,1);
vars = vars(~rm_cols);
d = d(:,~rm_cols);

end
