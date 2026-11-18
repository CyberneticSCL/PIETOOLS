function out = pie_compat_vars(PIE,val)
% OUT = PIE_COMPAT_VARS(PIE,VAL) returns the spatial variables of a PIE in
% the format used by earlier versions of the PIE structure, or stores
% spatial variables in the new format of the PIE structure.
%
% INPUT
% - PIE:    'pie_struct' or 'npie_struct' object, storing the names of the
%           N spatial variables as 1 x N cellstr PIE.space.vars;
% - val:    (optional) N x 2 or N x 1 object of type 'polynomial', with
%           the spatial variables in the first column, or a cellstr or
%           char specifying the names of the spatial variables, or [];
%
% OUTPUT
% - out:    N x 2 object of type 'polynomial' with the spatial variables in
%           the first column and the associated dummy variables in the
%           second, if only PIE is specified. Otherwise the input PIE, with
%           PIE.space.vars updated to the variables in val;
%
% NOTES
% The dummy variables are not stored, but named by adding '_dum' to the
% name of the associated spatial variable. If dummy variables are specified
% in val, they are ignored.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PIETOOLS - pie_compat_vars
%
% Copyright (C)2026  PIETOOLS Team
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
% DJ, 10/29/2026: Initial coding;

if nargin==1
    % % Extract the variables
    names = PIE.space.vars;
    N = numel(names);
    if N==0
        out = polynomial(zeros(0,2));
        return
    end
    dums = strcat(names,'_dum');
    out = polynomial([reshape(names,N,1),reshape(dums,N,1)]);
    return
end

% % Set the variables
if isempty(val)
    names = cell(1,0);
elseif ischar(val)
    names = {val};
elseif iscellstr(val)
    names = reshape(val,1,[]);
elseif isa(val,'polynomial') && ispvar(val)
    N = size(val,1);
    names = cell(1,N);
    for ii=1:N
        names{ii} = val(ii,1).varname{1};
    end
else
    error("Spatial variables should be specified as N x 2 'polynomial' object, or as cellstr of variable names.")
end
PIE.space.vars = names;
out = PIE;

end
