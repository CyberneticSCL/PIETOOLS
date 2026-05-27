function Kop_out = rename_vars(Kop,pvarname)
% KOP_OUT = RENAME_VARS(KOP,PVARNAME) takes a PI functional KOP and renames
% the dummy variables appearing in the kernels defining KOP.
%
% INPUTS
% - Kop:        m x n 'intop' object representing a functional on a
%               degree-d distributed monomial;
% - pvarname:   1 x p cell specifying the new names of the dummy
%               variables to use in the integral defining the functional;
%
% OUTPUTS
% - Kop_out:    m x n 'intop' object representing the same functional as
%               the input but now defined in terms of the dummy variables
%               pvarname, so that Kop_out.pvarname = pvarname;
%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PIETOOLS - rename_vars
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
% DJ, 05/26/2026: Initial coding


% Check the first input
if ~isa(Kop,'intop')
    error("First input must be of type 'intop.")
end
nvars = size(Kop.pvarname,2);

% Check the second input
if isa(pvarname,'polynomial')
    if ~ispvar(pvarname)
        error("Variable names must be specified as 1 x n cell or 'polynomial' array of pvar objects.")
    end
    pvarname = pvar2varname(pvarname);
elseif ~iscellstr(pvarname)
    error("Variable names must be specified as 1 x n cell or 'polynomial' array of pvar objects.")
end
pvarname = pvarname(:)';
if numel(pvarname)~=numel(unique(pvarname))
    error("Variable names must be unique.")
elseif numel(pvarname)~=nvars
    error("Number of variables must match number of dummy variables in integral.")
end

% Rename the variables appearing in the parameter
old_varname = Kop.pvarname;
Kparams = Kop.params;
param_varname = Kparams.varname;
for i=1:numel(param_varname)
    % Replace the old variable with the new variable name
    idx = ismember(old_varname,param_varname(i));
    param_varname(i) = pvarname(idx);
end
Kparams.varname = param_varname;

% Declare the output operator
Kop_out = Kop;
Kop_out.pvarname = pvarname;
Kop_out.params = Kparams;

end