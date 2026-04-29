function Kop_out = reorder_vars(Kop,var_order,pvarname)
% KOP_OUT = REORDER_VARS(KOP,VAR_ORDER,PVARNAME) takes a PI functional KOP
% and reorders the integrals defining this functional to match the order
% VAR_ORDER = [i1,...,id] of the state variables
%   KOP_OUT*(x_{i1} o ... o x_{id}) = KOP*(x1 o ... o xd)
% where o denotes the tensor product.
%
% INPUTS
% - Kop:        m x n 'intop' object representing a functional on a
%               degree-d distributed monomial
% - var_order:  1 x d array spcifying the new order of the state variables
%               in the distributed monomial;
% - pvarname:   (optional) 1 x p cell specifying the new names of the dummy
%               variables to use in the integral defining the functional;
%
% OUTPUTS
% - Kop_out:    m x n 'intop' object representing the same functional as
%               the input but now acting on the reordered monomial;
%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PIETOOLS - reorder_vars
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
% DJ, 04/21/2026: Initial coding

% % Check the inputs
if nargin<=1
    error("Not enough inputs.")
end
if ~isa(Kop,'intop')
    error("Functional must be specified as 'intop' object.")
end
nvars = numel(Kop.pvarname);
if ~isa(var_order,'double')
    error("Order of variables must be specified as  d x 1 array of integers.")
elseif numel(var_order)~=nvars || ~isempty(setdiff(1:nvars,var_order))
    error("A new index must be specified for each of the variables in the functional.")
end

% Reorder the names of the variables in the functional
Kop_out = Kop;
if ~isequal(var_order(:),(1:nvars)')
    if Kop.matdim(2)~=1
    error("Functionals with multiple columns are not currently supported.")
    end
    Kop_out.pvarname = Kop.pvarname(var_order);
    % Reorder columns of omat to match new ordering of variables
    for jj=1:numel(var_order)
        Kop_out.omat(Kop.omat==var_order(jj)) = jj;
    end
end


if nargin<=2
    return
end
% % If a desired set of variables is specified, replace the current
% % variable names
if isa(pvarname,'polynomial')
    if ~ispvar(pvarname)
        error("Variables must be specified as array of 'pvar' objects.")
    end
    pvarname = pvar2varname(pvarname);
elseif ~iscellstr(pvarname)
    error("Variables must be specified as 'cellstr' array.")
end

if isscalar(pvarname)
    % Declare vector of variable names based on single variable
    pvarname_new = cell(1,nvars);
    for j=1:nvars
        pvarname_new{j} = [pvarname{j},'_',num2str(j)];
    end
    pvarname = pvarname_new;
elseif numel(pvarname)~=nvars
    error("Number of variable names should match number of variables in the operator.")
end

% Check the variables appearing in the kernsl
Cparams = Kop_out.params;
Cvarname = cell(size(Cparams.varname));
for kk=1:numel(Cparams.varname)
    % Replace variables in kernel by new variable names
    old_var = Cparams.varname(kk);
    var_idx = strcmp(Kop_out.pvarname,old_var);
    Cvarname(kk) = pvarname(var_idx);
end
% Set the new variables of the functional
Cparams.varname = Cvarname;
Kop_out.params = Cparams;
Kop_out.pvarname = pvarname;


end