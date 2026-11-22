function vars_out = get_tensorvars(vars,degs,opts)
% VARS_OUT = GET_TENSORVARS(VARS,DEGS,OPTS) returns the names of the copies
% si_1,...,si_di of spatial variables si, the variables of a function of
% degree di in si, e.g. a tensor product of di factors depending on si.
%
% INPUTS
% - vars:   1 x N or N x 1 'cellstr', or 'polynomial' vector of N 'pvar'
%           objects, specifying the spatial variables si. May also be a
%           'char' for N=1;
% - degs:   1 x N or N x 1 array of nonnegative integers, the number of
%           copies di of each variable si. A scalar is used for all
%           variables;
% - opts:   (optional) 'char' or 'cellstr' of options:
%             'dum':          return the dummy variables
%                             si_dum1,...,si_dumdi instead;
%             'full':         return the copies of all variables in a
%                             single array, see OUTPUTS;
%             'bare':         return si (si_dum) instead of si_1
%                             (si_dum1) for di=1, the name used where
%                             there are no copies, see NOTES;
%             'pvar' or 'polynomial': return 'polynomial' objects;
%             'cellstr' or 'cell':    return 'cellstr' objects;
%           Without 'pvar' or 'cellstr', the output is of type 'polynomial'
%           if VARS is, and 'cellstr' otherwise;
%
% OUTPUTS
% - vars_out:   1 x N cell, with element i the 1 x di 'cellstr'
%               {'si_1',...,'si_di'}, or the di x 1 'polynomial' vector of
%               these variables. With option 'full', the 1 x sum(degs)
%               'cellstr' (sum(degs) x 1 'polynomial') concatenating these
%               elements, ordered by variable and then by copy;
%
% NOTES
% - A single copy is named si_1 (si_dum1), never si (si_dum), so that the
%   name of a copy does not depend on the number of copies. The bare name
%   si (option 'bare') is the name of a variable without copies.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PIETOOLS - get_tensorvars
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
% For support, contact M. Peet, Arizona State University at mpeet@asu.edu
% DJ, 11/19/2026: Initial coding;
% DJ, 11/22/2026: Build the names without num2str (faster);

% Check the options
if nargin<3
    opts = cell(1,0);
elseif ischar(opts)
    opts = {opts};
elseif ~iscellstr(opts)
    error("Options must be specified as 'char' or 'cellstr'.")
end
is_valid = ismember(opts,{'dum','full','bare','pvar','polynomial','cellstr','cell'});
if ~all(is_valid)
    error("Unknown option '%s'.",opts{find(~is_valid,1)})
end
% Flags of the options; 'pvar' and 'cellstr' exclude each other
is_dum = ismember('dum',opts);
is_full = ismember('full',opts);
is_bare = ismember('bare',opts);
is_pvar_opt = any(ismember({'pvar','polynomial'},opts));
is_cell_opt = any(ismember({'cellstr','cell'},opts));
if is_pvar_opt && is_cell_opt
    error("Options 'pvar' and 'cellstr' cannot both be specified.")
end

% Convert the variables to a 1 x N cellstr; the output is of the same type
% as the input, unless specified otherwise
is_pvar_out = (isa(vars,'polynomial') || is_pvar_opt) && ~is_cell_opt;
if ischar(vars)
    vars = {vars};
elseif isa(vars,'polynomial')
    if ~isempty(vars) && ~ispvar(vars)
        error("Variables must be specified as 'cellstr' or vector of 'pvar' objects.")
    end
    vars = pvar2varname(vars);
elseif ~iscellstr(vars)
    error("Variables must be specified as 'cellstr' or vector of 'pvar' objects.")
end
vars = reshape(vars,1,[]);
N = numel(vars);

% Check the degrees, using a scalar degree for all variables
if ~isnumeric(degs) || any(degs(:)<0) || any(degs(:)~=round(degs(:)))
    error("Degrees must be specified as nonnegative integers.")
elseif isscalar(degs)
    degs = repmat(degs,1,N);
elseif numel(degs)~=N
    error("Number of degrees must match the number of variables.")
end

% Declare the names si_j (or si_dumj) of the copies j=1,...,di of each si,
% or si (si_dum) for di=1 with option 'bare'
if is_dum
    sep = '_dum';       sep_bare = '_dum';
else
    sep = '_';          sep_bare = '';
end
vars_out = cell(1,N);
for i=1:N
    if is_bare && degs(i)==1
        vars_out{i} = {[vars{i},sep_bare]};
    else
        % A loop, with char(48+j) the digit j for j<10: num2str made this
        % function ~20x slower, and it runs on every Kop.varnames call
        names_i = cell(1,degs(i));
        for j=1:degs(i)
            if j<10
                names_i{j} = [vars{i},sep,char(48+j)];
            else
                names_i{j} = [vars{i},sep,int2str(j)];
            end
        end
        vars_out{i} = names_i;
    end
end
if is_full
    vars_out = [cell(1,0),vars_out{:}];
end

% Return column vectors of 'pvar' objects
if is_pvar_out && is_full
    vars_out = polynomial(vars_out(:));
elseif is_pvar_out
    vars_out = cellfun(@(a) polynomial(a(:)),vars_out,'UniformOutput',false);
end

end
