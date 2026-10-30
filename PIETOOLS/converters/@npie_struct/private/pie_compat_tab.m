function out = pie_compat_tab(PIE,type,val)
% OUT = PIE_COMPAT_TAB(PIE,TYPE,VAL) returns the table describing the
% state, input, or output components of a PIE in the format used by earlier
% versions of the PIE structure, or stores the information in such a table
% in the new format of the PIE structure.
%
% INPUT
% - PIE:    'pie_struct' or 'npie_struct' object;
% - type:   One of 'x', 'w', 'u', 'y', 'z', specifying whether to extract
%           or set the table for the state, exogenous input, actuator
%           input, observed output, or regulated output;
% - val:    (optional) ncomps x (2+N) or ncomps x (2+2*N) array, with N the
%           number of spatial variables in PIE.space.vars. Row i of the
%           table provides
%               val(i,1):       index of the associated PDE component,
%                               stored as PIE.pde2pie.(type).idx(i,1);
%               val(i,2):       size of the component,
%                               stored as PIE.space.(type).dims(i);
%               val(i,3:2+N):   binary values indicating whether the
%                               component depends on each spatial variable,
%                               stored as PIE.space.(type).depmat(i,:);
%               val(i,3+N:2+2*N):   order of the derivative taken of the
%                               PDE component with respect to each spatial
%                               variable, stored as
%                               PIE.pde2pie.(type).diff(i,:). This is
%                               only used for type = 'x', 'w', 'u';
%
% OUTPUT
% - out:    The table, if only PIE and type are specified. For type = 'x',
%           the table always has 2+2*N columns. For type = 'w' or 'u', it
%           only has 2+2*N columns if PIE.pde2pie.(type).diff is
%           specified. For type = 'y' or 'z', it has 2+N columns.
%           Otherwise the input PIE with updated fields;
%
% NOTES
% The variables of the PIE must be specified before a table can be set.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PIETOOLS - pie_compat_tab
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

N = numel(PIE.space.vars);
has_diff = any(strcmp(type,{'x','w','u'}));

if nargin==2
    % % Extract the table
    if ~isfield(PIE.space,type)
        out = zeros(0,2+(1+strcmp(type,'x'))*N);
        return
    end
    dims = PIE.space.(type).dims(:);
    n = numel(dims);
    depmat = PIE.space.(type).depmat;
    if isempty(depmat)
        depmat = false(n,N);
    end
    % Index of the associated PDE component, defaults to 1:n.
    idx = (1:n)';
    p2p = PIE.pde2pie;
    has_map = isfield(p2p,type) && isstruct(p2p.(type));
    if has_map && isfield(p2p.(type),'idx') && size(p2p.(type).idx,1)==n && n>0
        idx = p2p.(type).idx(:,1);
    end
    out = [idx,dims,double(depmat)];
    % Order of the derivative taken of the PDE component.
    if has_diff
        if has_map && isfield(p2p.(type),'diff') && ~isempty(p2p.(type).diff)
            out = [out,p2p.(type).diff];
        elseif strcmp(type,'x')
            out = [out,zeros(n,N)];
        end
    end
    return
end

% % Set the table
n = size(val,1);
nc = size(val,2);
if isempty(val)
    n = 0;
    nc = 2+N;
end
if nc<2+N || (N==0 && nc>2)
    error(['The number of columns of PIE.',type,'_tab does not match the number of spatial variables in PIE.vars. Please specify PIE.vars before setting PIE.',type,'_tab.'])
end
if ~isfield(PIE.space,type) || ~isstruct(PIE.space.(type))
    PIE.space.(type) = struct('dims',zeros(0,1),'depmat',false(0,N),'varname',{cell(0,1)});
end
if isempty(val)
    val = zeros(0,nc);
end
PIE.space.(type).dims = val(:,2);
PIE.space.(type).depmat = logical(val(:,3:2+N));
% The names of the components no longer match if the components change.
if ~isfield(PIE.space.(type),'varname') || numel(PIE.space.(type).varname)~=n
    PIE.space.(type).varname = cell(0,1);
end
if strcmp(type,'x')
    % The first column indicates a PDE state, the others a PDE input.
    PIE.pde2pie.(type).idx = [val(:,1),zeros(n,2)];
else
    PIE.pde2pie.(type).idx = val(:,1);
end
if has_diff
    if nc>=2+2*N
        PIE.pde2pie.(type).diff = val(:,3+N:2+2*N);
    elseif strcmp(type,'x')
        PIE.pde2pie.(type).diff = zeros(n,N);
    else
        PIE.pde2pie.(type).diff = [];
    end
end
out = PIE;

end
