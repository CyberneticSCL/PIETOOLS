function [logval,obj] = is_pde_var(PDE)
% LOGVAL = IS_PDE_VAR(PDE) Checks whether a pde_struct object PDE
% corresponds to a PDE variable (state variable, input variable, etc.).
%
% INPUT
% - PDE     pde_struct object;
%
% OUTPUT
% - logval  boolean variable specifying whether the pde_struct object
%           corresponds to just a PDE variable;
% - obj     char object specifying what type of PDE variable the passed
%           object is, either 'x', 'y', 'z', 'u', or 'w';
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Copyright (C)2024 PIETOOLS Team
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
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% If you modify this code, document all changes carefully and include date
% authorship, and a brief description of modifications
%
% DJ, 06/23/2024: Initial coding;
% DJ, 06/23/2025: Update to assume variable is specified as single term,
%                   see also update to "pde_var";

% Check that the input is indeed a 'pde_struct' object;
if ~isa(PDE,'pde_struct')
    error("Input argument must be of type 'pde_struct'.")
end

% Innocent until proven guilty;
logval = true;
obj = '';

% A PDE variable must correspond to precisely one type of object
objs = ['x';'w';'u';'y';'z'];
is_obj = false(5,1);
for jj=1:length(objs)
    is_obj(jj) = ~isempty(PDE.(objs(jj)));
end
if sum(is_obj)~=1
    logval = false;
    return
end
obj = objs(is_obj);

% A PDE variable cannot correspond to a boundary condition or multiple
% terms
if ~isempty(PDE.BC) || numel(PDE.free)>1
    logval = false;
    return
end

% A PDE variable cannot consist of multiple components
if numel(PDE.(obj))~=1
    logval = false;
    return
end

% A PDE variable cannot involve any terms as part of an equation
if isfield(PDE.(obj){1},'term') && ~isempty(PDE.(obj){1}.term)
    logval = false;
    return
end

% The only term in the PDE structure should just be the variable itself     % DJ, 06/23/2025
if isscalar(PDE.free) && isfield(PDE.free{1},'term')
    if ~isscalar(PDE.free{1}.term)
        % There can be only 1 term.
        logval = false;
        return
    end
    trm = PDE.free{1}.term{1};
    if ~isfield(trm,obj) || trm.(obj)~=1
        % Term must involve the current variable.
        logval = false;
        return
    elseif isfield(trm,'C') && (~isdouble(trm.C) ||...
            (~isequal(double(trm.C),1) && ~isequal(double(trm.C),eye(PDE.(obj){1}.size))))
        % Variable cannot be multiplied with anything.
        logval = false;
        return
    elseif isfield(trm,'D') && any(trm.D)
        % Variable cannot be differentiated.
        logval = false;
        return
    elseif isfield(trm,'loc') && ~isempty(trm.loc) ...
            && (any(size(trm.loc)~=size(PDE.(obj){1}.vars'))...
                || any(~isequal(polynomial(trm.loc),polynomial(PDE.(obj){1}.vars'))))
        % Variable cannot be evaluated at a particular position.
        logval = false;
        return
    elseif isfield(trm,'delay') && (~isdouble(trm.delay) || ~double(trm.delay)==0)
        % Variable cannot be delayed in time.
        logval = false;
        return
    elseif isfield(trm,'I') 
        % Variable cannot be integrated.
        for jj=1:numel(trm.I)
            if ~isempty(trm.I{jj})
                logval = false;
                return
            end
        end
    elseif isfield(trm,'tdiff') && ~trm.tdiff==0
        % Variable cannot be differentiated in time
        logval = false;
        return
    end
else
    logval = false;
    return
end

end