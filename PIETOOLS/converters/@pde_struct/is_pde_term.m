function logval = is_pde_term(PDE)
% LOGVAL = IS_PDE_TERM(PDE) Checks whether a pde_struct object PDE
% corresponds to a free term in a PDE.
%
% INPUT
% - PDE     pde_struct object;
%
% OUTPUT
% - logval  boolean variable specifying whether the pde_struct object
%           corresponds to just a term in a PDE;
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Copyright (C)2024  M. Peet, S. Shivakumar, D. Jagt
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
% Initial coding DJ - 06/23/2024

% Check that the input is indeed a 'pde_struct' object;
if ~isa(PDE,'pde_struct')
    error("Input argument must be of type 'pde_struct'.")
end

% Innocent until proven guilty;
logval = true;

% A PDE term must be specified in the field PDE.free
if isempty(PDE.free)
    logval = false;
    return
end

% Any other object in the PDE cannot involve any terms.
objs = {'x';'w';'u';'y';'z';'BC'};
for jj=1:length(objs)
    PDE_obj = PDE.(objs{jj});
    if isempty(PDE_obj)
        continue
    end
    for kk=1:numel(PDE_obj)
        if isfield(PDE_obj{kk},'term') && ~isempty(PDE_obj{kk}.term)
            logval = false;
            return
        end
    end
end

end