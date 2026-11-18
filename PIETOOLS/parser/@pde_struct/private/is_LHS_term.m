function [is_LHS,obj,idx,tdiff] = is_LHS_term(PDE_term)
% IS_LHS_TERM Checks if a 'struct' object "PDE_term" corresponds to a term
% in a 'pde_struct' object, belonging on the left-hand side of a PDE (i.e.
% involving an output signal or a temporal derivative of a state.
%
% INPUT
% - PDE_term:   'struct' object, representing a term in a 'pde_struct'
%               object.
% 
% OUTPUT
% - is_LHS:     Boolean variable specifying whether the term corresponds to
%               the left-hand side of the PDE, i.e. involves an output
%               signal or a temporal derivative of a state.
% - obj:        'char' object specifying whether the term involves a state
%               ('x'), a regulated output ('z'), or an observed output
%               ('y').
% - idx:        Integer specifying which component of type 'obj' the term
%               involves.
% - tdiff:      Integer specifying what order temporal derivative is taken
%               of the state component, if applicable.
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

if ~isa(PDE_term,'struct')
    error("Input term should be of type 'struct'")
end

% Assume the input corresponds to just a standard term.
is_LHS = false;
obj = '';
idx = 0;
tdiff = 0;

% Check if the term corresponds to a temporal derivative of a state, or to
% an output.
if isfield(PDE_term,'x') && isfield(PDE_term,'tdiff') && PDE_term.tdiff>0
    is_LHS = true;
    obj = 'x';
    idx = PDE_term.x;
    tdiff = PDE_term.tdiff;
elseif isfield(PDE_term,'y')
    is_LHS = true;
    obj = 'y';
    idx = PDE_term.y;
elseif isfield(PDE_term,'z')
    is_LHS = true;
    obj = 'z';
    idx = PDE_term.z;
end

end