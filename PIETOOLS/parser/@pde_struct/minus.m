function PDE_out = minus(PDE_1,PDE_2)
% PDE_OUT = MINUS(PDE_1,PDE_2) builds a PDE representing the difference of 
% the terms in the PDE objects PDE_1 and PDE_2.
%
% INPUT
% - PDE_1:      'pde_struct' object representing either a set of equations
%               ([d/dt x=... ; y=...; z=...; 0=...]), or a set of free 
%               terms to be used to declare an equation.
% - PDE_2:      'pde_struct' object representing a set of free terms to 
%               subtract from the equations or terms specified by PDE_1. 
%
% OUTPUT
% - PDE_out:    'pde_struct' object representing either the same equations
%               as in PDE_1 but now with the terms from PDE_2 subtracted,
%               or a new set of free terms (collected in PDE_out.free)
%               corresponding to the difference of the terms in PDE_1 and 
%               PDE_2;
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

% % Check that the first input is of suitable type.
if isa(PDE_1,'state')
    PDE_1 = state2pde_struct(PDE_1);
elseif isa(PDE_1,'terms')
    PDE_1 = terms2pde_struct(PDE_1);
end
if isa(PDE_1,'polynomial') && isdouble(PDE_1)
    % Convert possible polynomial to double.
    PDE_1 = double(PDE_1);
end
if ~isa(PDE_1,'pde_struct') && (~isnumeric(PDE_1) || ~all(all(PDE_1==0)))
    % We can subtract zeros from PDEs or PDEs from zero, but anything
    % else is not supported!
    error("Subtracting non-'pde_struct' objects from 'pde_struct' objects is not supported.")
end

% % Check that the second input is of suitable type.
if isa(PDE_2,'state')
    PDE_2 = state2pde_struct(PDE_2);
elseif isa(PDE_2,'terms')
    PDE_2 = terms2pde_struct(PDE_2);
end
if isa(PDE_2,'polynomial') && isdouble(PDE_2)
    PDE_2 = double(PDE_2);
end
if isnumeric(PDE_2) && all(all(PDE_2==0))
    % PDE_1 - 0 = PDE_1;
    PDE_out = PDE_1;
    return
elseif isa(PDE_2,'pde_struct')
    % % Replace PDE_2 with -PDE_2, by changing the sign of all the 
    % % coefficients of all the terms.
    PDE_2 = uminus(PDE_2);
else
    error("Subtracting non-'pde_struct' objects from 'pde_struct' objects is not supported.")
end

% % Take the sum    PDE_1 +(-PDE_2);
PDE_out = plus(PDE_1,PDE_2);

end