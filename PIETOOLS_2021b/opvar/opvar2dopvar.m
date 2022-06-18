function dopVar = opvar2dopvar(opVar)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This function converst a opvar class object to an dopvar class object
% provided the there are no product between the decision variables in the
% opvar properties
% 
% Input: 
% opVar: an opvar class object
% 
% Output:
% dopVar: a dopvar class object
%
% NOTES:
% For support, contact M. Peet, Arizona State University at mpeet@asu.edu
% or S. Shivakumar at sshivak8@asu.edu

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PIETools - opvar2dopvar
%
% Copyright (C)2021  M. Peet, S. Shivakumar
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
% Replaced poly2dpvar with just dpvar, DJ 12/30/2021.


dopvar dopVar;
dopVar.I = opVar.I; dopVar.dim = opVar.dim;
dopVar.var1 = opVar.var1; dopVar.var2 = opVar.var2;
dopVar.P = dpvar(opVar.P);
dopVar.Q1 = dpvar(opVar.Q1);
dopVar.Q2 = dpvar(opVar.Q2);
dopVar.R.R0 = dpvar(opVar.R.R0);
dopVar.R.R1 = dpvar(opVar.R.R1);
dopVar.R.R2 = dpvar(opVar.R.R2);
end