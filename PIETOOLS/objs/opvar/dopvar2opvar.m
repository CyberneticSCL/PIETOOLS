function opVar = dopvar2opvar(dopVar)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This function converst a dopvar class object to an opvar class object.
% 
% Input: 
% dopVar: a dopvar class object
% 
% Output:
% opVar: an opvar class object
%
% NOTES:
% For support, contact M. Peet, Arizona State University at mpeet@asu.edu
% or S. Shivakumar at sshivak8@asu.edu

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PIETools - dopvar2opvar
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



opvar opVar;
opVar.I = dopVar.I; opVar.dim = dopVar.dim;
opVar.var1 = dopVar.var1; opVar.var2 = dopVar.var2;
opVar.P = dpvar2poly(dopVar.P);
opVar.Q1 = dpvar2poly(dopVar.Q1);
opVar.Q2 = dpvar2poly(dopVar.Q2);
opVar.R.R0 = dpvar2poly(dopVar.R.R0);
opVar.R.R1 = dpvar2poly(dopVar.R.R1);
opVar.R.R2 = dpvar2poly(dopVar.R.R2);
end