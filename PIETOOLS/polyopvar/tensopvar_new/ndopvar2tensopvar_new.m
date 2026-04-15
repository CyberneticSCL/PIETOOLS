function Top = ndopvar2tensopvar_new(Pop)
% TOP = NDOPVAR2TENSOPVAR(POP) returns a 'tensopvar_new' object TOP
% representing the same operator as the input 'nopvar' or 'ndopvar' object
% POP
%
% INPUTS
% - Pop:    m x n 'nopvar' or 'ndopvar' object representing a PI operator;
%
% OUTPUTS
% - Top:    m x n 'tensopvar' object representing the same operator as the
%           input;
%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PIETOOLS - ndopvar2tensopvar
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
% DJ, 04/14/2026: Initial coding

% Check the input
if ~isa(Pop,'nopvar') && ~isa(Pop,'ndopvar')
    error("Operator must be specified as object of type 'nopvar' or 'ndopvar'.")
end

% Declare a tensor-PI operator involving just a single factor
Top = tensopvar_new();
Top.ops = {Pop};
% Set spatial variables of output and input
Top.var1 = Pop.var1;
Top.var2 = Pop.var2;
% Set spatial domains of output and input variables
Top.dom1 = Pop.dom;
Top.dom2 = Pop.dom;
% Set the dimension of the operator
Top.dim = Pop.dim;

end