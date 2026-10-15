function Top = ndopvar2tensopvar(Pop,type)
% TOP = NDOPVAR2TENSOPVAR(POP,TYPE) returns a 'tensopvar' object TOP
% representing the same operator as the input 'nopvar' object POP, or the
% tensor-PI operator defined by a cell of 'nopvar' objects POP.
%
% INPUTS
% - Pop:    m x n 'nopvar' object representing a PI operator, or a K x d
%           cell of 'nopvar' objects specifying the factors (columns) in
%           each term (row) of a tensor-PI operator;
% - type:   (optional) 1 x 2 array of integers 1 or 2, specifying the type
%           of product of the factors to take along the row (type(1)) and
%           column (type(2)) dimensions, where 1 means an elementwise
%           product and 2 means a Kronecker product. A scalar value
%           specifies type(1), and type(2) is set to 2. Logical values are
%           also accepted, with false and true corresponding to 1 and 2.
%           Defaults to [1,2];
%
% OUTPUTS
% - Top:    m x n 'tensopvar' object representing the same operator as the
%           input;
%
% NOTES
% - Each 'nopvar' object is converted to a 'sopvar' object using
%   `nopvar2sopvar`, which cannot represent operators between different
%   function spaces. 'ndopvar' objects, with decision variables, are not
%   supported, as the factors of a 'tensopvar' object must be 'sopvar'
%   objects.

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
% MMP, 08/23/2026: point the three remaining 'ndopvar2tensopvar_new' call
%                  sites here, left over from the tensopvar_new rename
% DJ, 10/15/2026: Updated to use new polyopvar functions;

% Check the input
if ~isa(Pop,'nopvar') && ~isa(Pop,'cell')
    error("Operator must be specified as object of type 'nopvar', or as a cell of such objects.")
end
% By default assume pointwise product in row dimension, and Kronecker in
% the column dimension
if nargin<=1
    type = [1,2];
else
    if islogical(type)
        type = double(type)+1;
    end
    if isscalar(type)
        type = [type,2];
    end
end

% Convert the operators to 'sopvar' objects
if isa(Pop,'cell')
    for k=1:numel(Pop)
        if ~isa(Pop{k},'nopvar')
            error("Elements of the cell must all be of type 'nopvar'.")
        end
        Pop{k} = nopvar2sopvar(Pop{k});
    end
else
    Pop = nopvar2sopvar(Pop);
end

% Declare the tensor-PI operator
Top = sopvar2tensopvar(Pop,type);

end
