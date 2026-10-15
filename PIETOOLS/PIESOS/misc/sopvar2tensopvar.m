function Top = sopvar2tensopvar(Pop,type)
% TOP = SOPVAR2TENSOPVAR(POP,TYPE) returns a 'tensopvar' object TOP
% representing the same operator as the input 'sopvar' object POP.
%
% INPUTS
% - Pop:    m x n 'sopvar' object representing a PI operator, or a cell
%           array of operators specifying the factors (columns) in each
%           term (row) of a tensor-PI operator;
% - type:   (optional) 1 x 2 array of integers 1 or 2, specifying the type
%           of product of the factors to take along the row (type(1)) and
%           column (type(2)) dimensions, where 1 means an elementwise
%           product and 2 means a Kronecker product. A scalar value
%           specifies type(1), and type(2) is set to 2. Defaults to
%           [1,2];
%
% OUTPUTS
% - Top:    m x n 'tensopvar' object representing the same operator as the
%           input;
%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PIETOOLS - sopvar2tensopvar
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
% DJ, 10/15/2026: Initial coding;

% Check the input
if ~isa(Pop,'sopvar') && ~isa(Pop,'cell')
    error("Operator must be specified as object of type 'sopvar'.")
end
% By default assume pointwise product in row dimension, and Kronecker in
% the column dimension
if nargin<=1
    type = [1,2];
elseif isscalar(type)
    type = [type,2];
end

if ~isa(Pop,'cell')
    % Declare a tensor-PI operator involving just a single factor
    Top = tensopvar();
    Top.ops = {Pop};
    Top.otype = type;
    % Set spatial variables and domain
    Top.vars = Pop.vars;
    Top.dom = Pop.dom;
    Top.dims = size(Pop);
    Top.depmat.in = ones(1,size(Top.vars.in,2));
    Top.depmat.out = ones(1,size(Top.vars.out,2));
    Top.chkval = true;
else
    % Assume a cell of terms and factors is specified
    Top = tensopvar();
    Top.ops = Pop;
    Top.otype = type;
    Top = initialize(Top);
end

end