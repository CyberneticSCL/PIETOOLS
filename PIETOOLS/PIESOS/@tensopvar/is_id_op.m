function logval = is_id_op(Rop,ztol)
% LOGVAL = IS_ID_OP(ROP,ZTOL) checks whether the tensor-PI operator ROP
% is an identity operator.
%
% INPUTS
% - Rop:    m x n 'tensopvar' object;
% - ztol:   (optional) scalar tolerance below which coefficients are taken
%           to be zero, passed to 'clean_op' and to the 'is_id_op' of each
%           factor. Defaults to 1e-12;
%
% OUTPUTS
% - logval: boolean specifying whether (true) or not (false) Rop is an
%           identity operator on the function space on which it is defined;
%
% NOTES
% - For 'tensopvar' objects defined by multiple terms, the function does
%   not check if the terms add up to an identity operator. It only checks
%   if one of the terms is identity, and the remainder are zero.
%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PIETOOLS - is_id_op
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

% Set default tolerance, matching clean_op and the factors' is_id_op
if nargin==1
    ztol = 1e-12;
end

% Assume the operator is not identity until we have verified
logval = false;

% Make sure the input and output domains of Rop match
if ~isequal(Rop.dims(:,1),Rop.dims(:,2)) ||...
      ~isequal(Rop.vars.in,Rop.vars.out) ||...
        ~isequal(Rop.depmat.in,Rop.depmat.out) ||...
          (size(Rop.ops,2)>1 && Rop.otype(1)~=Rop.otype(2))
    return
end

% Get rid of zero terms in the operator
[Rop,is_zero]  = clean_op(Rop,ztol);
if is_zero || size(Rop.ops,1)~=1
    return
end

% Check that each of the PI operators defining Rop is an identity operator
for j=1:size(Rop.ops,2)
    is_id_j = is_id_op(Rop.ops{1,j},ztol);
    if ~is_id_j
        return
    end
end

% At this point, we can conclude that Rop is an identity operator
logval = true;


end