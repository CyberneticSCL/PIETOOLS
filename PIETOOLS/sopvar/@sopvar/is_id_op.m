function logval = is_id_op(Pop,ztol)
% LOGVAL = IS_ID_OP(POP,ZTOL) checks whether the PI operator POP is an
% identity operator.
%
% INPUTS
% - Pop:    'sopvar' object representing a PI operator;
% - ztol:   scalar tolerance below which coefficients are taken to be
%           zero. Defaults to 1e-12;
%
% OUTPUTS
% - logval: boolean specifying whether (true) or not (false) Pop is the
%           identity operator on the function space on which it is defined.
%
% NOTES
% - Pop is identity iff it maps between the same space (same dimensions,
%   variables and domains), acts as a multiplier in every variable
%   (params{1} only; all other params zero), and params{1} represents the
%   multiplier I_p*1. The monomial bases may be larger than {1}, provided
%   all other coefficients are zero.
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
% DJ, 10/07/2026: Initial coding;

% Set default tolerance
if nargin==1
    ztol = 1e-12;
end

% Assume the operator is not identity until we have verified
logval = false;

% Input and output space must coincide. With equal variable lists S1 and
% S2 are empty, so all variables pass through (S3) and params has 3^N cells.
p = Pop.dims(1);
if p~=Pop.dims(2) || ~isequal(Pop.vars.in(:),Pop.vars.out(:)) ||...
        ~isequal(Pop.dom.in,Pop.dom.out)
    return
end

% Cells other than params{1} correspond to an integral in some variable
params = Pop.params;
for k=2:numel(params)
    pk = params{k};
    if ~isempty(pk) && any(abs(nonzeros(pk))>ztol)
        return
    end
end

% Locate the constant monomial in ZL, ZR.
nZ = [1,1];     idx0 = [0,0];
Zs = {Pop.ZL, Pop.ZR};
for j=1:2
    for k=1:numel(Zs{j})
        pos = find(Zs{j}{k}==0,1);
        if isempty(pos)
            return      % no constant monomial, so no identity
        end
        nZ(j) = nZ(j)*numel(Zs{j}{k});
        idx0(j) = idx0(j)*numel(Zs{j}{k}) + pos-1;
    end
end
idx0 = idx0+1;

% Compare params{1} with the identity operator in this basis
P1 = params{1};
if ~isequal(size(P1),p*nZ)
    return
end
rows = (0:p-1)'*nZ(1) + idx0(1);
cols = (0:p-1)'*nZ(2) + idx0(2);
D = P1 - sparse(rows,cols,1,p*nZ(1),p*nZ(2));
dv = nonzeros(D);
logval = isempty(dv) || max(abs(dv))<=ztol;

end
