function logval = is_id_op(Kop,ztol)
% LOGVAL = IS_ID_OP(KOP,ZTOL) checks whether the functional KOP is an
% identity operator.
%
% INPUTS
% - Kop:    m x n 'intop' object;
% - ztol:   (optional) scalar tolerance below which the difference between
%           the kernel and the identity matrix is taken to be zero.
%           Defaults to 1e-12;
%
% OUTPUTS
% - logval: boolean specifying whether (true) or not (false) Kop is the
%           identity operator on R^m;
%
% NOTES
% - An 'intop' maps functions of N spatial variables to R^m. For N>0, it is
%   never the identity, as it maps from L2 to R^m. Hence, Kop is the
%   identity iff Kop.vars.in is empty and the coefficients of Kop.params{1}
%   are those of the m x m identity matrix: constant term 1 on the
%   diagonal, and no other nonzeros, in particular none for the decision
%   variables.
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

% Set default tolerance
if nargin==1
    ztol = 1e-12;
end

% Assume the operator is not identity until we have verified
logval = false;

% Only a map from R^m to R^m (no spatial variables) can be identity
if ~isempty(Kop.vars.in) || ~isempty(Kop.vars.out) || numel(Kop.params)~=1
    return
end
nd = numel(Kop.Zd);
P = Kop.params{1};
m = size(P,1)/(nd+1);
if m==0 || size(P,2)~=m
    return
end

% The kernel must equal the identity matrix: coefficient 1 at the constant-term
% row of each diagonal entry. The rows associated with decision variables are
% included in D, so these must vanish as well.
rows_const = (0:m-1)'*(nd+1) + 1;
Pid = sparse(rows_const,(1:m)',1,size(P,1),size(P,2));
D = P - Pid;
Cval = nonzeros(D);
logval = isempty(Cval) || max(abs(Cval))<=ztol;

end
