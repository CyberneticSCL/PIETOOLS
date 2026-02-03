function C = zero_nopvar(A)
% O = zero_nopvar(A) returns a nopvar object with the same dom, and
% vars as A, but with all degree 0 monomials and zero scalar coeffients.
% Used for filling in empty cells when multiplying two polyopvars.
%
% INPUTS
% - A:     polyopvar object; 
%
% OUTPUTS
% - O:     nopvar object with zero coefficients and same dom as A.
%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PIETOOLS - zero_nopvar
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
% CR, 01/28/2026: Initial coding.
% CR,01/29/2026: Updated to just returning coefficients when N=1.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% % Define N variables and dummy variables.
% N=size(A.dom,1);
% var1_name=char(A.pvarname);
% var2_name=[var1_name,repmat('_dum',[N,1])];
% vars=polynomial(mat2cell(var1_name,ones(N,1),size(var1_name,2)));
% dvars=polynomial(mat2cell(var2_name,ones(N,1),size(var2_name,2)));

% Generate cell of 3x1 or 3^N dimensions.
N=1; % 01/29/2026.
if N==1
    C=cell(3,1);
else
    dims=3*ones(1,N);
    C=cell(dims);
end

% Fill C with zeros.
% m=1; n=1;
Z = 0.0;
[C{:}] = deal(Z);

% % create nopvar.
% O=nopvar();
% O.C=C;
% O.deg=0;
% O.dom=A.dom;
% O.dim=[m,n];
% O.vars=[vars,dvars];

end

