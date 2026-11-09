function P_out = kronI(P_in,k)
% P_OUT = KRONI(P_IN,K) returns the 'nopvar' object representing the
% Kronecker product of the operator defined by 'nopvar' object P_IN and
% an identity matrix of dimension K x K
%
% INPUTS
% - P_in:   m  x n 'nopvar' object representing a PI operator Pop;
% - k:      scalar integer specifying the dimension of the identity
%           matrix;
%
% OUTPUTS
% - P_out:  m*k x n*k 'nopvar' object representing the operator 
%           (Pop o I_{k});

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PIETOOLS - kronI
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
% DJ, 01/15/2026: Initial coding


% Construct the permutation matrices (I_{m(d+1)} o P)^T and 
% (I_{n(d+1)} o P) for P defined such that
%   Z_{d}(s) o I_{k} = P*(I_{k} o Z_{d}(s))
deg = P_in.deg;
N = numel(deg);
d = prod(deg+1);
m = P_in.dim(1);   n = P_in.dim(2);

% Declare the (k,d) commutation matrix, so that
%   Z_{d}(s)' o I_{k} = (I_{k} o Z_{d}(s)')P_{k,d}'
Pt = commat(k,d,'transpose');
IP_m = spIkron(m,Pt);


% Construct the coefficients representing P_in o I_{q}
C_cell = P_in.C;
D_cell = cell(size(C_cell));
sz_C = size(C_cell);
for ii=1:numel(D_cell)
    % Determine the index of element ii along each dimension of the cell C
    idcs = cell(1,N);
    [idcs{:}] = ind2sub(sz_C,ii);
    idcs = cell2mat(idcs);
    % If element ii corresponds to an integral, we need to account for the
    % monomial basis in the associated dummy variable
    is_int = logical(idcs-1);
    d_tmp = prod(deg(is_int)+1);
    r_idcs = 1:k*d_tmp;
    c_idcs = (0:k-1)'*d_tmp + (1:d_tmp);
    P = sparse(r_idcs,c_idcs,1,k*d_tmp,k*d_tmp);
    IP_n = spIkron(n,P);

    D_cell{ii} = IP_m*kron(C_cell{ii},speye(k))*IP_n;
end
P_out = P_in;
P_out.C = D_cell;

end

