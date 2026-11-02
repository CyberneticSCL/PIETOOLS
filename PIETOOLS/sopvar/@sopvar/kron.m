function C = kron(A,B)
% C = KRON(A,B) returns the Kronecker product of a numeric matrix and a
% 'sopvar' object, treating the operator as a matrix of operators, so that
%       C = kron(M,P)   has the blocks   M(i,j)*P,
%       C = kron(P,M)   has the blocks   P(a,b)*M,
% as for the Kronecker product of matrices, e.g. kron(speye(n),P) is the
% operator that applies P to each of n copies of the input, and kron(P,
% speye(n)) the operator that applies P to the components with the same
% index i of n copies of the components of the input.
%
% INPUTS
% - A,B:    one is an m x n numeric matrix M (full or sparse), the other an
%           r x s 'sopvar' object P, from L_2^s[S1,S3] to L_2^r[S2,S3];
%
% OUTPUTS
% - C:      'sopvar' object from L_2^{n*s}[S1,S3] to L_2^{m*r}[S2,S3] (the
%           same spaces for both orders of the arguments), such that for
%           kron(M,P) the component (i-1)*r+a of the output is
%           sum_j M(i,j)*(P u_j)_a, with u_j the components (j-1)*s+(1:s)
%           of the input u, and for kron(P,M) the component (a-1)*m+i
%           of the output is sum_j M(i,j)*(P v_j)_a, with v_j the
%           components (b-1)*n+j, b=1,...,s, of u;
%
% NOTES
% - By definition of the sopvar class, the parameters of C only depend on
%   those of P and on M: the monomial bases and variables are unchanged.
%   The parameters have the components of the output and input outermost
%   in their row and column indices, and the monomials innermost, so that
%   for kron(M,P) the parameters are kron(M,Cp), for Cp those of P, and
%   for kron(P,M) the nonzero entries of Cp and M are combined to place the
%   components of M between those of P and the monomials.
% - A zero block of the parameters, [] or scalar 0, remains zero.
%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PIETOOLS - kron (sopvar)
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
% DJ, 11/01/2026: Initial coding.

% % % Check the inputs
if nargin~=2
    error("Exactly two inputs are required.")
end
if isa(A,'sopvar') && isa(B,'sopvar')
    error("The Kronecker product of two 'sopvar' objects is not supported.")
elseif isa(A,'sopvar')
    P = A;      M = B;      M_left = false;
elseif isa(B,'sopvar')
    P = B;      M = A;      M_left = true;
else
    error("At least one input must be a 'sopvar' object.")
end
if ~isnumeric(M) || ~ismatrix(M)
    error("The Kronecker product of a 'sopvar' object is only supported with a numeric matrix.")
end
M = sparse(M);
if isscalar(M)
    % kron(M,P)=kron(P,M)=M*P, with the dimensions unchanged
    C = M*P;
    return
end

% % % Lift the parameters
% The row of a parameter is (component, monomial), with the monomial index
% fastest, and likewise for the columns, so that for kron(M,P) the
% parameters follow from the Kronecker product of M with the parameters
M_nrows = size(M,1);        M_ncols = size(M,2);
nnZL_arr = prod(cellfun(@numel,P.ZL));
nnZR_arr = prod(cellfun(@numel,P.ZR));
dims = [P.dims(1)*M_nrows,P.dims(2)*M_ncols];
params_new = P.params;
if ~M_left
    [ridcs_M,cidcs_M,vals_M] = find(M);
    ridcs_M = ridcs_M(:)';      cidcs_M = cidcs_M(:)';      vals_M = vals_M(:)';
end
for lidx=1:numel(params_new)
    param = params_new{lidx};
    if isempty(param) || (isscalar(param) && param==0)
        continue                            % zero block stays zero
    end
    if M_left
        params_new{lidx} = kron(M,param);
        continue
    end
    % kron(P,M): entry (a,ZL;b,ZR) of the parameter, with component a,
    % monomial index ZL, and component b, monomial index ZR, moves to
    % (a,i,ZL;b,j,ZR) for each nonzero M(i,j), scaled by M(i,j)
    [rridcs,ccidcs,vals] = find(param);
    rridcs = rridcs(:)-1;       ccidcs = ccidcs(:)-1;       vals = vals(:);
    ridcs = floor(rridcs/nnZL_arr);         ZL_idcs = mod(rridcs,nnZL_arr);
    cidcs = floor(ccidcs/nnZR_arr);         ZR_idcs = mod(ccidcs,nnZR_arr);
    rridcs_new = (ridcs*M_nrows+(ridcs_M-1))*nnZL_arr+ZL_idcs+1;
    ccidcs_new = (cidcs*M_ncols+(cidcs_M-1))*nnZR_arr+ZR_idcs+1;
    params_new{lidx} = sparse(rridcs_new(:),ccidcs_new(:),vals.*vals_M,dims(1)*nnZL_arr,dims(2)*nnZR_arr);
end
C = sopvar(params_new,P.vars,P.ZL,P.ZR,P.dom,dims);

end
