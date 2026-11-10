function [Kop] = quad2lin_term_v2(Pmat,Lmon,Rmon,dom,var1,var2)
% KOP = QUAD2LIN_TERM_V2(PMAT,LMON,RMON,DOM,VAR1,VAR2) converts a
% distributed polynomial in quadratic form
%   V(x) = <LMON(x), PMAT*RMON(x)>_{L2},
% for
%   LMON(x) = C1*(x1^d1*...*xm^dm)(s)  and  RMON = C2*(x1^p1*...*xm^pm)
% to the linear form, computing a functional operator KOP such that
%   V(x) = KOP*(Z1(x) o Z2(x)),
% where Z1(x) = x1^d1*...*xm^dm and Z2(x) = x1^p1*...*xm^pm are the
% distributed monomials of LMON and RMON, and o denotes the tensor product.
%
% INPUTS
% - Pmat:   m x n 'double', 'polynomial' or 'dpvar' object parameterizing
%           the inner product, for m and n the numbers of rows of LMON
%           and RMON. A scalar is taken to be a multiple of the identity,
%           if m=n. An empty array is taken to be 1;
% - Lmon:   'polyopvar' object with m rows, representing a single
%           distributed monomial, or a constant 'double', 'polynomial' or
%           'dpvar' object. An empty array is taken to be 1;
% - Rmon:   'polyopvar' object with n rows, representing a single
%           distributed monomial, or a constant 'double', 'polynomial' or
%           'dpvar' object;
% - dom:    (optional) 1 x 2 array specifying the interval, [a,b], over
%           which integration is performed. Not used;
% - var1:   (optional) 1 x 1 'pvar' object specifying the variable used
%           in the integral defining the inner product. Not used;
% - var2:   (optional) d x 1 'pvar' object specifying the variables used
%           in the multi-integral defining the inner product in linear
%           format. Not used;
%
% OUTPUTS
% - Kop:    'intop' object representing the functional operator defining
%           the inner product in linear format;
%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PIETOOLS - quad2lin_term_v2
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
% CR, 09/07/2026: Initial coding - Main diff is introduction of
%  is_tensor_L, is_tensor_R, and parameter_product function.
%  This replaces lines such as
% % Ltmp = Ltmp.*subs(Lop_params{trm,k}{3},var1,vars_m(j))
% by
% % parameter_product(Ltmp, ...
% %     subs(Lop_params{trm,k}{3},var1,vars_m(j)),is_tensor_L);
% %
% which replaces the standard product by a kronecker product if required
% as indicated by is_tensor_L.
%
% DJ, 10/05/2026: Make sure that type(2)==2 (no pointwise product of
%                   input);
% DJ, 10/15/2026: Updated to use new polyopvar functions;

% Process the left and right factors
if isa(Rmon,'polynomial') || isa(Rmon,'dpvar') || isa(Rmon,'double')
    Rmon = constant2polyopvar(Rmon);
end
if isempty(Lmon) || isa(Lmon,'polynomial') || isa(Lmon,'dpvar') || isa(Lmon,'double')
    % Constant on the left: <c,P*Rmon>, with c a matrix (empty: 1)
    if isempty(Lmon),   Lmon = 1;   end
    Lmon = constant2polyopvar(Lmon);
end
if ~isa(Lmon,'polyopvar') || ~isa(Rmon,'polyopvar')
    error("Monomials terms in the inner product should be specified as 'polyopvar' objects.")
end

if size(Lmon.Z,1)~=1 || size(Rmon.Z,1)~=1
    error("The left and right polynomial in the inner product must be defined by exactly one monomial.")
end

% Check that the other arguments are properly specified
if isempty(Pmat)
    Pmat = 1;
end

% Extract the operators acting on the left and right monomials
CL = Lmon.C.ops{1,1};
CR = Rmon.C.ops{1,1};
m = size(CL,1);     n = size(CR,1);

% Take the inner product of the coefficient operators: the weight is passed
% to int, so that it is integrated if it varies in space
if isscalar(Pmat) && m==n
    Kop = int(otimes(CL,CR,[1,2]),Pmat*ones(1,m));                          % DJ, 10/14/2026
else
    if size(Pmat,1)~=m || size(Pmat,2)~=n
        error("Dimensions of the weight must match those of the left and right operators.")
    end
    % Kronecker product along the rows: element (i,j) of Pmat, the entry k=i+(j-1)*m of its
    % vectorization, weights row i of CL and row j of CR, in the order of repmat(CL,n,1) and
    % repelem(CR,m,1) (as in piesos_quadvar).
    Wij = vec_dpvar(Pmat);
    Kop = int(otimes(repmat(CL,n,1),repelem(CR,m,1),[1,2]),Wij);            % DJ, 10/14/2026
end

end

%% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% 
function W = vec_dpvar(P)                                                   % DJ, 10/14/2026
% W = VEC_DPVAR(P) returns the 1 x numel(P) 'dpvar' W such that W = P(:).'
% in the order of the linear index (column-major), computed from the
% coefficients of P. Copy of the local function of `piesos_quadvar` (DJ, 10/10/2026), as
% indexing P(perm) of a dpvar needs an array of size (nd+1)*m*n x (m*n).
if ~isa(P,'dpvar')
    W = P(1:end);
    if size(W,1)>1,     W = W.';    end
    return
end
[m,n] = size(P);
nd = numel(P.dvarname)+1;       nZ = size(P.degmat,1);
[r,c,v] = find(P.C);
i = ceil(r/nd);                 d = r-(i-1)*nd;
j = ceil(c/nZ);                 z = c-(j-1)*nZ;
k = i+(j-1)*m;
W = compress(dpvar(sparse(d,(k-1)*nZ+z,v,nd,m*n*nZ),P.degmat,P.varname,P.dvarname,[1,m*n]));
end
