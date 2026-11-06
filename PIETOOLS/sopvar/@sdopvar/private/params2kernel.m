function K = params2kernel(P,k)
% K = PARAMS2KERNEL(P,K) returns the kernel defined by the parameter
% P.PARAMS.A{K}, P.PARAMS.B{K} of the 'sdopvar' object P, as a 'dpvar'
% object.
%
% INPUTS
% - P:      'sdopvar' object;
% - k:      integer, the linear index of the parameter in P.params.A;
%
% OUTPUTS
% - K:      m x n 'dpvar' object, in the output variables P.vars.out, and in
%           the dummy variables of the input variables P.vars.in, named as
%           the input variable with _dum appended, and in the decision
%           variables P.Zd, with the coefficients A + B'*d of the parameter.
%           Only the decision variables and monomials that appear in the
%           kernel are included. In a direction in which the parameter is a
%           multiplier, the kernel does not depend on the dummy variable;
%
% NOTES
% - The operator is the sum over the parameters of the integrals of the
%   kernels K_k(s,s') times the function, see 'display' for the integrals
%   that correspond to each parameter.
% - The kernel is constructed with 'sdvar2dpvar'. This assumes the
%   canonical multiplier form of the object, in which a multiplier direction
%   has no dependence on the dummy variable, so that the dummy variable need
%   not be identified with the output variable.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PIETOOLS - params2kernel
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
% DJ, 11/06/2026: Initial coding;

m = P.dims(1);
n = P.dims(2);
% The decision variables as a column cell of char, also if they are a string array
if ischar(P.Zd)
    dvars = {P.Zd};
else
    dvars = cellstr(reshape(P.Zd,[],1));
end
nd = numel(dvars);

% The coefficients, with a zero block specified as [] or 0 replaced by zeros
ncoef = m*n*prod(cellfun(@numel,P.ZL))*prod(cellfun(@numel,P.ZR));
A = P.params.A{k};
B = P.params.B{k};
if nnz(A)==0 && nnz(B)==0
    K = dpvar(zeros(m,n));
    return
end
if numel(A)~=ncoef
    A = sparse(ncoef,1);
end
if numel(B)~=nd*ncoef
    B = sparse(nd,ncoef);
end

% The dummy variables are named as the input variables, with _dum appended
vars = struct('out',{reshape(P.vars.out,1,[])},'in',{strcat(reshape(P.vars.in,1,[]),'_dum')});
K = sdvar2dpvar(struct('A',A(:),'B',B),[m,n],vars,P.ZL,P.ZR,dvars);
K = compress(K);

end
