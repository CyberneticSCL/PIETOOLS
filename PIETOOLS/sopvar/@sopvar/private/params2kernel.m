function K = params2kernel(P,k)
% K = PARAMS2KERNEL(P,K) returns the kernel defined by the parameter
% P.PARAMS{K} of the 'sopvar' object P, as a 'polynomial' object.
%
% INPUTS
% - P:      'sopvar' object;
% - k:      integer, the linear index of the parameter in P.params;
%
% OUTPUTS
% - K:      m x n 'polynomial' object, in the output variables P.vars.out,
%           and in the dummy variables of the input variables P.vars.in,
%           named as the input variable with _dum appended. It only
%           involves the variables that appear in it. In a direction in
%           which the parameter is a multiplier, the kernel does not depend
%           on the dummy variable;
%
% NOTES
% - The operator is the sum over the parameters of the integrals of the
%   kernels K_k(s,s') times the function, see 'display' for the integrals
%   that correspond to each parameter.

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

[i,j,v,E,names] = decode_params(P,k);
m = P.dims(1);
n = P.dims(2);
if isempty(v)
    K = polynomial(zeros(m,n));
    return
end
used = any(E,1);
[dg,~,ic] = unique(E(:,used),'rows');
C = sparse(ic,(j-1)*m+i,v,size(dg,1),m*n);
K = polynomial(C,dg,reshape(names(used),[],1),[m,n]);

end
