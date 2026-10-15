function Y = eval_polyopvar(P,xvals)
% Y = EVAL_POLYOPVAR(P,XVALS) evaluates the matrix-valued distributed
% polynomial P at the state values XVALS, returning P(x) = C*Z(x).
%
% TEST-ONLY evaluator, independent of 'subs': the tests of 'subs' (and
% of routines that call it) must not use apply_poly, which evaluates by
% 'subs'. The inputs of each monomial are built here and the operators
% applied by apply_op of 'tensopmat', so the only shared code is the
% action of the operators themselves.
%
% INPUTS
% - P:      'polyopvar' object;
% - xvals:  p x 1 cell, one entry per state variable (matching
%           P.varname, in order), each a 'polynomial'/'double' column
%           vector of size P.dims.in(i);
%
% OUTPUTS
% - Y:      polynomial array representing P(x) evaluated at xvals;
%
% NOTES
% - For operators with a Kronecker (otype(1)=2) output, apply_op renames
%   the output variables s to s_1,... and these are not renamed back
%   (unlike apply_poly, for depmat.out<=1).

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PIETOOLS - eval_polyopvar
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
% DJ (CC), 10/15/2026: Initial coding;

if ~P.chkval
    P = initialize(P);
end
Z = full(P.Z);
nZ = size(Z,1);
Xcell = cell(1,nZ);
nfac = num_factors(P);
for k=1:nZ
    idcs = Z(k,Z(k,:)~=0);
    if nfac(k)>numel(idcs)
        % More factors than the degree: factors at the zeros of Z are multipliers
        % from R (dummy input 1); includes degree-0 monomials.
        Xk = num2cell(ones(1,nfac(k)));
        Xk(find(Z(k,:))) = reshape(xvals(idcs),1,[]);
        Xcell{k} = Xk;
    elseif isempty(idcs)
        % Degree-0 monomial: block expects a single dummy "1" input, not zero factors.
        Xcell{k} = {1};
    else
        Xcell{k} = xvals(idcs);
    end
end
Y = apply_op(P.C,Xcell);

end
