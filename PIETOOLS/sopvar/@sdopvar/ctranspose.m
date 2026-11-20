function At = ctranspose(A)
% AT = CTRANSPOSE(A) returns the adjoint A' of the 'sdopvar' object A with
% respect to the L2 inner product, affine in the same decision variables:
% in a direction with an integral the lower and upper integrals exchange and
% the monomial bases swap; in a multiplier direction the bases stay, since
% delta(s-s') collapses the kernel onto the diagonal (canonical multiplier
% form).
%
% INPUTS
% - A:  'sdopvar' object representing a PI operator from L_2^p[S1,S3] to
%       L_2^q[S2,S3];
%
% OUTPUTS
% - At: 'sdopvar' object representing the adjoint operator, from
%       L_2^q[S2,S3] to L_2^p[S1,S3], with the input and output variables
%       and domains exchanged and the same decision variables A.Zd;
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PIETOOLS - ctranspose
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
% For support, contact M. Peet, Arizona State University at mpeet@asu.edu
%
% Initial coding DB - 08/17/2026

% Number of shared variables, one cell direction each
nvars3 = numel(intersect(A.vars.in,A.vars.out));
% Cell alpha of A becomes cell -alpha of At: the lower (2) and upper (3)
% integrals exchange in every direction, the multiplier (1) stays, so
% position [1 2 3] maps to [1 3 2]
idx = repmat({[1 3 2]},1,nvars3);
params = A.params;
if numel(params.A)~=3^nvars3
    error("An 'sdopvar' object should hold 3^n3 parameters, one per "...
          +"multi-index in {1,2,3}^n3.")
end
if nvars3>0
    lin = reshape(1:3^nvars3,[3*ones(1,nvars3),1]);
    adjidx = lin(idx{:});    % where the adjoint of parameter k belongs
else
    adjidx = 1;
end

% The coefficients are rearranged, not merely transposed: the bases swap in
% an integral direction and stay in a multiplier direction, which
% 'canonical_adjoint_map' encodes as one map Tcell{ii} per cell acting on
% vec(C); B maps by the transpose of the same matrix.
[Tcell,ZLt,ZRt] = canonical_adjoint_map(A.vars,A.ZL,A.ZR,A.dims);
nC_t = A.dims(2)*prod([cellfun(@numel,ZLt),1]) ...
       * A.dims(1)*prod([cellfun(@numel,ZRt),1]);
% q: the number of decision variables, read from a stored B parameter where
% there is one, so that an operator with none is handled the same way
q = numel(A.Zd);
for ii = 1:numel(A.params.B)
    if ~isempty(A.params.B{ii})
        q = size(A.params.B{ii},1);
        break
    end
end

% Map each cell of A to the cell adjidx(ii) of At; [] stands for a zero cell
pA = cell(size(params.A));      pB = cell(size(params.B));
for ii=1:numel(params.A)
    Aii = params.A{ii};         Bii = params.B{ii};
    if isempty(Aii)
        pA{adjidx(ii)} = sparse(nC_t,1);
    else
        pA{adjidx(ii)} = Tcell{ii}*sparse(Aii(:));
    end
    if isempty(Bii)
        pB{adjidx(ii)} = sparse(q,nC_t);
    else
        pB{adjidx(ii)} = Bii*Tcell{ii}.';
    end
end
params.A = pA;      params.B = pB;

% Exchange the input and output variables and domains; declare the adjoint
At_vars = struct('in',{A.vars.out},'out',{A.vars.in});
At_dom = struct('in',A.dom.out,'out',A.dom.in);
At = sdopvar(params,At_vars,A.Zd,ZLt,ZRt,At_dom,[A.dims(2),A.dims(1)]);
end

