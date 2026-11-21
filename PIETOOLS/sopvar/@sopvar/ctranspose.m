function At = ctranspose(A)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% AT = CTRANSPOSE(A) returns the adjoint A' of the 'sopvar' object A with
% respect to the L2 inner product. For A: L_2^p[S1,S3] -> L_2^q[S2,S3] of
% the form
%   sum_alpha int_S1 int_S3' I_alpha(S3-S3') Z_L(S2,S3)' C_alpha Z_R(S3',S1),
% the adjoint At: L_2^q[S2,S3] -> L_2^p[S1,S3] has the kernel of the cell
% -alpha transposed: in a direction with an integral the lower and upper
% integrals exchange and the monomial bases swap; in a multiplier direction
% the bases do not swap, since delta(s-s') collapses the kernel onto the
% diagonal (canonical multiplier form).
%
% INPUT
% - A:      q x p 'sopvar' object;
%
% OUTPUT
% - At:     p x q 'sopvar' object representing A', with vars.in = A.vars.out
%           and vars.out = A.vars.in.
%
% For support, contact M. Peet, Arizona State University at mpeet@asu.edu
% or S. Shivakumar at sshivak8@asu.edu

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PIETools - ctranspose
%
% Copyright (C)2026  M. Peet, S. Shivakumar
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
% Initial coding MMP, SS  - 1_16_2026

% Number of shared variables, one cell direction each
nvars3 = numel(A.vars_S3);

% Cell alpha of A becomes cell -alpha of At: with the stored labels
% {1,2,3} for {0,1,-1}, the lower (2) and upper (3) integrals exchange in
% every direction, so position [1 2 3] maps to [1 3 2]
idx = repmat({[1 3 2]}, 1, nvars3);
if numel(A.params)~=3^nvars3
    error("A 'sopvar' object should hold 3^n3 parameters, one per "...
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
% 'canonical_adjoint_map' encodes as one map Tcell{kk} per cell.
[Tcell,leftZ,rightZ] = canonical_adjoint_map(A.vars,A.ZL,A.ZR,A.dims);
nrow_t = A.dims(2)*prod([cellfun(@numel,leftZ),1]);
ncol_t = A.dims(1)*prod([cellfun(@numel,rightZ),1]);

% Map each cell of A to the cell adjidx(kk) of At
Ct = cell(size(A.params));
for kk = 1:numel(A.params)
    Ck = A.params{kk};
    % [] and scalar 0 are zero blocks of the size of the block
    if isempty(Ck) || (isscalar(Ck) && Ck==0)
        Ct{adjidx(kk)} = sparse(nrow_t,ncol_t);
    else
        Ct{adjidx(kk)} = reshape(Tcell{kk}*sparse(Ck(:)),nrow_t,ncol_t);
    end
end

% Exchange the input and output variables
At_vars = struct();
At_vars.in = A.vars.out;
At_vars.out = A.vars.in;

if ~isa(At_vars.in, 'cell')
    At_vars.in = {At_vars.in};
end
if ~isa(At_vars.out, 'cell')
    At_vars.out = {At_vars.out};
end

% Switch output and input domains/dimensions
At_dom  = struct('in',A.dom.out, 'out',A.dom.in);
At_dim = [A.dims(2), A.dims(1)];

% Declare the adjoint operator
At = sopvar(Ct, At_vars, leftZ, rightZ, At_dom, At_dim);

end
