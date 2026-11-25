function H = ctranspose(F)
%CTRANSPOSE Conjugate-transpose of quadPoly.
%
% If F(s,t) = (I_m ⊗ Zs(s)^T) * C * (I_n ⊗ Zt(t)),
% then H(s,t) = F(s,t)' is represented as
%   H(s,t) = (I_n ⊗ Zs(s)^T) * Cnew * (I_m ⊗ Zt(t)),
% using the SAME Zs/ns on the left and SAME Zt/nt on the right.
%
% Cnew is obtained by conjugate transposing only the outer (m x n) block
% structure of C, keeping the within-block (ds x dt) coordinates fixed.
%
% INPUTS
% - F:      m x n 'quadPoly' object;
%
% OUTPUTS
% - H:      n x m 'quadPoly' object, F' (coefficients conjugated).

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
% SS, 01/16/2026: Initial coding;

m  = F.dim(1);
n  = F.dim(2);

ds = prod(cellfun(@numel, F.Zs));  % tensor basis size on s-side, length(Zs)
dt = prod(cellfun(@numel, F.Zt));  % tensor basis size on t-side, length(Zt)

C = sparse(F.C);
[i, j, v] = find(C);   % find nonzero entries

% Within-block indices and block indices
a  = mod(i-1, ds) + 1;        % 1..ds
ib = floor((i-1) / ds);       % 0..m-1

b  = mod(j-1, dt) + 1;        % 1..dt
jb = floor((j-1) / dt);       % 0..n-1

% Swap ONLY the outer block indices (ib <-> jb)
i2 = jb*ds + a;             % rows: n blocks of size ds
j2 = ib*dt + b;             % cols: m blocks of size dt

Cnew = sparse(i2, j2, conj(v), n*ds, m*dt);

H = quadPoly(Cnew, F.Zs, F.Zt, [n, m], F.ns, F.nt);
end