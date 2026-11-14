function [G,H,p,J] = compute_canonical_BCs(n,E,F,dom,vars,opts)
% [G,H,p,J] = COMPUTE_CANONICAL_BCS(n,E,F,dom,vars,opts) writes the
% q conditions on u=(u_0,...,u_d) on [a,b]
%   sum_{j=1}^{d} E_{i,j}*[D_{j}u_{j}(a); D_{j}u_{j}(b)]
%       + sum_{j=0}^{d} int_{a}^{b} F_{i,j}(s)*u_{j}(s) ds = 0,
% with D_{j}u = [u; (d/ds)u; ...; (d/ds)^{j-1}u], in the canonical form
%   sum_{j=1}^{d} G_{i,j}*D_{j}u_{j}(a)
%       + sum_{j=0}^{d} int_{a}^{b} H_{i,j}(s)*(d/ds)^{j}u_{j}(s) ds = 0,
% in the boundary values at a and the highest derivatives (the
% fundamental state). The conditions are reordered by row operations J,
% so that G is in reduced row echelon form: the first r=rank(G) rows
% fix boundary values, the others (zero rows of G) constrain the
% fundamental state. Option 'nosplit' skips the reordering (J=I).
%
% INPUT
% - n:      1 x (d+1) array, with n(j+1) the size of u_{j};
% - E:      q x 2*n_tot array, n_tot=sum((1:d).*n(2:end)), with columns
%           [D_{j}u_{j}(a); D_{j}u_{j}(b)] for each j (in this order);
%           q>=n_tot, except for 'nosplit';
% - F:      q x sum(n) 'polynomial' or 'double', in at most one
%           variable, s or its dummy;
% - dom:    1 x 2 array, the interval [a,b];
% - vars:   1 x 2 'pvar' array, the variable s and a dummy variable;
% - opts:   (optional) 'nosplit', to not reorder the conditions;
%
% OUTPUT
% - G:      q x n_tot array, the boundary values in the canonical form;
% - H:      q x sum(n) 'polynomial' in s, the integral terms;
% - p:      1 x r array, the pivot columns of G (the boundary values
%           fixed by the conditions); empty for 'nosplit';
% - J:      q x q invertible array, such that G and H are the canonical
%           form of the conditions J*E, J*F;
%
% NOTES
% With E_{j} = [E_{j,a},E_{j,b}] and F_{j} the columns of u_{j},
%   G = [G_{1},...,G_{d}],      H = [F_{0},H_{1},...,H_{d}],
%   G_{j} = E_{j,a} + E_{j,b}*Q_{j}(b-a)
%                   + int_{a}^{b} F_{j}(s)*e_{j,1}(s-a) ds,
%   H_{j}(s) = E_{j,b}*e_{j,j}(b-s) + int_{s}^{b} p_{j-1}(t-s)*F_{j}(t) dt,
% for p_{k}(z) = z^k/k!, Q_{j} the block upper triangular matrix with
% blocks p_{k-i}(z)*I, and e_{j,1}, e_{j,j} its first block row and last
% block column. G and H are linear in (E,F).

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PIETOOLS - compute_canonical_BCs
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
% DJ, 11/02/2026: Initial coding;
% MMP, 11/04/2026: For a constant F, the variables s and s_dum are declared
%   as s1 and s1_dum directly instead of read off an empty 'opvar', whose
%   defaults they are: the conversion of a pde_struct PDE no longer needs
%   the opvar class (removal of the opvar classes);
% DJ, 11/06/2026: Allow q>=n_tot conditions, instead of exactly n_tot, so
%   that conditions beyond the total order of the states become constraints
%   on the highest derivatives (zero rows of G, rows of H). The row
%   reduction acts on the q rows (was n_tot);
% DJ (CC), 11/13/2026: Require vars; drop the input checks and the
%   output vars; rank tolerance relative to the terms summed into G;

do_split = nargin<6 || ~strcmpi(opts,'nosplit');
n = n(:)';
d = numel(n)-1;
n_tot = sum((1:d).*n(2:end));       % number of boundary values
q = size(E,1);                      % number of conditions
a = dom(1);     b = dom(2);
% Starting indices of the blocks in E, G, and F (and H), for each component
idx_E0 = cumsum([0,2*(1:d).*n(2:end)]);
idx_G0 = cumsum([0,(1:d).*n(2:end)]);
idx_F0 = cumsum([0,n]);
% The spatial variable s and the dummy variable s_dum: those of vars if F
% is constant or in vars(2) (which is replaced by vars(1)), else the
% variable of F and <s>_dum
if isa(F,'double') || isdouble(F)
    F = polynomial(double(F));
    s = vars(1);    s_dum = vars(2);
elseif strcmp(F.varname{1},vars(2).varname{1})
    F = subs(F,vars(2),vars(1));
    s = vars(1);    s_dum = vars(2);
else
    s = polynomial(F.varname);
    s_dum = pvar([s.varname{1},'_dum']);
end


% % % Construct G and H
% % Declare G = [G_1,...,G_d] and H = [F_0,H_1,...,H_d], component by
% % component. For u_j with j>=1,
% %     G_j = E_{j,a} + E_{j,b}*Q_j(b-a)
% %           + int_{a}^{b} F_j(s)*e_{j,1}(s-a) ds,
% %     H_j(s) = E_{j,b}*e_{j,j}(b-s) + int_{s}^{b} p_{j-1}(t-s)*F_j(t) dt,
% % where p_k(z) = z^k/k!, Q_j(z) has blocks p_{k-i}(z)*I, e_{j,1} is the
% % first block row of Q_j, and e_{j,j} its last block column. For u_0,
% % there are no boundary values, and F_0 enters H directly.
G = zeros(q,n_tot);
G_scale = 0;                        % largest term summed into G
Hblk = cell(1,d+1);
if n(1)>0
    Hblk{1} = F(:,1:n(1));
end
for j=1:d
    nj = n(j+1);
    if nj==0
        continue
    end
    cols_E = idx_E0(j)+1:idx_E0(j+1);
    cols_u = idx_F0(j+1)+1:idx_F0(j+2);
    Ea = E(:,cols_E(1:j*nj));           % E_{j,a}: acts on D_{j}u_{j}(a)
    Eb = E(:,cols_E(j*nj+1:end));       % E_{j,b}: acts on D_{j}u_{j}(b)
    Fj = F(:,cols_u);

    % G_j, with Q_j(b-a) the Toeplitz matrix of p_{k-i}(b-a), k>=i
    Qj = toeplitz([1,zeros(1,j-1)],(b-a).^(0:j-1)./factorial(0:j-1));
    EbQ = Eb*kron(Qj,eye(nj));
    Gj = Ea + EbQ;
    G_scale = max([G_scale;abs(Ea(:));abs(EbQ(:))]);
    for k=0:j-1
        cols_k = k*nj+1:(k+1)*nj;       % block column k of e_{j,1}
        Fk = double(int(Fj*(s-a)^k/factorial(k),s,a,b));
        Gj(:,cols_k) = Gj(:,cols_k) + Fk;
        G_scale = max([G_scale;abs(Fk(:))]);
    end
    G(:,idx_G0(j)+1:idx_G0(j+1)) = Gj;

    % H_j, where block i of e_{j,j}(b-s) is p_{j-1-i}(b-s)*I
    Hj = polynomial(zeros(q,nj));
    for i=0:j-1
        Hj = Hj + Eb(:,i*nj+1:(i+1)*nj)*(b-s)^(j-1-i)/factorial(j-1-i);
    end
    Hj = Hj + int((s_dum-s)^(j-1)/factorial(j-1)*subs(Fj,s,s_dum),s_dum,s,b);
    Hblk{j+1} = Hj;
end
H = horzcat(Hblk{:});

% % % Split the conditions
J = eye(q);
p = zeros(1,0);
if ~do_split
    return
end
% More conditions than boundary values (q>n_tot) are allowed. The rows of
% G beyond its rank are then zero, and the same rows of H constrain v1

% % Convert to reduced row echelon form, by row operations J*G. Entries
% % below tol are regarded as zero, so that round-off in the row operations
% % cannot create spurious pivots (and hence overestimate the rank r). The
% % tolerance is relative to the terms summed into G, so that the round-off
% % of terms that cancel (a zero column of G) is not a pivot either.
tol = 1e-12*G_scale;                                                        % DJ (CC), 11/13/2026
ii_strt = 1;
for jj=1:n_tot
    % Find row in G with nonzero element in column jj
    ii = find(abs(G(ii_strt:end,jj))>tol,1,'first');
    if isempty(ii)
        % If there is no such row, continue to the next column
        continue
    end
    ii = ii+ii_strt-1;
    % Construct a row permutation matrix Pjj such that Pjj*G has all zeros
    % in column jj, except in element (ii,jj), which must be equal to 1;
    Pjj = eye(q);                                                           % DJ, 11/06/2026
    Pii = -G(:,jj)/G(ii,jj);
    Pii(ii) = 1/G(ii,jj);
    Pjj(:,ii) = Pii;

    % Exchange rows to place row ii in position ii_strt
    Ujj = eye(q-ii_strt+1);
    Ujj = Ujj([(ii+1:q+1)-ii_strt,1:ii-ii_strt],:);
    Ujj = blkdiag(eye(ii_strt-1),Ujj);
    Pjj = Ujj*Pjj;

    % Perform the actual permutation, and keep track of full permutation
    % matrix J;
    G = Pjj*G;
    G(abs(G)<=tol) = 0;
    J = Pjj*J;
    p = [p,jj];

    % Row ii is now linearly independent of all other rows
    ii_strt = ii_strt+1;
    if ii_strt>q
        break
    end
end

% Apply the same row operations to H, which is linear in (E,F), so that it
% corresponds to the reordered BCs Enew = J*E, Fnew = J*F.
H = J*H;

end
