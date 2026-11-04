function [G,H,p,J,vars] = compute_canonical_BCs(n,E,F,dom,vars,opts)
% [G,H,p,J,vars] = COMPUTE_CANONICAL_BCS(n,E,F,dom,vars,opts) takes a
% matrix E and polynomial F defining BCs on domain dom=[a,b]
%   sum_{j=1}^{d} E_{i,j}*[D_{j}u_{j}(a); D_{j}u_{j}(b)]
%       + sum_{j=0}^{d} int_{a}^{b} F_{i,j}(s)*u_{j}(s) ds = 0,
% for i in {1,...,d}, where
%   D_{j} u = [u; (d/ds) u; ...; (d/ds)^{j-1} u]
% and converts them to canonical format
%   sum_{j=1}^{d} G_{i,j}*D_{j}u_{j}(a)
%       + sum_{j=0}^{d} int_{a}^{b} H_{i,j}(s)*(d/ds)^{j} u_{j}(s) ds = 0,
% for i in {1,...,d}, where now G and H are such that the first r rows of
% both matrices specify constraints on the boundary values of u directly
% (r=rank(G)), and the remaining rows specify constraints on the
% highest-order derivatives (d/ds)^{j} u_{j} for each j in {0,...,d}, i.e.
% on the fundamental state.
%
% [G,H,p,J,vars] = COMPUTE_CANONICAL_BCS(n,E,F,dom,vars,'nosplit') only
% expresses the conditions in terms of the boundary values and
% highest-order derivatives, without splitting them: the order of the
% conditions is not changed, so that J=I and p is empty, and the number q
% of conditions (rows of E and F) is not required to equal
% n_tot=sum((1:d).*n(2:end)).
%
% INPUT
% - n:  1 x d+1 array specifying the vector dimension of the different
%       state variables, n = (n0,...,nd);
% - E:  n_tot x 2*n_tot array of type 'double', n_tot=sum((1:d).*n(2:end)),
%       defining the contribution of the boundary states to the boundary
%       conditions. The columns of the block E_{i,j} correspond to
%       [D_{j}u_{j}(a); D_{j}u_{j}(b)], and the blocks are ordered by j.
%       For the option 'nosplit', E can have any number q of rows;
% - F:  n_tot x sum(n) array of type 'polynomial' or 'double' specifying
%       the integral term in the boundary conditions. Can depend on at most
%       1 variable. For the option 'nosplit', F must have q rows;
% - dom:    1x2 array specifying the upper and lower boundaries [a,b] of
%           the interval on which u(s) is defined;
% - vars:   (optional) 1x2 array of type 'polynomial', specifying the
%           spatial variable s and a corresponding dummy variable theta.
%           Can be set to [] to leave the variables unspecified, when opts
%           is specified;
% - opts:   (optional) string or cell of strings specifying options. The
%           only option is 'nosplit', setting the function to return G and
%           H without splitting the conditions;
%
% OUTPUT
% - G:  n_tot x n_tot array of type 'double' representing the canonical
%       form of the boundary conditions (q x n_tot for 'nosplit');
% - H:  n_tot x sum(n) array of type 'polynomial' representing the integral
%       term in the canonical form of the boundary conditions. Will be
%       defined in terms of the same variable as the input F (q x sum(n)
%       for 'nosplit');
% - p:  1 x r array for r<=n_tot specifying the pivots of matrix G,
%       corresponding to elements of D_{j}u_{j}(a) which are actually
%       constrained by the boundary conditions. Empty for 'nosplit';
% - J:  n_tot x n_tot invertible array of type 'double' representing the
%       row operations such that Enew = J*E and Fnew = J*F. The identity
%       for 'nosplit';
% - vars:   1x2 array of type 'polynomial' specifying the spatial variable
%           s and dummy variable theta, as input, or as inferred from F (or
%           the default variables of 'opvar') if not specified;
%
% NOTES
% For E_{j} = [E_{j,a},E_{j,b}] the columns of E associated to u_{j}, and
% F_{j} the columns of F associated to u_{j}, the boundary conditions
% defined by E and F can be equivalently expressed as
%   G(E,F)*[D_{1}u_{1}(a);...;D_{d}u_{d}(a)]
%         = -int_{a}^{b} H(E,F)(s)*[u_{0};u_{1}';...;(d/ds)^{d}u_{d}](s) ds
% where
%     G(E,F) = [G_{1},...,G_{d}],
%     H(E,F) = [F_{0},H_{1},...,H_{d}],
%     G_{j} = E_{j,a} + E_{j,b}*Qj(b-a) + int_{a}^{b} F_{j}(s)*ej1(s-a) ds,
%     H_{j}(s) = E_{j,b}*ejj(b-s) + int_{s}^{b} pj_1(t-s)*F_{j}(t) dt,
% for the known Taylor matrices Qj (block upper triangular, with blocks
% p_{k-i}(z)*I) and ej1, ejj (first block row and last block column of Qj).
% These maps are linear in (E,F), and defined for any number q of
% conditions. The canonical form of the boundary conditions is such that
% G(Enew,Fnew) is in row-echelon form, with pivots in columns specified in
% p, and such that the last rows of G(Enew,Fnew) are zero.

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


% % Process the options
if nargin<6
    opts = {};
end
if ischar(opts) || isstring(opts)
    opts = {char(opts)};
end
if ~iscellstr(opts)
    error("Options must be specified as a string or cell of strings.")
end
do_split = true;
for k=1:numel(opts)
    switch lower(opts{k})
        case 'nosplit'
            do_split = false;
        case 'split'
            do_split = true;
        otherwise
            error("Unknown option '"+string(opts{k})+"'; the only supported option is 'nosplit'.")
    end
end

% % Process the input arguments
% The spatial variables are optional
if nargin<5
    vars = [];
end
has_vars = ~isempty(vars);% Check the state dimensions
if ~isnumeric(n) || ~isvector(n) || numel(n)<2 || any(n<0) || any(n~=round(n))
    error("Input 'n' must be a 1 x (d+1) array of nonnegative integers, with d>=1.")
end
n = n(:)';
d = numel(n)-1;
n_tot = sum((1:d).*n(2:end));       % number of boundary values
n_all = sum(n);                     % dimension of the state [u_0;...;u_d]
if n_tot==0
    error("At least one of the differentiable state components must be nonempty.")
end
% Starting indices of the blocks in E, G, and F (and H), for each component
idx_E0 = cumsum([0,2*(1:d).*n(2:end)]);
idx_G0 = cumsum([0,(1:d).*n(2:end)]);
idx_F0 = cumsum([0,n]);
% Check the E matrix
if ~isa(E,'double')
    error("Input 'E' must be an array of 'type 'double'.")
end
q = size(E,1);                      % number of conditions
if size(E,2)~=2*n_tot
    error("Input E should be a matrix with 2*n_tot columns, for n_tot = sum((1:d).*n(2:end)).")
end
% Check the F matrix
if isa(F,'double') || isa(F,'polynomial') && isdouble(F)
    F = double(F);
%     tmp = opvar();                                                        % MMP, 11/04/2026 (was)
%     s = tmp.var1;    s_dum = tmp.var2;                                    % MMP, 11/04/2026 (was)
    % the defaults of an empty opvar, without the class                     % MMP, 11/04/2026
    s = pvar('s1');    s_dum = pvar('s1_dum');                              % MMP, 11/04/2026
elseif isa(F,'polynomial')
    if ~isscalar(F.varname)
        error("Input 'F' can depend on at most one variable.")
    end
    s = polynomial(F.varname);
    s_dum = pvar([s.varname{1},'_dum']);
else
    error("Input 'F' must be of type 'double' or polynomial")
end
if ~all(size(F)==[q,n_all])
    error("Input F should be a matrix or polynomial of dimension q x sum(n), for q the number of rows of E.")
end
% Check the domain
if numel(dom)~=2
    error("Spatial domain must be specified as 1x2 array.")
end
a = dom(1);     b = dom(2);
if b<=a
    error("Upper boundary of spatial interval must be strictly greater than lower boundary.")
end
% Check the variables
if has_vars
    if ~isa(vars,'polynomial') || ~ispvar(vars)
        error("Spatial variables must be specified as 1x2 array of 'pvar' objects.")
    elseif numel(vars)~=2
        error("Spatial variables must be specified as 1x2 array.")
    end
    if isa(F,'polynomial')
        if strcmp(s.varname{1},vars(2).varname{1})
            F = subs(F,vars(2),vars(1));
            s = vars(1);
            s_dum = vars(2);
        elseif ~strcmp(s.varname{1},vars(1).varname{1})
            error("Variable appearing in F must match specified spatial variable.")
        end
    else
        s = vars(1);
        s_dum = vars(2);
    end
else
    vars = [s,s_dum];
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
F = polynomial(F);
G = zeros(q,n_tot);
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
    Gj = Ea + Eb*kron(Qj,eye(nj));
    for k=0:j-1
        cols_k = k*nj+1:(k+1)*nj;       % block column k of e_{j,1}
        Gj(:,cols_k) = Gj(:,cols_k) + double(int(Fj*(s-a)^k/factorial(k),s,a,b));
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
q = size(G,1);
n_tot = size(G,2);
J = eye(q);
p = zeros(1,0);
if ~do_split
    return
end
if q~=n_tot
    error("Input E should be a matrix of dimension n_tot x 2n_tot, for n_tot = sum((1:d).*n(2:end)).")
end

% % Convert to reduced row echelon form, by row operations J*G. Entries
% % below tol are regarded as zero, so that round-off in the row operations
% % cannot create spurious pivots (and hence overestimate the rank r).
tol = 1e-12*max(abs(G(:)));
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
    Pjj = eye(n_tot);
    Pii = -G(:,jj)/G(ii,jj);
    Pii(ii) = 1/G(ii,jj);
    Pjj(:,ii) = Pii;

    % Exchange rows to place row ii in position ii_strt
    Ujj = eye(n_tot-ii_strt+1);
    Ujj = Ujj([(ii+1:n_tot+1)-ii_strt,1:ii-ii_strt],:);
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
    if ii_strt>n_tot
        break
    end
end

% Apply the same row operations to H, which is linear in (E,F), so that it
% corresponds to the reordered BCs Enew = J*E, Fnew = J*F.
H = J*H;

end
