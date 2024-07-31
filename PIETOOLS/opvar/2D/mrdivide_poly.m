function Xp = mrdivide_poly(Bp,Ap,res_tol,max_deg)
% X = MRDIVIDE_POLY(P2,P1)
% computes an (approximate) polynomial inverse Xp = Bp/Ap such that X*Ap=Bp
%
% INPUT:
% - Bp          nxm object of type 'double' or 'polynomial';
% - Ap          pxm object of type 'double' or 'polynomial';
% - res_tol     1x1 object of type 'double' representing the allowed error
%               in the coefficients of the residual Xp*Ap - Bp.
%               Defaults to 1e-8;
% - max_deg     Maximal cumulative degree of the monomials defining X.
%               Defaults to 10 times the maximum degree of A;
%
% OUTPUT:
% - Xp          nxp object of type 'polynomial'. If there exists a
%               polynomial P of cumulative degree at most max_deg in all
%               variables such that P*Ap=Bp, then Xp=P. Otherwise, Xp is
%               the polynomial P of lowest degree that satisfies
%                   ||Xp*Ap-Bp||_{L2} <= res_tol,
%               or of cumulative degree at most max_deg such that
%                   ||Xp*Ap-Bp||_{L2} = min_{P} ||P*Ap-Bp||_{L2},
%               whichever comes first. The maximal degree of the polynomial
%               Xp will be increased until either the res_tol condition is
%               satisfied, or the degree matches the value of "deg_max".
%
% NOTES:
% The polynomial Xp is computed as follows. Suppose Ap is of dimension nxm
% and Bp of dimension 1xm, so that Xp is of dimension 1xn.
% We expand Xp as
%   Xp = c_X * kron(eye(n),Z_X)
% for some vector of monomials Z_X, and vector of coefficients c_X. Then
%   Xp*Ap = c_X * kron(eye(n),Z_X) * Ap = c_X * Qp,
% where Qp is again some polynomial. We expand
%   Qp = C_Q*kron(eye(n),Z)    and     Bp = c_B*kron(eye(n),Z)
% for some vector of (unique) monomials Z, and coefficients C_Q and c_B.
% Then we must have
%   c_X * C_Q * kron(eye(nc),Z) = Xp*Ap = Bp = c_B*kron(eye(n),Z).
% Thus, we can compute the unknown coefficients c_X by solving 
%   c_X*C_Q = c_B,
% so by calling c_X = C_Q\c_B;
% Note that this equation may not be solvable exactly, in which case Matlab
% returns a least-squares approximation. This least-squares approximation
% does not necesseraly yield an optimal polynomial approximation of Xp that
% minimizes e.g. the L2-norm ||Xp*Ap-Bp||_{L_{2}}.
% However, if we were to expand the polynomials using an L2 orthogonal
% basis (rather than vector of monomials), we will get an optimal 
% approximation in terms of the L2 norm. 
% Such a modification is left for future updates.
%
% For support, contact M. Peet, Arizona State University at mpeet@asu.edu
% or D. Jagt at djagt@asu.edu

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Copyright (C)2022  M. Peet, S. Shivakumar, D. Jagt
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
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% If you modify this code, document all changes carefully and include date
% authorship, and a brief description of modifications
%
% Initial coding DJ - 07/31/2024


% Make sure inputs are polynomial
if isa(Ap,'double')
    Ap = polynomial(Ap);
elseif ~isa(Ap,'polynomial')
    error('Input "A" in call "B/A" must be of type polynomial')
end
if isa(Bp,'double')
    Bp = polynomial(Bp);
elseif ~isa(Bp,'polynomial')
    error('Input "B" in call "B/A" must be of type polynomial')
end
% Allow inverse to be computed by setting Bp=1
if all(size(Bp)==1) && isequal(Bp,1)
    Bp = polynomial(eye(size(Ap,2)));
end

if nargin<=2
    res_tol = 1e-8;
    max_deg = 10*max(max(sum(Ap.degmat,2)),max(sum(Bp.degmat,2))-max(sum(Ap.degmat,2)));
end

% Deal with case that denominator is just a matrix
if isdouble(Ap)
    Xp = polynomial(Bp/double(Ap));
    return
end


% Extract variables and dimensions of the polynomials Ap and Bp
varnames = unique([Ap.varname;Bp.varname]);
nvars = length(varnames);
[~,Avar_idx] = ismember(Ap.varname,varnames);
[~,Bvar_idx] = ismember(Bp.varname,varnames);

[nr_A,nc] = size(Ap);
[nr_B,nc_B] = size(Bp);
if nc_B~=nc
    error('Objects "A" and "B" in call "B/A" should have the same number of columns.')
end

% Set degree of inverse in each variable:
% --> start with just degree of the polynomial itself
A_deg = zeros(1,nvars);
A_deg(Avar_idx) = max(Ap.degmat,[],1);
% Extract degmat of B in terms of new variables
degmat_B = zeros(max(size(Bp.degmat,1),1),nvars);
if ~isempty(Bp.varname)
    degmat_B(:,Bvar_idx) = Bp.degmat;
end
B_deg = max(degmat_B,[],1);
X_deg = max(A_deg,B_deg-A_deg);


while any(X_deg<max_deg)

% Determine degrees of monomials in each variable
degmat_X = (0:X_deg(1))';
for j=2:nvars
    degmat_X = [repmat(degmat_X,[X_deg(j)+1,1]),kron((0:X_deg(j))',ones(size(degmat_X,1),1))];
end
% Get rid of rows that exceed allowed maximal degree
degmat_X(sum(degmat_X,2)>max_deg,:) = [];
nZ_X = size(degmat_X,1);

Z_X = polynomial(eye(nZ_X),degmat_X,varnames,[nZ_X,1]);
Z_X = kron(eye(nr_A),Z_X);


% Now, we multiply Xp*Ap = C_X * Z_X * Ap
% Here, Q=Z_X*Ap will be polynomial:
Qp = Z_X*Ap;

% % % Expand Qp = C_Q*kron(eye(nc),Z) and  Bp = C_B*kron(eye(nc),Z)
% First, make sure the degmat of Q is expressed in terms of all variables
[~,Qvar_idx] = ismember(Qp.varname,varnames);
degmat_Q = zeros(size(Qp.degmat,1),nvars);
degmat_Q(:,Qvar_idx) = Qp.degmat;

% Next, combine monomials of Q and B, so that Pmat'*Z_QB = [Z_Q;Z_B];
[Pmat,~] = uniquerows_integerTable([degmat_Q;degmat_B],'transpose');

% Extract coefficients of Q and B associated to combined monomial vector
Q_coeffs = Pmat(:,1:size(Qp.C,1))*Qp.C;
B_coeffs = Pmat(:,size(Qp.C,1)+1:end)*Bp.C;

% Then, reorganize coefficients into matrix C so that Q = C*kron(eye(nc),Z)
nZ_QB = size(Q_coeffs,1);
C_Q = zeros(nr_A*nZ_X,nc*nZ_QB);
% Fill in the elements of C_Q
for k=1:nr_A*nZ_X*nc
    % Determine which row and column of Qp the element k defines
    cidx = ceil(k/(nr_A*nZ_X));       % Which column of Qp are we considering?
    ridx = k-((cidx-1)*nr_A*nZ_X);    % Which row of Qp are qe considering?

    % Fill in coefficients defining this element Qp
    C_Q(ridx,(cidx-1)*nZ_QB+1:cidx*nZ_QB) = Q_coeffs(:,k);  
end

C_B = zeros(nr_B,nc*nZ_QB);
for k=1:nr_B*nc
    % Determine which row and column of Bp the element k defines
    cidx = ceil(k/nr_B);        % Which column of Bp are we considering?
    ridx = k-((cidx-1)*nr_B);   % Which row of Bp are qe considering?

    % Fill in coefficients defining this element Bp
    C_B(ridx,(cidx-1)*nZ_QB+1:cidx*nZ_QB) = B_coeffs(:,k);  
end



% Finally, we have C_X*C_Q*Z_QB = C_B*Z_QB
% --> C_X = C_B/C_Q;
C_X = C_B/C_Q;
Xp = cleanpoly(C_X*Z_X,1e-14);

% Check accuracy of the result
res = Xp*Ap - Bp;
C_res = max(max(abs(res.C)));
% If the residual exceeds the admitted tolerance, repeat with greater
% degree
if C_res<res_tol
    break
else
    X_deg = min(X_deg+1,max_deg);
end

end


% % % Some code to test...
% Z_QB = polynomial(eye(nZ_QB),degmat_QB,varnames,[nZ_QB,1]);
% Z_QB = kron(eye(nc_A),Z_QB);
% err_Q = Qp - C_Q*Z_QB;
% err_B = Bp - C_B*Z_QB;

end