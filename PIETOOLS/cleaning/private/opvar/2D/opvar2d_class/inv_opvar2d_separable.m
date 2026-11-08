function [Phat] = inv_opvar2d_separable(Pop,inv_tol)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% [Pinv] = inv(P) takes in one operator, 
% P: R^m0 x L2^mx x L2^my x L2^m2 --> R^n0 x L2^nx x L2^ny x L2^n2
% It returns the (left) inverse Pinv: 
% Phat: R^n0 x L2^nx x L2^ny x L2^n2 --> R^m0 x L2^mx x L2^my x L2^m2
% 
% INPUT
% P: opvar2d class object representing a separable PI operator mapping 
%       R^m0 x L2^mx x L2^my to R^n0 x L2^nx x L2^ny.
%   That is, P.dim = [n0,m0; nx,mx; ny,my; 0,0] and P.Rxx{2} = P.Rxx{3} and
%   P.Ryy{2} = P.Ryy{3}.
%
% OUTPUT
% Pinv: the matlab structure such that subsequent application Pinv*P*u for
% an arbitrary polynomial u:R^n0 x L2^nx(s1) x L2^ny(s2) returns u.
%
% NOTES:
% For support, contact M. Peet, Arizona State University at mpeet@asu.edu
% or D. Jagt at djagt@asu.edu

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PIETools - inv_opvar2d_separable
%
% Copyright (C)2024 PIETOOLS Team
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
% Initial coding DJ - 01/28/2024
% DJ, 01/16/2025: Bugfix in "expand_opvar2d_quadratic"

% Deal with empty operator...
if isempty(Pop)
    Phat = Pop;
    return
end

if nargin<=1
    inv_tol = 1e-10;
end

% Convert all operator parameters to class 'polynomial'
Pop = poly_opvar2d(Pop);


% % % Error checking: Make sure inversion is supported.

% % First, check if the operator is separable.
[is_sep,is_full_int,is_no_mult] = is_separable(Pop,1e-12);
% If the operator is not separable, return with a warning of what might
% be wrong.
if ~is_sep
    if ~is_full_int(1)
        % Check that the map L2[x]-->L2[x] is defined by full integral operator
        error('Inversion only supported for full integral operators; make sure that P.Rxx{2}=P.Rxx{3}')
    elseif ~is_full_int(2)
        % Check that the map L2[y]-->L2[y] is defined by full integral operator
        error('Inversion only supported for full integral operators; make sure that P.Ryy{2}=P.Ryy{3}')
    elseif ~is_full_int(3)
        % Check that the map L2[x,y]-->L2[x,y] is defined by full integral operator
        error('Inversion only supported for full integral operators; make sure that P.R22{2,2}=P.R22{3,2}=P.R22{2,3}=P.R22{3,3}')
    elseif ~is_full_int(4)
        % Check that the maps L2[x,y]-->L2[x] and L2[x]-->L2[x,y] are defined by
        % full integral operators
        error('Inversion only supported for full integral operators; make sure that P.Rx2{2}=P.Rx2{3} and P.R2x{2}=P.R2x{3}')
    elseif ~is_full_int(5)
        % Check that the maps L2[x,y]-->L2[y] and L2[y]-->L2[x,y] are defined by
        % full integral operators
        error('Inversion only supported for full integral operators; make sure that P.Ry2{2}=P.Ry2{3} and P.R2y{2}=P.R2y{3}')
    end

    % Check that off-diagonal multipliers are zero
    if ~is_no_mult(1)
        error('Inversion only supported for full integral operators; make sure that P.Rx2{1}=0 and P.R2x{1}=0')
    elseif ~is_no_mult(2)
        error('Inversion only supported for full integral operators; make sure that P.Ry2{1}=0 and P.R2y{1}=0')
    elseif ~is_no_mult(3)
        error('Inversion only supported for full integral operators; make sure that P.R22{2,1}=P.R22{3,1}=0 and P.R22{1,2}=P.R22{1,3}=0')
    end
end

% % Also check that the matrix R00:\R-->\R is indeed invertible
inv_tol_mat = 1e-14;
if Pop.dim(1,1)>0 && Pop.dim(1,2)>0
    %error('P.R00 matrix is empty; inversion not supported')
    sgm = svd(double(Pop.R00));
    if sgm(end)/sgm(1) <= inv_tol_mat
        error('P.R00 matrix is (close to) singular; operator is not invertible')
    elseif length(sgm) < size(Pop.R00,2)
        error('P.R00 matrix is of insufficient column rank; operator is not (left-)invertible')
    end
end



% % % Now then, let's actually compute the inverse:
% First, extract the diagonal multiplier elements from the operator,
%   M = diag(R00,Rxx{1},Ryy{1},R22{1,1});
Pop_alt = Pop;
Pop_alt.R00 = 0*Pop.R00;        Pop_alt.Rxx{1} = 0*Pop.Rxx{1};
Pop_alt.Ryy{1} = 0*Pop.Ryy{1};  Pop_alt.R22{1,1} = 0*Pop.R22{1,1};

% Next, decompose the remaining separable integral operator as
%   Pop - M = ZL_op*H*ZR_op, 
% where ZL_op and ZR_op are monomial multiplier and integral operators,
% respectively.
[H,ZL_op,ZR_op] = expand_opvar2d_quadratic(Pop_alt);


% Now, we have decomposed   Pop = M + ZL_op*H*ZR_op.
% Assume that similarly     Phat = Minv + Minv*ZL_op*Hhat*ZR_op*Minv;
% Then
% Phat*Pop = Minv*M + Minv*ZL_op*H*ZR_op + Minv*ZL_op*Hhat*ZR_op*Minv*M
%               + Minv*ZL_op*Hhat*ZR_op*Minv*ZL_op*H*ZR_op
%          = I +Minv*ZL_op*(H +Hhat +Hhat*K*H)*ZR_op
% where K = ZR_op*Minv*ZL_op;.
% Thus, we need H +Hhat*(I+K*H) = 0, and therefore Hhat = -H/(I+K*H);

% Construct the inverse of the multiplier component
R00_hat = pinv(double(Pop.R00));
Rxx_hat = mrdivide_poly(eye(Pop.dim(2,2)),Pop.Rxx{1,1},inv_tol,20);
Ryy_hat = mrdivide_poly(eye(Pop.dim(3,2)),Pop.Ryy{1},inv_tol,20);
R22_hat = mrdivide_poly(eye(Pop.dim(4,2)),Pop.R22{1,1},inv_tol,15);

Minv = opvar2d(fliplr(Pop.dim));
Minv.var1 = Pop.var1;   Minv.var2 = Pop.var2;   Minv.I = Pop.I;
Minv.R00 = R00_hat;     
Minv.Rxx{1} = Rxx_hat;
Minv.Ryy{1} = Ryy_hat;
Minv.R22{1,1} = R22_hat;

% Construct the full inverse Hhat = -H\(I+K*H)
K = ZR_op*(Minv*ZL_op);
K = double(K.R00);
Hhat_R00 = -H/(eye(size(K,1),size(H,2)) +K*H);
Hhat = opvar2d();
Hhat.dim = [size(Hhat_R00,1),size(Hhat_R00,2);0,0;0,0;0,0];
Hhat.var1 = Pop.var1;   Hhat.var2 = Pop.var2;   Hhat.I = Pop.I;
Hhat.R00 = Hhat_R00;

% Finally, build the inverse Phat = Minv + Minv*ZL_op*Hhat*ZR_op*Minv;
Phat = (Minv*ZL_op)*Hhat*(ZR_op*Minv);
Phat.R00 = Phat.R00 + R00_hat;
Phat.Rxx{1} = Phat.Rxx{1} + Rxx_hat;
Phat.Ryy{1} = Phat.Ryy{1} + Ryy_hat;
Phat.R22{1,1} = Phat.R22{1,1} + R22_hat;

end



%%
function [H,ZL_op,ZR_op] = expand_opvar2d_quadratic(Pop)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This function decompose the operator Pop = ZL_op*H*ZR_op,
% where ZL_op is a multiplier operator and ZR_op an integral operator. 
% This requires the operator Pop to be separable, expanding the different
% parameters as
%
% [(),  R0x, R0y, R02] = [I            ]^T [(),  H0x, H0y, H02] [I            ]
% [Rx0, Rxx, Rxy, Rx2]   [  Zx0        ]   [Hx0, Hxx, Hxy, Hx2] [  Z0x        ]
% [Ry0, Ryx, Ryy, Ry2]   [      Zy0    ]   [Hy0, Hyx, Hyy, Hy2] [      Z0y    ]
% [R20, R2x, R2y, R22]   [          Z20]   [H20, H2x, H2y, H22] [          Z02]
%
% Here matrices Z.. each comprise a basis of vector-valued polynomials in
% (s1,s2,t1,t2)
% 
% INPUT
%  - Pop:   mxn 'opvar2d' class object representing a separable 2D PI operator.
%           Must satisfy P.R00=0, P.Rxx{1}=0, P.Ryy{1}=0, P.R22{1}=0,
%           as well as P.Rx2{1}=0, P.R2x{1}=0, P.Ry2{1}=0, P.R2y{1}=0,
%           and P.R22{2,1}=P.R22{3,1}=P.R22{1,2}=P.R22{1,3}=0.
%           In addition, only full integral operator is supported so that
%           P.Rxx{3}=P.Rxx{2}, P.Ryy{3}=P.Ryy{2},
%           P.R{3,2}=P.R{2,3}=P.R{3,3}=P.R{2,2}, as well as
%           P.Rx2{3}=P.Rx2{2}, P.R2x{3}=P.R2x{2}, P.Ry2{3}=P.R2y{2},
%           P.R2y{3}=P.R2y{2}.
% 
% OUTPUT
% - H:      nZ_L x nZ_R object of type 'double', representing a matrix of
%           coefficients.
% - ZL_op:  m x nZ_L 'opvar2d' class object, representing a multiplier
%           operator    \R-->[\R;L2[x];L2[y];L2[x,y]]
% - ZR_op:  nZ_R x n 'opvar2d' class object, representing an integral
%           operator    [\R;L2[x];L2[y];L2[x,y]] --> \R
%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% If you modify this code, document all changes carefully and include date
% authorship, and a brief description of modifications
%
% Initial coding DJ - 01/28/2024  
% DJ, 01/16/2025: Bugfix expanding multiplier terms, dimensions of Pmat_R
%                   in call to "expand_polynomial_quadratic" now match;

% Extract the involved parameters from the operator
% Note that we assume full integral operator, so that e.g. Pop.Rxx{1}=0 and
% Pop.Rxx{2}=Pop.Rxx{3};
                R0x = Pop.R0x;     R0y = Pop.R0y;     R02 = Pop.R02;
Rx0 = Pop.Rx0;  Rxx = Pop.Rxx{2};  Rxy = Pop.Rxy;     Rx2 = Pop.Rx2{2};
Ry0 = Pop.Ry0;  Ryx = Pop.Ryx;     Ryy = Pop.Ryy{2};  Ry2 = Pop.Ry2{2};
R20 = Pop.R20;  R2x = Pop.R2x{2};  R2y = Pop.R2y{2};  R22 = Pop.R22{2,2};

s1 = Pop.var1(1).varname{1};  s1_dum = Pop.var2(1).varname{1};
s2 = Pop.var1(2).varname{1};  s2_dum = Pop.var2(2).varname{1};

nr = Pop.dim(:,1);  nc = Pop.dim(:,2);

% Extract the degrees of the monomials in each variable defining the
% different parameters R
R0x_dmat = extract_degs(R0x,{s1});           nZ_0x = size(R0x_dmat,1);
R0y_dmat = extract_degs(R0y,{s2});           nZ_0y = size(R0y_dmat,1);
R02_dmat = extract_degs(R02,{s1;s2});     nZ_02 = size(R02_dmat,1);

Rx0_dmat = extract_degs(Rx0,{s1});           nZ_x0 = size(Rx0_dmat,1);
Rxx_dmat = extract_degs(Rxx,{s1;s1_dum});     nZ_xx = size(Rxx_dmat,1);
Rxy_dmat = extract_degs(Rxy,{s1;s2});     nZ_xy = size(Rxy_dmat,1);
Rx2_dmat = extract_degs(Rx2,{s1;s1_dum;s2}); nZ_x2 = size(Rx2_dmat,1);

Ry0_dmat = extract_degs(Ry0,{s2});           nZ_y0 = size(Ry0_dmat,1);
Ryx_dmat = extract_degs(Ryx,{s2;s1});     nZ_yx = size(Ryx_dmat,1);
Ryy_dmat = extract_degs(Ryy,{s2;s2_dum});     nZ_yy = size(Ryy_dmat,1);
Ry2_dmat = extract_degs(Ry2,{s2;s1;s2_dum}); nZ_y2 = size(Ry2_dmat,1);

R20_dmat = extract_degs(R20,{s1;s2});         nZ_20 = size(R20_dmat,1);
R2x_dmat = extract_degs(R2x,{s1;s2;s1_dum});     nZ_2x = size(R2x_dmat,1);
R2y_dmat = extract_degs(R2y,{s1;s2;s2_dum});     nZ_2y = size(R2y_dmat,1);
R22_dmat = extract_degs(R22,{s1;s2;s1_dum;s2_dum}); nZ_22 = size(R22_dmat,1);

% Collect all degrees of monomials that we need to expand the operator in
% the proposed quadratic form
Zx_dmat = [Rx0_dmat; Rxx_dmat(:,1); Rxy_dmat(:,1); Rx2_dmat(:,1)];
Zth_dmat = [R0x_dmat; Rxx_dmat(:,2); Ryx_dmat(:,2); R2x_dmat(:,3)];
Zy_dmat = [Ry0_dmat; Ryx_dmat(:,1); Ryy_dmat(:,1);  Ry2_dmat(:,1)];
Znu_dmat = [R0y_dmat; Rxy_dmat(:,2); Ryy_dmat(:,2); R2y_dmat(:,3)];
Zxy_dmat = [R20_dmat; R2x_dmat(:,[1,2]); R2y_dmat(:,[1,2]); R22_dmat(:,[1,2])];
Ztn_dmat = [R02_dmat; Rx2_dmat(:,[2,3]); Ry2_dmat(:,[2,3]); R22_dmat(:,[3,4])];

% Keep track of which rows in each dmat correspond to which parameter in R
nZ_x = [nZ_x0; nZ_xx; nZ_xy; nZ_x2];    nZ_th = [nZ_0x; nZ_xx; nZ_yx; nZ_2x];
nZ_y = [nZ_y0; nZ_yx; nZ_yy; nZ_y2];    nZ_nu = [nZ_0y; nZ_xy; nZ_yy; nZ_2y];
nZ_xy = [nZ_20; nZ_2x; nZ_2y; nZ_22];   nZ_tn = [nZ_02; nZ_x2; nZ_y2; nZ_22];

nnZ_x = cumsum([0;nZ_x]);       nnZ_th = cumsum([0;nZ_th]);
nnZ_y = cumsum([0;nZ_y]);       nnZ_nu = cumsum([0;nZ_nu]);
nnZ_xy = cumsum([0;nZ_xy]);     nnZ_tn = cumsum([0;nZ_tn]);

% Get rid of redundant monomials, keeping track of how to map new monomials
% back to old ones:     Pmat'*Zx_new = Zx_old;
[Pmat_x,Zx_dmat] = uniquerows_integerTable(Zx_dmat,'transpose');
[Pmat_th,Zth_dmat] = uniquerows_integerTable(Zth_dmat,'transpose');
[Pmat_y,Zy_dmat] = uniquerows_integerTable(Zy_dmat,'transpose');
[Pmat_nu,Znu_dmat] = uniquerows_integerTable(Znu_dmat,'transpose');
[Pmat_xy,Zxy_dmat] = uniquerows_integerTable(Zxy_dmat,'transpose');
[Pmat_tn,Ztn_dmat] = uniquerows_integerTable(Ztn_dmat,'transpose');

% Compute actual monomials associated to each degree matrix
Zx0 = polynomial(eye(size(Zx_dmat,1)),Zx_dmat,{s1},[size(Zx_dmat,1),1]);
Zy0 = polynomial(eye(size(Zy_dmat,1)),Zy_dmat,{s2},[size(Zy_dmat,1),1]);
Z20 = polynomial(eye(size(Zxy_dmat,1)),Zxy_dmat,{s1;s2},[size(Zxy_dmat,1),1]);
Z0x = polynomial(eye(size(Zth_dmat,1)),Zth_dmat,{s1},[size(Zth_dmat,1),1]);
Z0y = polynomial(eye(size(Znu_dmat,1)),Znu_dmat,{s2},[size(Znu_dmat,1),1]);
Z02 = polynomial(eye(size(Ztn_dmat,1)),Ztn_dmat,{s1;s2},[size(Ztn_dmat,1),1]);

% Take kronecker products to match number of rows and columns of the input
Zx0 = kron(eye(nr(2)),Zx0);     Z0x = kron(eye(nc(2)),Z0x);
Zy0 = kron(eye(nr(3)),Zy0);     Z0y = kron(eye(nc(3)),Z0y);
Z20 = kron(eye(nr(4)),Z20);     Z02 = kron(eye(nc(4)),Z02);


% Finally, compute coefficients H such that
% [(),  R0x, R0y, R02] = [I            ]^T [(),  H0x, H0y, H02] [I            ]
% [Rx0, Rxx, Rxy, Rx2]   [  Zx0        ]   [Hx0, Hxx, Hxy, Hx2] [  Z0x        ]
% [Ry0, Ryx, Ryy, Ry2]   [      Zy0    ]   [Hy0, Hyx, Hyy, Hy2] [      Z0y    ]
% [R20, R2x, R2y, R22]   [          Z20]   [H20, H2x, H2y, H22] [          Z02]
H0x = expand_polynomial_quadratic(R0x,1,Pmat_th(:,nnZ_th(1)+1:nnZ_th(2)));
H0y = expand_polynomial_quadratic(R0y,1,Pmat_nu(:,nnZ_nu(1)+1:nnZ_nu(2)));
H02 = expand_polynomial_quadratic(R02,1,Pmat_tn(:,nnZ_tn(1)+1:nnZ_tn(2)));

Hx0 = expand_polynomial_quadratic(Rx0,Pmat_x(:,nnZ_x(1)+1:nnZ_x(2)),ones(1,size(Rx0.degmat,1)));        % DJ, 01/16/2025: 1 --> ones(1,size(Rx0.degmat,1))
Hxx = expand_polynomial_quadratic(Rxx,Pmat_x(:,nnZ_x(2)+1:nnZ_x(3)),Pmat_th(:,nnZ_th(2)+1:nnZ_th(3)));
Hxy = expand_polynomial_quadratic(Rxy,Pmat_x(:,nnZ_x(3)+1:nnZ_x(4)),Pmat_nu(:,nnZ_nu(2)+1:nnZ_nu(3)));
Hx2 = expand_polynomial_quadratic(Rx2,Pmat_x(:,nnZ_x(4)+1:nnZ_x(5)),Pmat_tn(:,nnZ_tn(2)+1:nnZ_tn(3)));

Hy0 = expand_polynomial_quadratic(Ry0,Pmat_y(:,nnZ_y(1)+1:nnZ_y(2)),ones(1,size(Ry0.degmat,1)));        % DJ, 01/16/2025
Hyx = expand_polynomial_quadratic(Ryx,Pmat_y(:,nnZ_y(2)+1:nnZ_y(3)),Pmat_th(:,nnZ_th(3)+1:nnZ_th(4)));
Hyy = expand_polynomial_quadratic(Ryy,Pmat_y(:,nnZ_y(3)+1:nnZ_y(4)),Pmat_nu(:,nnZ_nu(3)+1:nnZ_nu(4)));
Hy2 = expand_polynomial_quadratic(Ry2,Pmat_y(:,nnZ_y(4)+1:nnZ_y(5)),Pmat_tn(:,nnZ_tn(3)+1:nnZ_tn(4)));

H20 = expand_polynomial_quadratic(R20,Pmat_xy(:,nnZ_xy(1)+1:nnZ_xy(2)),ones(1,size(R20.degmat,1)));     % DJ, 01/16/2025
H2x = expand_polynomial_quadratic(R2x,Pmat_xy(:,nnZ_xy(2)+1:nnZ_xy(3)),Pmat_th(:,nnZ_th(4)+1:nnZ_th(5)));
H2y = expand_polynomial_quadratic(R2y,Pmat_xy(:,nnZ_xy(3)+1:nnZ_xy(4)),Pmat_nu(:,nnZ_nu(4)+1:nnZ_nu(5)));
H22 = expand_polynomial_quadratic(R22,Pmat_xy(:,nnZ_xy(4)+1:nnZ_xy(5)),Pmat_tn(:,nnZ_tn(4)+1:nnZ_tn(5)));



% % Collect the outputs such that Pop = ZL_op*H*ZR_op;
% Combine coefficients into a single matrix
H = [zeros(nr(1),nc(1)), H0x, H0y, H02;
     Hx0,                Hxx, Hxy, Hx2;
     Hy0,                Hyx, Hyy, Hy2;
     H20,                H2x, H2y, H22];

% Declare multiplier operator ZL_op
ZL_op = opvar2d([nr,[size(H,1);0;0;0]]);
ZL_op.var1 = Pop.var1;      ZL_op.var2 = Pop.var2;  ZL_op.I = Pop.I;
ZL_pars = blkdiag(eye(nr(1)),Zx0',Zy0',Z20');
nnr = cumsum([0;nr]);
ZL_op.R00 = ZL_pars(nnr(1)+1:nnr(2),:);
ZL_op.Rx0 = ZL_pars(nnr(2)+1:nnr(3),:);
ZL_op.Ry0 = ZL_pars(nnr(3)+1:nnr(4),:);
ZL_op.R20 = ZL_pars(nnr(4)+1:nnr(5),:);

% Declare integral operator ZR_op
ZR_op = opvar2d([[size(H,1);0;0;0],nc]);    
ZR_op.var1 = Pop.var1;      ZR_op.var2 = Pop.var2;  ZR_op.I = Pop.I;
ZR_pars = blkdiag(eye(nc(1)),Z0x,Z0y,Z02);
nnc = cumsum([0;nc]);
ZR_op.R00 = ZR_pars(:,nnc(1)+1:nnc(2));
ZR_op.R0x = ZR_pars(:,nnc(2)+1:nnc(3));
ZR_op.R0y = ZR_pars(:,nnc(3)+1:nnc(4));
ZR_op.R02 = ZR_pars(:,nnc(4)+1:nnc(5));

end


%%
function P_dmat = extract_degs(P,varnames)
% Extract degree matrix of the 'polynomial' class object 'P' with respect
% to the variables in "varnames".

nvars = length(varnames);
nZ = size(P.degmat,1);

P_dmat = sparse(nZ,nvars);
[is_var,var_idx] = ismember(varnames,P.varname);
for j=1:nvars
    if is_var(j)
        P_dmat(:,j) = P.degmat(:,var_idx(j));
    end
end

end

%%
function C = expand_polynomial_quadratic(R,Pmat_L,Pmat_R)
% Decompose the 'polynomial' class object "R" in the form
% R = kron(eye(nr),Z_L)' * C * kron(eye(nc),Z_R);
% Here, if "R.degmat" defines a monomial vector "Z_F", we require
% Pmat_L'*Z_L = Z_F and Pmat_R'*Z_R = Z_F;

nZ_F = size(R.degmat,1);
nZ_L = size(Pmat_L,1);      nZ_R = size(Pmat_R,1);
[nr,nc] = size(R);
C = sparse(nr*nZ_L,nc*nZ_R);
for j=1:nZ_F
    % Check which monomial in Z_L and Z_R correspond to monomial Z_F(j)
    idx_L = find(Pmat_L(:,j));
    idx_R = find(Pmat_R(:,j));

    % Determine which rows and columns in C correspond to these monomials
    r_idcs = idx_L:nZ_L:nr*nZ_L;
    c_idcs = idx_R:nZ_R:nc*nZ_R;

    % Set the appropriate elements of C
    C(r_idcs,c_idcs) = C(r_idcs,c_idcs) +reshape(R.C(j,:),[nr,nc]);

end

end