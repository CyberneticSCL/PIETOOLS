function [Phat] = inv_opvar2d(P)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% [Pinv] = inv(P) takes in one operator, 
% P: R^m0 x L2^mx x L2^my to R^n0 x L2^nx x L2^ny
% It returns the (left) inverse Pinv: R^n0 x L2^nx x L2^ny to R^m0 x L2^mx x L2^my
% Version 1.0
% Date: 02/04/21
% 
% INPUT
% P: opvar2d class object
%
% OUTPUT
% Pinv: the matlab structure such that subsequent application Pinv*P*y for
% an arbitrary polynomial y:R^n0 x L2^nx(s1) x L2^ny(s2) returns y.
%
% NOTES:
% For support, contact M. Peet, Arizona State University at mpeet@asu.edu
% or D. Jagt at djagt@asu.edu

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PIETools - inv
%
% Copyright (C)2021  M. Peet, D. Jagt
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
% Initial coding DJ - 02_04_2021  

if isempty(P)
    Phat = P;
    return
end

% Collect necessary polynomials
R00 = double(P.R00);

Rxx0 = polynomial(P.Rxx{1,1});
Rxx1 = polynomial(P.Rxx{2,1});
Rxx2 = polynomial(P.Rxx{3,1});

Ryy0 = polynomial(P.Ryy{1,1});
Ryy1 = polynomial(P.Ryy{1,2});
Ryy2 = polynomial(P.Ryy{1,3});

a = P.I(:,1);
b = P.I(:,2);
ds = P.var1;
dt = P.var2;

m0 = P.dim(1,2);

% Check whether inversion is supported
tol = 1e-15;

if any(size(R00)==0)
    error('P.R00 matrix is empty; inversion not supported')
end
sgm = svd(R00);
if sgm(end)/sgm(1) <= tol
    error('P.R00 matrix is (close to) singular; operator is not invertible')
elseif length(sgm) < size(R00,2)
    error('P.R00 matrix is of insufficient column rank; operator not (left-)invertible')
end

if any(P.dim(4,:)>0)
    error('Operator maps from or to 2-dimensional functions; inversion not supported')
end

Rxxdiff = polynomial(Rxx1-Rxx2);
if any(any(Rxxdiff.coefficient >= tol))
    error('P.Rxx{2,1} does not equal P.Rxx{3,1}; inversion not supported')
else
    Rxx1 = Rxx2;
end

Ryydiff = polynomial(Ryy1-Ryy2);
if any(any(Ryydiff.coefficient >= tol))
    error('P.Ryy{1,2} does not equal P.Ryy{1,3}; inversion not supported')
else
    Ryy1 = Ryy2;
end

Rxx0_degmax = max(Rxx0.degmat);
Ryy0_degmax = max(Ryy0.degmat);
if any(Rxx0_degmax~=0) || any(Ryy0_degmax~=0)
    error('P.Rxx{1,1} or P.Ryy{1,1} is not constant; inversion not supported')
else
    Rxx0 = double(Rxx0);
    sgm = svd(Rxx0);
    if sgm(end)/sgm(1) <= tol
        error('P.Rxx{1,1} matrix is (close to) singular; operator is not invertible')
    elseif length(sgm) < size(Rxx0,2)
        error('P.Rxx{1,1} matrix is of insufficient column rank; operator not (left-)invertible')
    end
    
    Ryy0 = double(Ryy0);
    sgm = svd(Ryy0);
    if sgm(end)/sgm(1) <= tol
        error('P.Ryy{1,1} matrix is (close to) singular; operator is not invertible')
    elseif length(sgm) < size(Ryy0,2)
        error('P.Ryy{1,1} matrix is of insufficient column rank; operator not (left-)invertible')
    end
end


% Perform the separation into constant and polynomial parts
P.Rxx = Rxx1;   P.Ryy = Ryy1;
[Z,H] = poly_separate(P);

Z0x = Z.Z0x;    Z0y = Z.Z0y;
Zx0 = Z.Zx0;    Zy0 = Z.Zy0;

                H0x = H.H0x;  H0y = H.H0y;
Hx0 = H.Hx0;    Gxx = H.Gxx;  Gxy = H.Gxy;
Hy0 = H.Hy0;    Gyx = H.Gyx;  Gyy = H.Gyy;

q0x = size(Z0x,1);  q0y = size(Z0y,1);


% Construct the inverse
Ryy0Hat = pinv(Ryy0);
Rxx0Hat = pinv(Rxx0);
R00inv = pinv(R00);

Kxx = double(int(Z0x*Rxx0Hat*Zx0',ds(1),a(1),b(1)));
Kyy = double(int(Z0y*Ryy0Hat*Zy0',ds(2),a(2),b(2)));

Pixx = Gxx - Hx0*R00inv*H0x;
Pixy = Gxy - Hx0*R00inv*H0y;
Piyx = Gyx - Hy0*R00inv*H0x;
Piyy = Gyy - Hy0*R00inv*H0y;

M1 = pinv(eye(q0x) + Kxx*Pixx - (Kxx*Pixy/(eye(q0y) + Kyy*Piyy)) * Kyy*Piyx);
GyxHat = -(Piyx - (Piyy/(eye(q0y) + Kyy*Piyy)) * Kyy*Piyx) * M1;
GyyHat = -(Piyy + GyxHat*Kxx*Pixy)/(eye(q0y) + Kyy*Piyy);
Hy0Hat = -(Hy0 + GyyHat*Kyy*Hy0 + GyxHat*Kxx*Hx0) * R00inv;

M2 = pinv(eye(q0y) + Kyy*Piyy - (Kyy*Piyx/(eye(q0x) + Kxx*Pixx)) * Kxx*Pixy);
GxyHat = -(Pixy - (Pixx/(eye(q0x) + Kxx*Pixx)) * Kxx*Pixy) * M2;
GxxHat = -(Pixx + GxyHat*Kyy*Piyx)/(eye(q0x) + Kxx*Pixx);
Hx0Hat = -(Hx0 + GxxHat*Kxx*Hx0 + GxyHat*Kyy*Hy0) * R00inv;

H0yHat = -(R00inv*H0y - (R00inv*H0x/(eye(q0x) + Kxx*Pixx)) * Kxx*Pixy) * M2;
H0xHat = -(R00inv*H0x + H0yHat*Kyy*Piyx)/(eye(q0x) + Kxx*Pixx);
R00Hat = (eye(m0) - H0xHat*Kxx*Hx0 - H0yHat*Kyy*Hy0) * R00inv;
        
Phat = opvar2d();
Phat.I = P.I;
%Phat.dim = P.dim;

Phat.R00 = R00Hat;
Phat.R0x = H0xHat*Z0x*Rxx0Hat;
Phat.R0y = H0yHat*Z0y*Ryy0Hat;

Phat.Rx0 = Rxx0Hat*Zx0'*Hx0Hat;
Phat.Rxx{1,1} = Rxx0Hat;
Phat.Rxx{2,1} = Rxx0Hat*Zx0'*GxxHat*subs(Z0x*Rxx0Hat,ds(1),dt(1));
Phat.Rxx{3,1} = Rxx0Hat*Zx0'*GxxHat*subs(Z0x*Rxx0Hat,ds(1),dt(1));
Phat.Rxy = Rxx0Hat*Zx0'*GxyHat*Z0y*Ryy0Hat;

Phat.Ry0 = Ryy0Hat*Zy0'*Hy0Hat;
Phat.Ryx = Ryy0Hat*Zy0'*GyxHat*Z0x*Rxx0Hat;
Phat.Ryy{1,1} = Ryy0Hat;
Phat.Ryy{1,2} = Ryy0Hat*Zy0'*GyyHat*subs(Z0y*Ryy0Hat,ds(2),dt(2));
Phat.Ryy{1,3} = Ryy0Hat*Zy0'*GyyHat*subs(Z0y*Ryy0Hat,ds(2),dt(2));

Phat.dim = Phat.dim;

end



function [Z,H]=poly_separate(Pxy)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This function splits the matrix-valued functions R0x, R0y, Rx0, Ry0, Rxy, Ryx, Rxx, Ryy as:
%                   R0x = H0x*Z0x(s1);              R0y = H0y*Z0y(s2);
% Rx0 = Zx0'*Hx0;   Rxx = Zx0'(s1)*Gxx*Z0x(t1);     Rxy = Zx0'(s1)*Gxy*Z0y(s2)
% Ry0 = Zy0'*Hy0;   Ryx = Zy0'(s2)*Gyx*Z0x(s1);     Ryy = Zy0'(s2)*Gyy*Z0y(t2)

% Here matrices Z.. each comprise a basis of vector-valued polynomials in
% (s1,s2,t1,t2)
% Matrices H0x, H0y, Hx0, Hy0, Gxx, Gxy, Gyx, and Gyy are constant
% Version 1.0
% Date: 02/04/21
% 
% Inputs:
% Pxy: A struct containing the polynomial functions of the 0xy subblock of
% a 2d PI operator
% 
% Output:
% Separated matrices

% Collext the inputs
                            R0x = polynomial(Pxy.R0x);  R0y = polynomial(Pxy.R0y);
Rx0 = polynomial(Pxy.Rx0);  Rxx = polynomial(Pxy.Rxx);  Rxy = polynomial(Pxy.Rxy);
Ry0 = polynomial(Pxy.Ry0);  Ryx = polynomial(Pxy.Ryx);  Ryy = polynomial(Pxy.Ryy);

ds = Pxy.var1;      dt = Pxy.var2;

% Compute the necessary degrees of monomial matrices Z..
if size(R0x.degmat,2)==0
    s1_degR0x = 0;
else
    s1_degR0x = R0x.degmat;         
end
if size(R0y.degmat,2)==0
    s2_degR0y = 0;
else
    s2_degR0y = R0y.degmat;
end
if size(Rx0.degmat,2)==0
    s1_degRx0 = 0;
else
    s1_degRx0 = Rx0.degmat;
end
if size(Ry0.degmat,2)==0
    s2_degRy0 = 0;
else
    s2_degRy0 = Ry0.degmat;
end

if size(Rxx.degmat,2)==0
    s1_degRxx = 0;                  t1_degRxx = 0;
else
    s1_degRxx = Rxx.degmat(:,1);    t1_degRxx = Rxx.degmat(:,2); 
end
if size(Rxy.degmat,2)==0
    s1_degRxy = 0;                  s2_degRxy = 0; 
else
    s1_degRxy = Rxy.degmat(:,1);    s2_degRxy = Rxy.degmat(:,2); 
end
if size(Ryx.degmat,2)==0
    s1_degRyx = 0;                  s2_degRyx = 0; 
else
    s1_degRyx = Ryx.degmat(:,1);    s2_degRyx = Ryx.degmat(:,2); 
end
if size(Ryy.degmat,2)==0
    s2_degRyy = 0;                  t2_degRyy = 0; 
else
    s2_degRyy = Ryy.degmat(:,1);    t2_degRyy = Ryy.degmat(:,2); 
end

deg_min0x = min([min(s1_degR0x),min(t1_degRxx),min(s1_degRyx)]);
deg_max0x = max([max(s1_degR0x),max(t1_degRxx),max(s1_degRyx)]);

deg_min0y = min([min(s2_degR0y),min(s2_degRxy),min(t2_degRyy)]);
deg_max0y = max([max(s2_degR0y),max(s2_degRxy),max(t2_degRyy)]);

deg_minx0 = min([min(s1_degRx0),min(s1_degRxx),min(s1_degRxy)]);
deg_maxx0 = max([max(s1_degRx0),max(s1_degRxx),max(s1_degRxy)]);

deg_miny0 = min([min(s2_degRy0),min(s2_degRyx),min(s2_degRyy)]);
deg_maxy0 = max([max(s2_degRy0),max(s2_degRyx),max(s2_degRyy)]);



% Construct vectors of monomials
degvec0x = deg_min0x:deg_max0x;
ndeg0x = length(degvec0x);
degvec0y = deg_min0y:deg_max0y;
ndeg0y = length(degvec0y);
degvecx0 = deg_minx0:deg_maxx0;
ndegx0 = length(degvecx0);
degvecy0 = deg_miny0:deg_maxy0;
ndegy0 = length(degvecy0);

z0x = monomials(ds(1),degvec0x);
z0y = monomials(ds(2),degvec0y);
zx0 = monomials(ds(1),degvecx0);
zy0 = monomials(ds(2),degvecy0);



% Build matrices Z0x, Z0y, Zx0 and Zy0 
nx = size(Rxx,1);
mx = size(Rxx,2);
ny = size(Ryy,1);
my = size(Ryy,2);

Z0x = polynomial(zeros(length(z0x)*mx,mx));
for j=0:mx-1
    Z0x(j*length(z0x)+1:(j+1)*length(z0x),j+1) = z0x;
end

Z0y = polynomial(zeros(length(z0y)*my,my));
for j=0:my-1
    Z0y(j*length(z0y)+1:(j+1)*length(z0y),j+1) = z0y;
end

Zx0 = polynomial(zeros(length(zx0)*nx,nx));
for j=0:nx-1
    Zx0(j*length(zx0)+1:(j+1)*length(zx0),j+1) = zx0;
end

Zy0 = polynomial(zeros(length(zy0)*ny,ny));
for j=0:ny-1
    Zy0(j*length(zx0)+1:(j+1)*length(zy0),j+1) = zy0;
end



% Construct the constant coefficient matrices
H0x = zeros(size(R0x,1),mx*ndeg0x);
for d=1:size(R0x.coeff,1)
    indx = s1_degR0x(d) - deg_min0x + 1;
    coeffR0x = reshape(R0x.coeff(d,:),size(R0x));
    for j=0:mx-1
        H0x(:,indx+j*ndeg0x) = H0x(:,indx+j*ndeg0x) + coeffR0x(:,j+1);
    end
end

H0y = zeros(size(R0y,1),my*ndeg0y);
for d=1:size(R0y.coeff,1)
    indx = s2_degR0y(d) - deg_min0y + 1;
    coeffR0y = reshape(R0y.coeff(d,:),size(R0y));
    for j=0:my-1
        H0y(:,indx+j*ndeg0y) = H0y(:,indx+j*ndeg0y) + coeffR0y(:,j+1);
    end
end

Hx0 = zeros(nx*ndegx0,size(Rx0,2));
for d=1:size(Rx0.coeff,1)
    indx = s1_degRx0(d) - deg_minx0 + 1;
    coeffRx0 = reshape(Rx0.coeff(d,:),size(Rx0));
    for i=0:nx-1
        Hx0(indx+i*ndegx0,:) = Hx0(indx+i*ndegx0,:) + coeffRx0(i+1,:);
    end
end

Hy0 = zeros(ny*ndegy0,size(Ry0,2));
for d=1:size(Ry0.coeff,1)
    indx = s2_degRy0(d) - deg_miny0 + 1;
    coeffRy0 = reshape(Ry0.coeff(d,:),size(Ry0));
    for i=0:ny-1
        Hy0(indx+i*ndegy0,:) = Hy0(indx+i*ndegy0,:) + coeffRy0(i+1,:);
    end
end



Gxx = zeros(nx*ndegx0,mx*ndeg0x);
for d=1:size(Rxx.coeff,1)
    indx_s = s1_degRxx(d) - deg_minx0 + 1;
    indx_t = t1_degRxx(d) - deg_min0x + 1;
    coeffRxx = reshape(Rxx.coeff(d,:),size(Rxx));
    for i=0:nx-1
        for j=0:mx-1
            Gxx(indx_s+i*ndegx0,indx_t+j*ndeg0x) = Gxx(indx_s+i*ndegx0,indx_t+j*ndeg0x) + coeffRxx(i+1,j+1);
        end
    end
end

Gxy = zeros(nx*ndegx0,my*ndeg0y);
for d=1:size(Rxy.coeff,1)
    indx_s = s1_degRxy(d) - deg_minx0 + 1;
    indx_t = s2_degRxy(d) - deg_min0y + 1;
    coeffRxy = reshape(Rxy.coeff(d,:),size(Rxy));
    for i=0:nx-1
        for j=0:my-1
            Gxy(indx_s+i*ndegx0,indx_t+j*ndeg0y) = Gxy(indx_s+i*ndegx0,indx_t+j*ndeg0y) + coeffRxy(i+1,j+1);
        end
    end
end

Gyx = zeros(ny*ndegy0,mx*ndeg0x);
for d=1:size(Ryx.coeff,1)
    indx_s = s2_degRyx(d) - deg_miny0 + 1;
    indx_t = s1_degRyx(d) - deg_min0x + 1;
    coeffRyx = reshape(Ryx.coeff(d,:),size(Ryx));
    for i=0:ny-1
        for j=0:mx-1
            Gyx(indx_s+i*ndegy0,indx_t+j*ndeg0x) = Gyx(indx_s+i*ndegy0,indx_t+j*ndeg0x) + coeffRyx(i+1,j+1);
        end
    end
end

Gyy = zeros(ny*ndegy0,my*ndeg0y);
for d=1:size(Ryy.coeff,1)
    indx_s = s2_degRyy(d) - deg_miny0 + 1;
    indx_t = t2_degRyy(d) - deg_min0y + 1;
    coeffRyy = reshape(Ryy.coeff(d,:),size(Ryy));
    for i=0:ny-1
        for j=0:my-1
            Gyy(indx_s+i*ndegy0,indx_t+j*ndeg0y) = Gyy(indx_s+i*ndegy0,indx_t+j*ndeg0y) + coeffRyy(i+1,j+1);
        end
    end
end


% Collect the outputs
Z = struct();
Z.Z0x = Z0x;
Z.Z0y = Z0y;
Z.Zx0 = Zx0;
Z.Zy0 = Zy0;

H = struct();
                H.H0x = H0x;    H.H0y = H0y;
H.Hx0 = Hx0;    H.Gxx = Gxx;    H.Gxy = Gxy;
H.Hy0 = Hy0;    H.Gyx = Gyx;    H.Gyy = Gyy;

end
