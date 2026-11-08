function logval = eq(P1,P2,tol)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% logval = eq(P1,P2,tol) tests equality of two opvar2ds P1 and P2 with tolerance tol
% Date: 02/04/21
% Version: 1.0
% 
% INPUT
% P1, P2: opvar2d class objects
% tol: acceptable tolerance value. If max(P1-P2)<tol, then P1=P2
% 
% OUTPUT
% logval: returns 1 if the objects are equal, 0 if not equal
%
% NOTES:
% For support, contact M. Peet, Arizona State University at mpeet@asu.edu
% or D. Jagt at djagt@asu.edu

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PIETools - eq
%
% Copyright (C)2021  M. Peet, S. Shivakumar, D. Jagt
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
% Initial coding DJ - 02_04_2021 
%   ^ Based heavily on "@opvar"-eq code by SS ^
% DJ, 12/11/2024: Bug fix Pop==0, make sure independent variables match;

if nargin<3
    tol=1e-14;
end    

if ~isa(P1,'opvar2d')&&(P1==0) 
    opvar2d P1; P1.I = P2.I; P1.dim = P2.dim; P1.var1 = P2.var1; P1.var2 = P2.var2;     % DJ, 12/11/2024
elseif ~isa(P2,'opvar2d')&&(P2==0)
    opvar2d P2; P2.I = P1.I; P2.dim = P1.dim; P2.var1 = P1.var1; P2.var2 = P1.var2;     % DJ, 12/11/2024
elseif ~isa(P1,'opvar2d')|| ~isa(P2,'opvar2d')
    error('To check equality either both values must be opvar2d objects, or one of them has to be zero');
end


if any(P1.dim~=P2.dim)
    disp('Opvars have different dimensions and hence cannot be equal');
    logval=0;
    return
end

diff = P1-P2; 

diff.R00 = double(diff.R00);
diff.R0x = polynomial(diff.R0x);
diff.R0y = polynomial(diff.R0y);
diff.R02 = polynomial(diff.R02);
diff.Rx0 = polynomial(diff.Rx0);
diff.Rxy = polynomial(diff.Rxy);
diff.Ry0 = polynomial(diff.Ry0);
diff.Ryx = polynomial(diff.Ryx);
diff.R20 = polynomial(diff.R20);

diff.R00(abs(diff.R00)<tol)=0;
diff.R0x.coefficient(abs(diff.R0x.coefficient)<tol)=0;
diff.R0y.coefficient(abs(diff.R0y.coefficient)<tol)=0;
diff.R02.coefficient(abs(diff.R02.coefficient)<tol)=0;
diff.Rx0.coefficient(abs(diff.Rx0.coefficient)<tol)=0;
diff.Rxy.coefficient(abs(diff.Rxy.coefficient)<tol)=0;
diff.Ry0.coefficient(abs(diff.Ry0.coefficient)<tol)=0;
diff.Ryx.coefficient(abs(diff.Ryx.coefficient)<tol)=0;
diff.R20.coefficient(abs(diff.R20.coefficient)<tol)=0;

for i=1:3
    diff.Rxx{i,1} = polynomial(diff.Rxx{i,1});
    diff.Rxx{i,1}.coefficient(abs(diff.Rxx{i,1}.coefficient)<tol)=0;
    diff.Rx2{i,1} = polynomial(diff.Rx2{i,1});
    diff.Rx2{i,1}.coefficient(abs(diff.Rx2{i,1}.coefficient)<tol)=0;
    diff.R2x{i,1} = polynomial(diff.R2x{i,1});
    diff.R2x{i,1}.coefficient(abs(diff.R2x{i,1}.coefficient)<tol)=0;
    
    diff.Ryy{1,i} = polynomial(diff.Ryy{1,i});
    diff.Ryy{1,i}.coefficient(abs(diff.Ryy{1,i}.coefficient)<tol)=0;
    diff.Ry2{1,i} = polynomial(diff.Ry2{1,i});
    diff.Ry2{1,i}.coefficient(abs(diff.Ry2{1,i}.coefficient)<tol)=0;
    diff.R2y{1,i} = polynomial(diff.R2y{1,i});
    diff.R2y{1,i}.coefficient(abs(diff.R2y{1,i}.coefficient)<tol)=0;
    
    for j=1:3
        diff.R22{i,j} = polynomial(diff.R22{i,j});
        diff.R22{i,j}.coefficient(abs(diff.R22{i,j}.coefficient)<tol)=0;
    end
end


try
    Req = 1;
    diff.R00 = double(diff.R00);
    Req = all(all(diff.R00==0)) && Req;
    diff.R0x = double(diff.R0x);
    Req = all(all(diff.R0x==0)) && Req;
    diff.R0y = double(diff.R0y);
    Req = all(all(diff.R0y==0)) && Req;
    diff.R02 = double(diff.R02);
    Req = all(all(diff.R02==0)) && Req;
    diff.Rx0 = double(diff.Rx0);
    Req = all(all(diff.Rx0==0)) && Req;
    diff.Rxy = double(diff.Rxy);
    Req = all(all(diff.Rxy==0)) && Req;
    diff.Ry0 = double(diff.Ry0);
    Req = all(all(diff.Ry0==0)) && Req;
    diff.Ryx = double(diff.Ryx);
    Req = all(all(diff.Ryx==0)) && Req;
    diff.R20 = double(diff.R20);
    Req = all(all(diff.R20==0)) && Req;

    for i=1:3
        diff.Rxx{i,1} = double(diff.Rxx{i,1});
        Req = all(all(diff.Rxx{i,1}==0)) && Req;
        diff.Rx2{i,1} = double(diff.Rx2{i,1});
        Req = all(all(diff.Rx2{i,1}==0)) && Req;
        diff.R2x{i,1} = double(diff.R2x{i,1});
        Req = all(all(diff.R2x{i,1}==0)) && Req;

        diff.Ryy{1,i} = double(diff.Ryy{1,i});
        Req = all(all(diff.Ryy{1,i}==0)) && Req;
        diff.Ry2{1,i} = double(diff.Ry2{1,i});
        Req = all(all(diff.Ry2{1,i}==0)) && Req;
        diff.R2y{1,i} = double(diff.R2y{1,i});
        Req = all(all(diff.R2y{1,i}==0)) && Req;

        for j=1:3
            diff.R22{i,j} = double(diff.R22{i,j});
            Req = all(all(diff.R22{i,j}==0)) && Req;
        end
    end
    
    if Req
        logval=1;
    else
        logval=0;
    end
catch 
    logval=0;
end
end
