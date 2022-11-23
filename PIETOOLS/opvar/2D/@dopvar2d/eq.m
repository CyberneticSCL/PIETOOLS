function logval = eq(P1,P2,tol)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% logval = eq(P1,P2,tol) tests equality of two opvar2ds P1 and P2 with tolerance tol
% Date: 07/12/21
% Version: 1.0
% 
% INPUT
% P1, P2: dopvar2d class objects
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
% Initial coding DJ - 07_12_2021

if nargin<3
    tol=1e-14;
end    

if ~isa(P1,'dopvar2d')&&(P1==0) 
    dopvar2d P1; P1.I = P2.I; P1.dim = P2.dim;
elseif ~isa(P2,'dopvar2d')&&(P2==0)
    dopvar2d P2; P2.I = P1.I; P2.dim = P1.dim;
elseif ~isa(P1,'dopvar2d')|| ~isa(P2,'dopvar2d')
    error('To check equality either both values must be dopvar2d objects, or one of them has to be zero');
end


if any(P1.dim~=P2.dim)
    disp('Dopvars have different dimensions and hence cannot be equal');
    logval=0;
    return
end

diff = P1-P2; 

% Convert all components to dpvars
diff.R00 = dpvar(diff.R00);
diff.R0x = dpvar(diff.R0x);
diff.R0y = dpvar(diff.R0y);
diff.R02 = dpvar(diff.R02);
diff.Rx0 = dpvar(diff.Rx0);
diff.Rxy = dpvar(diff.Rxy);
diff.Ry0 = dpvar(diff.Ry0);
diff.Ryx = dpvar(diff.Ryx);
diff.R20 = dpvar(diff.R20);

% Set all differences that are sufficiently small to zero
diff.R00.C(abs(diff.R00.C)<tol)=0;
diff.R0x.C(abs(diff.R0x.C)<tol)=0;
diff.R0y.C(abs(diff.R0y.C)<tol)=0;
diff.R02.C(abs(diff.R02.C)<tol)=0;
diff.Rx0.C(abs(diff.Rx0.C)<tol)=0;
diff.Rxy.C(abs(diff.Rxy.C)<tol)=0;
diff.Ry0.C(abs(diff.Ry0.C)<tol)=0;
diff.Ryx.C(abs(diff.Ryx.C)<tol)=0;
diff.R20.C(abs(diff.R20.C)<tol)=0;

for i=1:3
    diff.Rxx{i,1} = dpvar(diff.Rxx{i,1});
    diff.Rxx{i,1}.C(abs(diff.Rxx{i,1}.C)<tol)=0;
    diff.Rx2{i,1} = dpvar(diff.Rx2{i,1});
    diff.Rx2{i,1}.C(abs(diff.Rx2{i,1}.C)<tol)=0;
    diff.R2x{i,1} = dpvar(diff.R2x{i,1});
    diff.R2x{i,1}.C(abs(diff.R2x{i,1}.C)<tol)=0;
    
    diff.Ryy{1,i} = dpvar(diff.Ryy{1,i});
    diff.Ryy{1,i}.C(abs(diff.Ryy{1,i}.C)<tol)=0;
    diff.Ry2{1,i} = dpvar(diff.Ry2{1,i});
    diff.Ry2{1,i}.C(abs(diff.Ry2{1,i}.C)<tol)=0;
    diff.R2y{1,i} = dpvar(diff.R2y{1,i});
    diff.R2y{1,i}.C(abs(diff.R2y{1,i}.C)<tol)=0;
    
    for j=1:3
        diff.R22{i,j} = dpvar(diff.R22{i,j});
        diff.R22{i,j}.C(abs(diff.R22{i,j}.C)<tol)=0;
    end
end

% Check whether all components are zero
try
    Req = 1;
    diff.R00 = double(DPcompress(diff.R00));
    Req = all(all(diff.R00==0)) && Req;
    diff.R0x = double(DPcompress(diff.R0x));
    Req = all(all(diff.R0x==0)) && Req;
    diff.R0y = double(DPcompress(diff.R0y));
    Req = all(all(diff.R0y==0)) && Req;
    diff.R02 = double(DPcompress(diff.R02));
    Req = all(all(diff.R02==0)) && Req;
    diff.Rx0 = double(DPcompress(diff.Rx0));
    Req = all(all(diff.Rx0==0)) && Req;
    diff.Rxy = double(DPcompress(diff.Rxy));
    Req = all(all(diff.Rxy==0)) && Req;
    diff.Ry0 = double(DPcompress(diff.Ry0));
    Req = all(all(diff.Ry0==0)) && Req;
    diff.Ryx = double(DPcompress(diff.Ryx));
    Req = all(all(diff.Ryx==0)) && Req;
    diff.R20 = double(DPcompress(diff.R20));
    Req = all(all(diff.R20==0)) && Req;

    for i=1:3
        diff.Rxx{i,1} = double(DPcompress(diff.Rxx{i,1}));
        Req = all(all(diff.Rxx{i,1}==0)) && Req;
        diff.Rx2{i,1} = double(DPcompress(diff.Rx2{i,1}));
        Req = all(all(diff.Rx2{i,1}==0)) && Req;
        diff.R2x{i,1} = double(DPcompress(diff.R2x{i,1}));
        Req = all(all(diff.R2x{i,1}==0)) && Req;

        diff.Ryy{1,i} = double(DPcompress(diff.Ryy{1,i}));
        Req = all(all(diff.Ryy{1,i}==0)) && Req;
        diff.Ry2{1,i} = double(DPcompress(diff.Ry2{1,i}));
        Req = all(all(diff.Ry2{1,i}==0)) && Req;
        diff.R2y{1,i} = double(DPcompress(diff.R2y{1,i}));
        Req = all(all(diff.R2y{1,i}==0)) && Req;

        for j=1:3
            diff.R22{i,j} = double(DPcompress(diff.R22{i,j}));
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
