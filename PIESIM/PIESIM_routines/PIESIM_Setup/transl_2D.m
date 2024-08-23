%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% transl_2D.m     PIETOOLS 2024
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% If you modify this code, document all changes carefully and include date
% authorship, and a brief description of modifications
%
% Initial coding YP - 4_16_2024

function T = transl_2D(T,I)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% T = transl_2D(T,I) takes in a full 2D PI operator T that acts on
% functions on the original domain and changes it to a 2D PI operator that
% acts on functions on the domain I.
%
% INPUT:
%
% T: a 2D PI operator
% I: new domain
% 
% NOTES:
% For support, contact M. Peet, Arizona State University at mpeet@asu.edu,
% or Y. Peet at ypeet@asu.edu
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PIETools - transl_2D
%
% Copyright (C)2021  M. Peet, Y. Peet,
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
% Initial coding YP  - 4_16_2024

I_init = T.I;
ax = I_init(1,1);
bx = I_init(1,2);
ay = I_init(2,1);
by = I_init(2,2);

cx = I(1,1);
dx = I(1,2);
cy = I(2,1);
dy = I(2,2);

% First column operators
% R00 - scalar- no need to translate


if ~isempty(T.Rx0.varname)
T.Rx0=transl_1D(T.Rx0,ax,bx,cx,dx,'x',1);
end
if ~isempty(T.Ry0.varname)
T.Ry0=transl_1D(T.Ry0,ay,by,cy,dy,'y',1);
end


if ~isempty(T.R20.varname)
T.R20=transl_1D(T.R20,ay,by,cy,dy,'y',1);
T.R20=transl_1D(T.R20,ax,bx,cx,dx,'x',1);
end

% Second column operators

if ~isempty(T.R0x.varname)
T.R0x=transl_1D(T.R0x,ax,bx,cx,dx,'x',2);
else
T.R0x=((bx-ax)/(dx-cx))*T.R0x;
end

T.Rxx=transl_3PI(T.Rxx,ax,bx,cx,dx,'x');

if ~isempty(T.Ryx.varname)
T.Ryx=transl_1D(T.Ryx,ay,by,cy,dy,'y',1);
T.Ryx=transl_1D(T.Ryx,ax,bx,cx,dx,'x',2);
else
   T.Ryx=((bx-ax)/(dx-cx))*T.Ryx; 
end

T.R2x=transl_3PI(T.R2x,ax,bx,cx,dx,'x');

for k=1:3
if ~isempty(T.R2x{k}.varname)
T.R2x{k}=transl_1D(T.R2x{k},ay,by,cy,dy,'y',1);
end
end

% Third column operators


if ~isempty(T.R0y.varname)
T.R0y=transl_1D(T.R0y,ay,by,cy,dy,'y',2);
else
T.R0y=((by-ay)/(dy-cy))*T.R0y; 
end

if ~isempty(T.Rxy.varname)
T.Rxy=transl_1D(T.Rxy,ax,bx,cx,dx,'x',1);
T.Rxy=transl_1D(T.Rxy,ay,by,cy,dy,'y',2);
else
T.Rxy=((by-ay)/(dy-cy))*T.Rxy; 
end

T.Ryy=transl_3PI(T.Ryy,ay,by,cy,dy,'y');

T.R2y=transl_3PI(T.R2y,ay,by,cy,dy,'y');

for k=1:3
if ~isempty(T.R2y{k}.varname)
T.R2y{k}=transl_1D(T.R2y{k},ax,bx,cx,dx,'x',1);
end
end

% Fourth column operators

if~isempty(T.R02.varname)
T.R02=transl_1D(T.R02,ay,by,cy,dy,'y',2);
T.R02=transl_1D(T.R02,ax,bx,cx,dx,'x',2);
else
T.R02=((bx-ax)/(dx-cx))*((by-ay)/(dy-cy))*T.R02; 
end

T.Rx2=transl_3PI(T.Rx2,ax,bx,cx,dx,'x');
for k=1:3
if ~isempty(T.Rx2{k}.varname)
T.Rx2{k}=transl_1D(T.Rx2{k},ay,by,cy,dy,'y',2);
else
T.Rx2{k}=((by-ay)/(dy-cy))*T.Rx2{k};     
end
end

T.Ry2=transl_3PI(T.Ry2,ay,by,cy,dy,'y');
for k=1:3
if ~isempty(T.Ry2{k}.varname)
T.Ry2{k}=transl_1D(T.Ry2{k},ax,bx,cx,dx,'x',2);
else
T.Ry2{k}=((bx-ax)/(dx-cx))*T.Ry2{k};  
end
end

if ~isempty(T.R22{1,1}.varname)
T.R22{1,1}=transl_1D(T.R22{1,1},ay,by,cy,dy,'y',1);
T.R22{1,1}=transl_1D(R22{1,1},ax,bx,cx,dx,'x',1);
end

for k=2:3
if ~isempty(T.R22{1,k}.varname)
T.R22{1,k}=transl_1D(T.R22{1,k},ay,by,cy,dy,'y',3);
T.R22{1,k} = transl_1D(T.R22{1,k},ax,bx,cx,dx,'x',1);
else
T.R22{1,k}=((by-ay)/(dy-cy))*T.R22{1,k}; 
end

if ~isempty(T.R22{k,1}.varname)
T.R22{k,1}=transl_1D(T.R22{k,1},ax,bx,cx,dx,'x',3);
T.R22{k,1} = transl_1D(T.R22{k,1},ay,by,cy,dy,'y',1);
else
T.R22{k,1}=((bx-ax)/(dx-cx))*T.R22{k,1}; 
end

end

for n=2:3
    for k=2:3
    if ~isempty(T.R22{k,n}.varname)
        T.R22{k,n}=transl_1D(T.R22{k,n},ay,by,cy,dy,'y',3);
        T.R22{k,n}=transl_1D(T.R22{k,n},ax,bx,cx,dx,'x',3);
    else
    T.R22{k,n}=((bx-ax)/(dx-cx))*((by-ay)/(dy-cy))*T.R22{k,n}; 
    end

    end
end

