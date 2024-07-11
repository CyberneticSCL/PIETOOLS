function T = transl_2D(T,dom,I)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% T = transl(T,dom,I) takes in a full 2D PI operator T that acts on
% functions on the domain 'dom' and changes it to a 2D PI operator that
% acts on functions on the domain I.
%
% INPUT:
%
% T: a 2D PI operator
% dom: original domain
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

I_init = dom;
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

T.Rx0=transl_1D(T.Rx0,ax,bx,cx,dx,'x',1);
T.Ry0=transl_1D(T.Ry0,ay,by,cy,dy,'y',1);

T_inter=transl_1D(T.R20,ay,by,cy,dy,'y',1);
T.R20=transl_1D(T_inter,ax,bx,cx,dx,'x',1);

% Second column operators

T.R0x=transl_1D(T.R0x,ax,bx,cx,dx,'x',2);
T.Rxx=transl_3PI(T.Rxx,ax,bx,cx,dx,'x');

if ~isempty(T.Ryx.varname)
T_inter=transl_1D(T.Ryx,ay,by,cy,dy,'y',1);
T.Ryx=transl_1D(T_inter,ax,bx,cx,dx,'x',2);
else
   T.Ryx=((bx-ax)/(dx-cx))*T.Ryx; 
end

T_inter=transl_3PI(T.R2x,ax,bx,cx,dx,'x');
for k=1:3
if ~isempty(T.R2x{k}.varname)
T.R2x{k}=transl_1D(T_inter{k},ay,by,cy,dy,'y',1);
end
end

% Third column operators

T.R0y=transl_1D(T.R0y,ay,by,cy,dy,'y',2);


if ~isempty(T.Rxy.varname)
T_inter=transl_1D(T.Rxy,ax,bx,cx,dx,'x',1);
T.Rxy=transl_1D(T_inter,ay,by,cy,dy,'y',2);
else
   T.Rxy=((by-ay)/(dy-cy))*T.Rxy; 
end

T.Ryy=transl_3PI(T.Ryy,ay,by,cy,dy,'y');

T_inter=transl_3PI(T.R2y,ay,by,cy,dy,'y');
for k=1:3
if ~isempty(T.R2y{k}.varname)
T.R2y{k}=transl_1D(T_inter{k},ax,bx,cx,dx,'x',1);
end
end

% Fourth column operators

if~isempty(T.R02.varname)
T_inter=transl_1D(T.R02,ay,by,cy,dy,'y',2);
T.R02=transl_1D(T_inter,ax,bx,cx,dx,'x',2);
end

T_inter=transl_3PI(T.Rx2,ax,bx,cx,dx,'x');
for k=1:3
if ~isempty(T.Rx2{k}.varname)
T.Rx2{k}=transl_1D(T_inter{k},ay,by,cy,dy,'y',2);
else
T.Rx2{k}=((by-ay)/(dy-cy))*T.Rx2{k};     
end
end

T_inter=transl_3PI(T.Ry2,ay,by,cy,dy,'y');
for k=1:3
if ~isempty(T.Ry2{k}.varname)
T.Ry2{k}=transl_1D(T_inter{k},ax,bx,cx,dx,'x',2);
else
T.Ry2{k}=((bx-ax)/(dx-cx))*T.Ry2{k};  
end
end

T_inter=transl_1D(T.R22{1,1},ay,by,cy,dy,'y',1);
T.R22{1,1}=transl_1D(T_inter,ax,bx,cx,dx,'x',1);

for k=2:3
T_inter=transl_1D(T.R22{1,k},ay,by,cy,dy,'y',3);
T.R22{1,k} = transl_1D(T_inter,ax,bx,cx,dx,'x',1);
T_inter=transl_1D(T.R22{k,1},ax,bx,cx,dx,'x',3);
T.R22{k,1} = transl_1D(T_inter,ay,by,cy,dy,'y',1);
end

for n=2:3
    for k=2:3
        T_inter=transl_1D(T.R22{k,n},ay,by,cy,dy,'y',3);
        T.R22{k,n}=transl_1D(T_inter,ax,bx,cx,dx,'x',3);
    end
end

