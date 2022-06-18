function Psop = op_slice(Pbop,indr,indc)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Psop=op_slice(Pbop,indr,indc) extracts slice from a decision operator 
% P: R^m0 x L2^mx x L2^my x L2^m2 to R^n0 x L2^nx x L2^ny x L2^n2
% Date: 07/12/2021
% Version: 1.0
% 
% INPUT
% Pbop: opvar2d class object to slice
% indr: row indices of the slice
% indc: column indices of the slice
% 
% OUTPUT
% Psop: slice of the dopvar2d object
%  
% NOTES:
% For support, contact M. Peet, Arizona State University at mpeet@asu.edu
% or D. Jagt at djagt@asu.edu

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PIETools - op_slice
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

% Decompose the input opvar
dim_old = Pbop.dim; 
R00old = Pbop.R00;  R0xold = Pbop.R0x;  R0yold = Pbop.R0y;  R02old = Pbop.R02;
Rx0old = Pbop.Rx0;  Rxxold = Pbop.Rxx;  Rxyold = Pbop.Rxy;  Rx2old = Pbop.Rx2;
Ry0old = Pbop.Ry0;  Ryxold = Pbop.Ryx;  Ryyold = Pbop.Ryy;  Ry2old = Pbop.Ry2;
R20old = Pbop.R20;  R2xold = Pbop.R2x;  R2yold = Pbop.R2y;  R22old = Pbop.R22;


% Determine sliced indices for each subcomponent
n0_old = dim_old(1,1);      m0_old = dim_old(1,2);
nx_old = dim_old(2,1);      mx_old = dim_old(2,2);
ny_old = dim_old(3,1);      my_old = dim_old(3,2);
n2_old = dim_old(4,1);      m2_old = dim_old(4,2);

n0_new = 0;                 m0_new = 0;
nx_new = 0;                 mx_new = 0;
ny_new = 0;                 my_new = 0;
n2_new = 0;                 m2_new = 0;

indr0 = [];                 indc0 = [];
indrx = [];                 indcx = [];
indry = [];                 indcy = [];
indr2 = [];                 indc2 = [];

for i=indr %rows
    if i<(n0_old+1)
        n0_new = n0_new+1;
        indr0 = [indr0 i];
    elseif i<(n0_old+nx_old+1)
        nx_new = nx_new+1;
        indrx = [indrx i-n0_old];
    elseif i<(n0_old+nx_old+ny_old+1)
        ny_new = ny_new+1;
        indry = [indry i-nx_old-n0_old];
    elseif i<(n0_old+nx_old+ny_old+n2_old+1)
        n2_new = n2_new+1;
        indr2 = [indr2 i-ny_old-nx_old-n0_old];
    else
        error('index exceeds the number of rows in the original dopvar2d')
    end
end

for i=indc %columns
    if i<(m0_old+1)
        m0_new = m0_new+1;
        indc0 = [indc0 i];
    elseif i<(m0_old+mx_old+1)
        mx_new = mx_new+1;
        indcx = [indcx i-m0_old];
    elseif i<(m0_old+mx_old+my_old+1)
        my_new = my_new+1;
        indcy = [indcy i-mx_old-m0_old];
    elseif i<(m0_old+mx_old+my_old+m2_old+1)
        m2_new = m2_new+1;
        indc2 = [indc2 i-my_old-mx_old-m0_old];
    else
        error('index exceeds the number of columns in the original dopvar2d')
    end
end


% Construct the sliced dopvar
dopvar2d Psop;
Psop.I = Pbop.I;
Psop.var1 = Pbop.var1;
Psop.var2 = Pbop.var2;
Psop.dim = [n0_new m0_new;
            nx_new mx_new;
            ny_new my_new
            n2_new m2_new];
        
Psop.R00 = R00old(indr0,indc0);
Psop.R0x = R0xold(indr0,indcx);
Psop.R0y = R0yold(indr0,indcy);
Psop.R02 = R02old(indr0,indc2);

Psop.Rx0 = Rx0old(indrx,indc0);
Psop.Rxy = Rxyold(indrx,indcy);

Psop.Ry0 = Ry0old(indry,indc0);
Psop.Ryx = Ryxold(indry,indcx);

Psop.R20 = R20old(indr2,indc0);

for i=1:3
    Psop.Rxx{i} = Rxxold{i}(indrx,indcx);
    Psop.Rx2{i} = Rx2old{i}(indrx,indc2);
    Psop.R2x{i} = R2xold{i}(indr2,indcx);
    
    Psop.Ryy{i} = Ryyold{i}(indry,indcy);
    Psop.Ry2{i} = Ry2old{i}(indry,indc2);
    Psop.R2y{i} = R2yold{i}(indr2,indcy);
    
    for j=1:3
        Psop.R22{i,j} = R22old{i,j}(indr2,indc2);
    end
end

end



