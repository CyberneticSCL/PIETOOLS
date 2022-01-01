function Psop=op_slice(Pbop,indr,indc)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Psop=op_slice(Pbop,indr,indc) transposes an operator P: R^p x L2^q to R^m x L2^n
% Date: 7/1/2020
% Version: 1.0
% 
% INPUT
% Pbop: opvar class object to slice
% indr: row index of the slice
% indc: column index of the slice
% 
% OUTPUT
% Psop: slice of the opvar object
%  
% NOTES:
% For support, contact M. Peet, Arizona State University at mpeet@asu.edu
% or S. Shivakumar at sshivak8@asu.edu

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PIETools - ctranspose
%
% Copyright (C)2019  M. Peet, S. Shivakumar
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
% Initial coding MMP, SS  - 7_26_2019
%
dim_old=Pbop.dim;
Pold=Pbop.P;
Q1old=Pbop.Q1;
Q2old=Pbop.Q2;
R0old=Pbop.R.R0;
R1old=Pbop.R.R1;
R2old=Pbop.R.R2;

nr1_old=dim_old(1,1);
nr2_old=dim_old(2,1);
nc1_old=dim_old(1,2);
nc2_old=dim_old(2,2);
nr1_new=0;
nr2_new=0;
nc1_new=0;
nc2_new=0;
indr1=[];
indr2=[];
indc1=[];
indc2=[];

for i=indr %rows
    if i<(nr1_old+1)
        nr1_new=nr1_new+1;
        indr1=[indr1 i];
    elseif i<(nr1_old+nr2_old+1)
        nr2_new=nr2_new+1;
        indr2=[indr2 i-nr1_old];
    else
        error('index exceeds the number of rows in the original opvar')
    end
end
for i=indc %columns
    if i<(nc1_old+1)
        nc1_new=nc1_new+1;
        indc1=[indc1 i];
    elseif i<(nc1_old+nc2_old+1)
        nc2_new=nc2_new+1;
        indc2=[indc2 i-nc1_old];
    else
        error('index exceeds the number of columns in the original opvar')
    end
end
% Construct the sliced opvar
opvar Psop;
Psop.I=Pbop.I;
Psop.var1=Pbop.var1;
Psop.var2=Pbop.var2;
Psop.dim=[nr1_new nc1_new;
    nr2_new nc2_new];
Psop.P=Pold(indr1,indc1);
Psop.Q1=Q1old(indr1,indc2);
Psop.Q2=Q2old(indr2,indc1);
Psop.R.R0=R0old(indr2,indc2);
Psop.R.R1=R1old(indr2,indc2);
Psop.R.R2=R2old(indr2,indc2);

