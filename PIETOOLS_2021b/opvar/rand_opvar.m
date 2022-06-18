function L=rand_opvar(nmat,pvar_deg,var1,var2,I)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% L=rand_opvar(nmat,pvar_deg,var1,var2,I) function generates a random opvar object
% 
% INPUT
%   nmat: dimensions of the opvar object L
%   pvar_deg: max deg of polynomials in components of L
%   var1: pvar in Q1, Q2, R0
%   var2: pvar in R1, R2
%   I: Interval [a,b] of L
%   
% OUTPUT
%   L: random opvar object
% 
% NOTES:
% For support, contact M. Peet, Arizona State University at mpeet@asu.edu
% or S. Shivakumar at sshivak8@asu.edu

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PIETools - rand_opvar
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
% Initial coding MMP - 6_29_2020
 

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% RANDOM opvar GENERATOR %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
opvar L;
L.var1=var1; 
L.var2=var2;
L.I=I;
n1=nmat(1,1); n2=nmat(1,2);
n3=nmat(2,1); n4=nmat(2,2);

L.P=rand([n1,n2]);
L.Q1=rand([n1,n4]);
L.Q2=rand([n3,n2]);
L.R.R0=rand([n3,n4]);
L.R.R1=rand([n3,n4]);
L.R.R2=rand([n3,n4]);
for i=1:pvar_deg
    L.Q1=L.Q1+rand([n1,n4])*var1^i;
    L.Q2=L.Q2+rand([n3,n2])*var1^i;
    L.R.R0=L.R.R0+rand([n3,n4])*var1^i;
    for j=1:pvar_deg
        L.R.R1=L.R.R1+rand([n3,n4])*var1^i*var2*j;
        L.R.R2=L.R.R2+rand([n3,n4])*var1^i*var2*j;
    end
end