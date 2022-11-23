function logval = opvar_postest(P)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% logval = opvar_postest(P) tests, numerically, if P is a positive definite opvar 
% 
% INPUT
%   P: PI opvar variable
% 
% OUTPUT 
%   logval: -1 if negative definite, 0 if indefinite, 1 if positive
%   definite
% 
% NOTES:
% For support, contact M. Peet, Arizona State University at mpeet@asu.edu
% or S. Shivakumar at sshivak8@asu.edu

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PIETools - opvar_postest
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

if ~isa(P,'opvar')
    error('Input must be an opvar object');
end

dim = P.dim;
if any(dim(:,1)~=dim(:,2))|| ~(P==P')
    error('Opvar object can be sign definite if and only if it is self-adjoint')
end

dim = dim(:,1);
a = P.I(1); b=P.I(2); s = P.var1; theta = P.var2;
for i=1:100
    x = -10+20*rand(dim(1),1);
    xmbf = (-10+20*rand(dim(2),10))*monomials(s,[0:9]);
    inprod(i) = double(x'*P.P*x + int(x'*P.Q1*xmbf,s,a,b) + int(xmbf'*P.Q2*x,s,a,b)...
        + int(xmbf'*P.R.R0*xmbf,s,a,b)+ int(int(xmbf'*P.R.R1*subs(xmbf,s,theta),theta,a,s),s,a,b)...
        + int(int(xmbf'*P.R.R1*subs(xmbf,s,theta),theta,s,b),s,a,b));
end
if all(inprod>=0)
    logval=1;
elseif all(inprod<=0)
    logval=-1;
else
    logval=0;
end
end