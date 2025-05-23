function logval = eq(P1,P2,tol)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% logval = eq(P1,P2,tol) tests equality of two opvars P1 and P2 with tolerance tol
% Date: 6/13/19
% Version: 1.0
% 
% INPUT
% P1, P2: opvar class objects
% tol: acceptable tolerance value. If max(P1-P2)<tol, then P1=P2
% 
% OUTPUT
% logval: returns 1 if the objects are equal, 0 if not equal
%
% NOTES:
% For support, contact M. Peet, Arizona State University at mpeet@asu.edu
% or S. Shivakumar at sshivak8@asu.edu

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PIETools - eq
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
% SS, MP, 2020: Initial coding
% DJ, 12/11/2024: Bug fix Pop==0, make sure independent variables match;

if nargin<3
    tol=1e-14;
end    

if ~isa(P1,'opvar')&&(P1==0) 
    opvar P1; P1.I = P2.I; P1.dim = P2.dim;  P1.var1 = P2.var1; P1.var2 = P2.var2;  % DJ, 12/11/2024:
elseif ~isa(P2,'opvar')&&(P2==0)
    opvar P2; P2.I = P1.I; P2.dim= P1.dim;  P2.var1 = P1.var1;  P2.var2 = P1.var2;  % DJ, 12/11/2024:
elseif ~isa(P1,'opvar')|| ~isa(P2,'opvar')
    error('To check equality either both values must be opvar objects, or one of them have to be zero');
end



if any(P1.dim~=P2.dim)
    disp('Opvars have different dimensions and hence cannot be equal');
    logval=0;
end

diff = P1-P2; 
diff.P = double(diff.P);
diff.Q1 = polynomial(diff.Q1); diff.Q2 = polynomial(diff.Q2); diff.R.R0 = polynomial(diff.R.R0); 
diff.R.R1 = polynomial(diff.R.R1); diff.R.R2 = polynomial(diff.R.R2); 

diff.P(abs(diff.P)<tol)=0;
if isa(diff.Q1,'polynomial')
    diff.Q1.coefficient(abs(diff.Q1.coefficient)<tol)=0;
else
    diff.Q1(abs(diff.Q1)<tol)=0;
end
if isa(diff.Q2,'polynomial')
    diff.Q2.coefficient(abs(diff.Q2.coefficient)<tol)=0;
else
    diff.Q2(abs(diff.Q2)<tol)=0;
end
if isa(diff.R.R0,'polynomial')
    diff.R.R0.coefficient(abs(diff.R.R0.coefficient)<tol)=0;
else
    diff.R.R0(abs(diff.R.R0)<tol)=0;
end
if isa(diff.R.R1,'polynomial')
    diff.R.R1.coefficient(abs(diff.R.R1.coefficient)<tol)=0;
else
    diff.R.R1(abs(diff.R.R1)<tol)=0;
end
if isa(diff.R.R2,'polynomial')
    diff.R.R2.coefficient(abs(diff.R.R2.coefficient)<tol)=0;
else
    diff.R.R2(abs(diff.R.R2)<tol)=0;
end


try
    diff.P = double(diff.P);
    diff.Q1 = double(diff.Q1); diff.Q2 = double(diff.Q2);
    diff.R.R0 = double(diff.R.R0); diff.R.R1 = double(diff.R.R1);
    diff.R.R2 = double(diff.R.R2);
    if all(all(diff.P==0))&&all(all(diff.Q1==0))&&all(all(diff.Q2==0))&&all(all(diff.R.R0==0))&&all(all(diff.R.R1==0))&&all(all(diff.R.R2==0))
        logval=1;
    else
        logval=0;
    end
catch 
    logval=0;
end
end
