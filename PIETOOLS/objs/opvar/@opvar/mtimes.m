function [Pcomp] = mtimes(P1,P2) 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% [Pcomp] = mtimes(P1,P2)  performs takes in two operators, 
% P1:R^p x L2^q to R^m x L2^n and P2:R^m x L2^n to R^i x L2^j. 
% It returns the composition P1 o P2: R^p x L2^q to R^i x L2^j
% Version 1.0
% Date: 6/13/19
% 
% INPUT
% P1, P2: opvar class objects with same inner dimension
%
% OUTPUT
% Pcomp: the matlab structure which is the equivalent operator of
% P1 o P2, where 'o' represents composition of operators.
%
% NOTES:
% For support, contact M. Peet, Arizona State University at mpeet@asu.edu
% or S. Shivakumar at sshivak8@asu.edu

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PIETools - mtimes
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
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% If you modify this code, document all changes carefully and include date
% authorship, and a brief description of modifications
%
% Initial coding MMP, SS  - 7_26_2019
% Changed dummy polynomials names from s, theta to ds, dtheta - SS 6/29/20
% Adjusted so that dpvar times opvar returns dopvar.



if isa(P1,'opvar') && isa(P2,'opvar')
    P1.dim = P1.dim; P2.dim = P2.dim;
if any(P1.dim(:,2)~=P2.dim(:,1))
    error("Composition requires inner dimensions of the operators to match");
    return
end
if any(P1.I~=P2.I)
    error('Operators act on different intervals and cannot be composed');
end
a=P1.I(1);
b=P1.I(2);
ds = P1.var1;
dtheta = P1.var2;
I = P1.I;

opvar Pcomp; Pcomp.I = I; Pcomp.var1 = ds; Pcomp.var2 = dtheta;

Pcomp.P = P1.P*P2.P + int(P1.Q1*P2.Q2,ds,a,b);
Pcomp.Q1 = P1.P*P2.Q1 + P1.Q1*P2.R.R0+ int(var_swap(P1.Q1*P2.R.R1,ds,dtheta),dtheta,ds,b)+...
            int(var_swap(P1.Q1*P2.R.R2,ds,dtheta),dtheta,a,ds);
Pcomp.Q2 = P1.Q2*P2.P + P1.R.R0*P2.Q2 + int(P1.R.R1*subs(P2.Q2,ds,dtheta),dtheta,a,ds)+...
                    int(P1.R.R2*subs(P2.Q2,ds,dtheta),dtheta,ds,b);
Ptemp = PL2L_compose(P1,P2);
Pcomp.R.R0 = Ptemp.R.R0;
Pcomp.R.R1 = P1.Q2*subs(P2.Q1,ds,dtheta)+Ptemp.R.R1;
Pcomp.R.R2 = P1.Q2*subs(P2.Q1,ds,dtheta)+Ptemp.R.R2;

elseif ~isa(P2,'opvar') %multiplication of operator times matrix
    ds = P2.var1;
    dtheta = P2.var2;
    I = P2.I;
    opvar Pcomp; Pcomp.I = I; Pcomp.var1 = ds; Pcomp.var2 = dtheta;
    if all(size(P2)==[1,1]) %scalar multiplication
        Pcomp.P = P2*P1.P;
        Pcomp.Q1 = P2*P1.Q1;
        Pcomp.Q2 = P2*P1.Q2;
        Pcomp.R.R0 = P2*P1.R.R0;
        Pcomp.R.R1 = P2*P1.R.R1;
        Pcomp.R.R2 = P2*P1.R.R2;
    else
        if size(P2,1)~=sum(P1.dim(:,2))
            error("Multiplication requires inner dimensions of the operators to match");
            return
        end

        r = P1.dim(1,2); p = P1.dim(2,2);
        idxr = 1:r; idxp = r+1:r+p;
        P2r = P2(idxr,:); P2p = P2(idxp,:);
        
        Pcomp.P = P1.P*P2r;
        Pcomp.Q2 = P1.Q2*P2r;
        
        Pcomp.Q1 = P1.Q1*P2p;
        Pcomp.R.R0 = P1.R.R0*P2p;
        Pcomp.R.R1 = P1.R.R1*P2p;
        Pcomp.R.R2 = P1.R.R2*P2p;
    end
    if isa(P2,'dpvar')  % DJ, 12/30/2021
        Pcomp = opvar2dopvar(Pcomp);
    end
else %multiplication of matrix times the operator
    ds = P2.var1;
    dtheta = P2.var2;
    I = P2.I;
    opvar Pcomp; Pcomp.I = I; Pcomp.var1 = ds; Pcomp.var2 = dtheta;
    if all(size(P1)==[1,1]) %scalar times operator
        Pcomp.P = P1*P2.P;
        Pcomp.Q1 = P1*P2.Q1;
        Pcomp.Q2 = P1*P2.Q2;
        Pcomp.R.R0 = P1*P2.R.R0;
        Pcomp.R.R1 = P1*P2.R.R1;
        Pcomp.R.R2 = P1*P2.R.R2;
    else
        if size(P1,2)~=sum(P2.dim(:,1))
            error('Multiplication requires inner dimensions of the operators to match');
            return
        end
        r = P2.dim(1,1); p = P2.dim(2,1);
        idxr = 1:r; idxp = r+1:r+p;
        P1r = P1(:,idxr); P1p = P1(:,idxp);
        if isempty(P1r)
            P1r = zeros(0,r);
        end
        if isempty(P1p)
            P1p = zeros(0,p);
        end
        Pcomp.P = P1r*P2.P;
        Pcomp.Q1 = P1r*P2.Q1;
        Pcomp.Q2 = P1p*P2.Q2;
        Pcomp.R.R0 = P1p*P2.R.R0;
        Pcomp.R.R1 = P1p*P2.R.R1;
        Pcomp.R.R2 = P1p*P2.R.R2;
    end
    if isa(P1,'dpvar')  % DJ, 12/30/2021
        Pcomp = opvar2dopvar(Pcomp);
    end
end 
end
function [Ph] = PL2L_compose(P1,P2)                                 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This function performs the composition of two operators on L2 and returns
% an equivalent operator on L2
% Version 1.0
% Date: 6/2/19
% 
% Inputs:
% P1: Matrices and kernel of second operator
% P2: Matrices and kernel of first operator
% 
% Output:
% Returns the composition operator structure.

a=P1.I(1);
b=P1.I(2);
ds = P1.var1;
dtheta = P2.var2;
A = P1.R.R0; B1 = P1.R.R1; B2 = P1.R.R2;
M = P2.R.R0; N1 = P2.R.R1; N2 = P2.R.R2;
Mh = A*M;
pvar dbeta;
B1sb = subs(B1,dtheta,dbeta); B2sb = subs(B2,dtheta,dbeta);
N1bt = subs(N1,ds,dbeta); N2bt = subs(N2,ds,dbeta);
N1h = A*N1+B1*subs(M,ds,dtheta)+...
                               int(B1sb*N1bt,dbeta,dtheta,ds)+...
                               int(B1sb*N2bt,dbeta,a,dtheta)+...
                               int(B2sb*N1bt,dbeta,ds,b);
N2h = A*N2+B2*subs(M,ds,dtheta)+...
                                int(B2sb*N1bt,dbeta,dtheta,b)+...
                                int(B2sb*N2bt,dbeta,ds,dtheta)+...
                                int(B1sb*N2bt,dbeta,a,ds);
Ph.R.R0 = Mh; Ph.R.R1 = N1h; Ph.R.R2 = N2h;
end