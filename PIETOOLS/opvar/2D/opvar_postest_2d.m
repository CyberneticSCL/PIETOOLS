function logval = opvar_postest_2d(P)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% logval = opvar_postest(P) tests, numerically, if P is a positive definite opvar2d 
% 
% INPUT
%   P: PI opvar2d variable
% 
% OUTPUT 
%   logval: -1 if negative definite, 0 if indefinite, 1 if positive
%   definite
% 
% NOTES:
% For support, contact M. Peet, Arizona State University at mpeet@asu.edu
% or D. Jagt at djagt@asu.edu

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PIETools - opvar_postest_2d
%
% Copyright (C)2019  M. Peet, S. Shivakumar, D. Jagt
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
% Initial coding DJ - 02_10_2021 
% ^ Based heavily on "opvar_postest" script by MMP, SS

if ~isa(P,'opvar2d')
    error('Input must be an opvar object');
end

dim = P.dim;
if any(dim(:,1)~=dim(:,2))|| ~(P==P')
    error('Opvar2d object can be sign definite if and only if it is self-adjoint')
end

m0 = dim(1,1);
mx = dim(2,1);
my = dim(3,1);
m2 = dim(4,1);
mdim = [m0;mx;my;m2];

% Define the maximal degree of the test polynomials
sz = 3;


a1 = P.I(1,1); b1 = P.I(1,2); 
a2 = P.I(2,1); b2 = P.I(2,2); 
s1 = P.var1(1);
s2 = P.var1(2);

inprod = zeros(100,1);
for i=1:100
    x_0 = 2*rand([m0,1]) - 1;       
    if mx>0
        x_x = rndm_plnml(mx,1,s1.varname,sz);
    else
        x_x = zeros(0,1);
    end
    if my>0
        x_y = rndm_plnml(my,1,s2.varname,sz);
    else
        x_y = zeros(0,1);
    end
    if m2>0
        x_2 = rndm_plnml(m2,1,[s1.varname;s2.varname],sz);
    else
        x_2 = zeros(0,1);
    end
    x = [x_0;x_x;x_y;x_2];
    
    [P1x,ndim] = opvar_apply(P,x,mdim);
    
    nval = cumsum(ndim);

    yP1x_0 = x(1:nval(1),1)'*P1x(1:nval(1),1);
    yP1x_x = int(x(nval(1)+1:nval(2),1)'*P1x(nval(1)+1:nval(2),1),s1,a1,b1);
    yP1x_y = int(x(nval(2)+1:nval(3),1)'*P1x(nval(2)+1:nval(3),1),s2,a2,b2);
    yP1x_2 = int(int(x(nval(3)+1:nval(4),1)'*P1x(nval(3)+1:nval(4),1),s1,a1,b1),s2,a2,b2);

    inprod(i) = double(yP1x_0 + yP1x_x + yP1x_y + yP1x_2);
    
end
if all(inprod>=0)
    logval=1;
elseif all(inprod<=0)
    logval=-1;
else
    logval=0;
end
end