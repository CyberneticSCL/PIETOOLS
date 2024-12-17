%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% transl_1D.m     PIETOOLS 2024
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% If you modify this code, document all changes carefully and include date
% authorship, and a brief description of modifications
%
% Initial coding YP - 6_1_2021

function R = transl_1D(R,a,b,c,d,dir,flag)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% R = transl_1D(R,a,b,c,d,dir,flag) takes in a PI operator R that acts on
% functions on the interval [a,b] and changes it to a PI operator that
% acts on functions on the inteval [c,d].
%
% INPUT:
%
% R: a 3-PI operator
% [a,b]: original interal
% [c,d]: new interval, default [-1,1]
% dir - 'x' or 'y' - variables on which 3-PI operator is acting 
% flag - user-defined option: 1, 2, or 3.
% flag = 1 - multiplicative operator in one variable R(s)
% flag = 2 - integrative operator in one variable R(s)
% flag = 3 - integrative operator in two variables R(s, theta)
% 
% NOTES:
% For support, contact M. Peet, Arizona State University at mpeet@asu.edu,
% or Y. Peet at ypeet@asu.edu
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PIETools - transl_1D
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


pvar s sbar theta thetabar s1 s1_dum s2 s2_dum

if (dir=='x')
    s=s1;
    theta=s1_dum;
else
    s=s2;
    theta=s2_dum;
end

if (flag==1)
R = subs(R, s, ((b-a)/(d-c))*sbar+(a*d-b*c)/(d-c));
R = subs(R, sbar, s);
elseif (flag==2)
R = ((b-a)/(d-c))*subs(R, s, ((b-a)/(d-c))*sbar+(a*d-b*c)/(d-c));
R = subs(R, sbar, s);
else
R = ((b-a)/(d-c))*subs(subs(R, s, ((b-a)/(d-c))*sbar+(a*d-b*c)/(d-c)), theta, ((b-a)/(d-c))*thetabar+(a*d-b*c)/(d-c));
R = subs(subs(R, sbar, s),thetabar, theta);
end

end