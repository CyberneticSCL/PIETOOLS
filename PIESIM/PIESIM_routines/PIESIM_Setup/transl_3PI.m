function R = transl_3PI(R,a,b,c,d,dir)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% R = transl_3PI(T,R,a,b,c,d,dir) takes in a 3-PI operator R that acts on
% functions in the interval [a,b] and changes it to a 3-PI operator that
% acts on functions on the interval [c,d].
%
% INPUT:
%
% R: a 3-PI operator
% [a,b]: original interal
% [c,d]: new interval, default [-1,1]
% dir - 'x' or 'y' - variables on which 3-PI operator is acting 
% 
% NOTES:
% For support, contact M. Peet, Arizona State University at mpeet@asu.edu,
% or Y. Peet at ypeet@asu.edu

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PIETools - transl_3PI
%
% Copyright (C)2024  M. Peet, Y. Peet
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


pvar s theta sbar thetabar s1 s1_dum s2 s2_dum theta1 theta2;

if (dir=='x')
    s=s1;
    theta=s1_dum;
    var=R{2}.varname;
    if length(var)==2
    if (var{2}(1)=='t')
        theta=theta1;
    end
    end
else
    s=s2;
    theta=s2_dum;
    var=R{2}.varname;
    if length(var)==2
    if (var{2}(1)=='t')
        theta=theta2;
    end
    end
end


R{1} = subs(R{1}, s, ((b-a)/(d-c))*sbar+(a*d-b*c)/(d-c));
R{1} = subs(R{1}, sbar, s);
    
R{2} = ((b-a)/(d-c))*subs(subs(R{2}, s, ((b-a)/(d-c))*sbar+(a*d-b*c)/(d-c)), theta, ((b-a)/(d-c))*thetabar+(a*d-b*c)/(d-c));
R{2} = subs(subs(R{2}, sbar, s),thetabar, theta);

R{3} = ((b-a)/(d-c))*subs(subs(R{3}, s, ((b-a)/(d-c))*sbar+(a*d-b*c)/(d-c)), theta, ((b-a)/(d-c))*thetabar+(a*d-b*c)/(d-c));
R{3} = subs(subs(R{3}, sbar, s),thetabar,theta);

