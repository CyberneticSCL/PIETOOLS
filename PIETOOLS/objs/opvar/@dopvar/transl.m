function T = transl(T,I)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% [Pcat] = transl(T,[c,d]) takes in a 4-PI operator T that acts on
% functions in the interval [a,b] and changes it to a 4-PI operator that
% acts of functions on the interval [c,d].
%
% INPUT:
%
% T: a dopvar object
% [c,d]: interval, default [-1,1]
% 
% NOTES:
% For support, contact M. Peet, Arizona State University at mpeet@asu.edu
% or S. Shivakumar at sshivak8@asu.edu

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PIETools - transl
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



if nargin==1
    I = [-1,1];
elseif nargin>2
    error("Incorrect number of inputs; only 2 inputs are allowed");
end
if ~isa(T,'dopvar')
    error("First argument must be a opvar class object");
end
if length(I)~=2
    error("Second argument must be an array of length 2");
end

I_init = T.I;
a = I_init(1);
b = I_init(2);

c = I(1);
d = I(2);

pvar s theta sbar thetabar;

T.Q1 = ((b-a)/(d-c))*subs(T.Q1, s, ((b-a)/(d-c))*sbar+(a*d-b*c)/(d-c));
T.Q1 = subs(T.Q1, sbar, s);

T.Q2 = subs(T.Q2, s, ((b-a)/(d-c))*sbar+(a*d-b*c)/(d-c));
T.Q2 = subs(T.Q2, sbar, s);

T.R.R0 = subs(T.R.R0, s, ((b-a)/(d-c))*sbar+(a*d-b*c)/(d-c));
T.R.R0 = subs(T.R.R0, sbar, s);
    
T.R.R1 = subs(subs(T.R.R1, s, ((b-a)/(d-c))*sbar+(a*d-b*c)/(d-c)), theta, ((b-a)/(d-c))*thetabar+(a*d-b*c)/(d-c));
T.R.R1 = subs(subs(T.R.R1, sbar, s),thetabar, theta);

T.R.R2 = subs(subs(T.R.R2, s, ((b-a)/(d-c))*sbar+(a*d-b*c)/(d-c)), theta, ((b-a)/(d-c))*thetabar+(a*d-b*c)/(d-c));
T.R.R2 = subs(subs(T.R.R2, sbar, s),thetabar,theta);

T.I = I;
end