function logval = isvalidPolynomial(exp)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% isvalidPolynomial.m     PIETOOLS 2022
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This function tests whether a mathematical expression exp is a valid
% polynomial in s and theta
% if logval is 0, then exp is a valid polynomial
% logval is 2, exp has unidentified polynomial variable different from s or theta
% logval is 3, exp has polynomial variable with non-integer power
% logval is 4, exp has other error such as mismatched brackets, incorrect
% operands etc.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Copyright (C)2022  M. Peet, S. Shivakumar, D. Jagt
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
% DEVELOPER LOGS:
% If you modify this code, document all changes carefully and include date
% authorship, and a brief description of modifications
pvar s theta;

exp = convertCharsToStrings(exp); % convert to string

logval=0;

try
    tmp=eval(exp);
catch err
    if strcmp(err.message,"B must be a constant.")
        logval=3;    
    elseif contains(err.message,"Unrecognized function or variable")
        logval=2;
    elseif strcmp(err.message,"Invalid use of power for polynomials")
        logval=3;
    else
        logval=4;
    end
end
end