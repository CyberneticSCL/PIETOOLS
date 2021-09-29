function [logval, msg] = isvalid(P)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% [logval, msg] = isvalid(P) tests if operator P: R^p x L2^q to R^m x L2^n has valid
% properties. The functions tests dimension compatibility and object
% compatibility of each component.
% 
% INPUT
% P: dopvar class object
% 
% OUTPUT
% logval: returns 0 if the object is a valid opvar
%                 1 if the object has incorrect dimensions
%                 2 if component P is not a matrix
%                 3 if R0, Q1 or Q2 are not polynomials in s
%                 4 if R1, R2 are not polynomials in s, theta
% msg: type of error
% 
% NOTES:
% For support, contact M. Peet, Arizona State University at mpeet@asu.edu
% or S. Shivakumar at sshivak8@asu.edu

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% PIETools - isvalid
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

if ~isa(P,'dopvar')
    error('To check validity input must be opvar object');
end

dim = P.dim;
Q1 = P.Q1; Q2 = P.Q2; R0 = P.R.R0; R1=P.R.R1; R2 = P.R.R2; 
opvarlist = {P.var1.varname{1}, P.var2.varname{1}}; 

if isa(Q1,'dpvar')
    Q1list = Q1.varname;
    Q1diff = setdiff(Q1list,P.var1.varname{1});
else
    Q1diff = {};
end
if isa(Q2,'dpvar')
    Q2list = Q2.varname;
    Q2diff = setdiff(Q2list,P.var1.varname{1});
else
    Q2diff = {};
end
if isa(R0,'dpvar')
    R0list = R0.varname;
    R0diff = setdiff(R0list,P.var1.varname{1});
else
    R0diff = {};
end
if isa(R1,'dpvar')
    R1list = R1.varname;
    R1diff = setdiff(R1list,opvarlist);
else
    R1diff = {};
end
if isa(R2,'polynomial')
    R2list = R2.varname;
    R2diff = setdiff(R2list,opvarlist);
else
    R2diff = {};
end

if any(isnan(dim(:)))
    logval=1;
    msg = 'Components have incompatible dimensions';
elseif ~isa(P.P,'double')
    logval=2;
    msg = 'P is not a matrix';
elseif ~isa(Q1,'double')&&~isempty(Q1diff)
    logval=3;
    msg= 'Q1 has pvars different from var1';
elseif ~isa(Q2,'double')&&~isempty(Q2diff)
    logval=3;
    msg= 'Q2 has pvars different from var1';
elseif ~isa(R0,'double')&&~isempty(R0diff)
    logval=3;
    msg= 'R0 has pvars different from var1';
elseif ~isa(R1,'double')&&~isempty(R1diff)
    logval=4;
    msg= 'R1 has pvars different from var1 and var2';
elseif ~isa(R2,'double')&&~isempty(R2diff)
    logval=4;
    msg= 'R2 has pvars different from var1 and var2';
else
    logval=0;
    msg = 'Valid Opvar';
end
end
