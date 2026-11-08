function dP = diff_opvar(P)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% DP = DIFF_OPVAR(P) computes the composition DP of a differential operator
% d/ds with a PI operator P, so that d/ds (P*x)(s) = (dP*[x;dx/ds])(s).
% 
% INPUT 
%   P: opvar object
% 
% OUTPUT 
%   dP: derivative of P with respect to variable P.var1
% 
% NOTES:
% For support, contact M. Peet, Arizona State University at mpeet@asu.edu

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PIETOOLS - diff_opvar
%
% Copyright (C)2024 PIETOOLS Team
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
% Initial coding MMP, SS  - 7_7_2020
% DJ, 03/24/2025: Make sure parameters are indeed polynomial.

if ~isa(P,'opvar')
    error("Input diff function must be an opvar object");
end

opvar dP;
dP.I = P.I;
dP.dim = [P.dim(1,1) P.dim(1,2); P.dim(2,1) 2*P.dim(2,2)];
ss = P.var1;    st = P.var2;
dP.var1 = ss;   dP.var2 = st;                                               % DJ, 03/24/2025

dP.P = zeros(size(P.P)); 
dP.Q1 = polynomial(zeros(dP.dim(1,1),dP.dim(2,2)));
if isa(P.Q2,'polynomial') && ~isempty(P.Q2)                                 % DJ, 03/24/2025
    dP.Q2 = diff(P.Q2,ss); 
else
    dP.Q2 = polynomial(zeros(size(P.Q2)));
end
if isa(P.R.R0,'polynomial') && ~isempty(P.R.R0)                             % DJ, 03/24/2025
    dP.R.R0 = [diff(P.R.R0,ss)+subs(P.R.R1-P.R.R2,st,ss) P.R.R0];
else
    dP.R.R0 = [subs(polynomial(P.R.R1-P.R.R2),st,ss),polynomial(zeros(size(P.R.R0)))];
end
if isa(P.R.R1,'polynomial') && ~isempty(P.R.R1)                             % DJ, 03/24/2025
    dP.R.R1 = [diff(P.R.R1,ss) zeros(size(P.R.R1))];
else
    dP.R.R1 = polynomial(zeros(size(P.R.R1,1),2*size(P.R.R1,2)));
end
if isa(P.R.R2,'polynomial') && ~isempty(P.R.R2)                             % DJ, 03/24/2025
    dP.R.R2 = [diff(P.R.R2,ss) zeros(size(P.R.R2))];
else
    dP.R.R2 = polynomial(zeros(size(P.R.R2,1),2*size(P.R.R2,2)));
end