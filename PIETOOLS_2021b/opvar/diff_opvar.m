function [dP] = diff_opvar(P)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% [[dP] = diff_opvar(P) calculates spatial derivative of (Px)(s) with respect to
% s
% 
% INPUT 
%   P: opvar object
% 
% OUTPUT 
%   dP: derivative of P with respect to pvar s
% 
% NOTES:
% For support, contact M. Peet, Arizona State University at mpeet@asu.edu

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PIETools - diff_opvar
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
% Initial coding MMP, SS  - 7_7_2020

if ~isa(P,'opvar')
    error("Input diff function must be an opvar object");
end

opvar dP;
dP.I = P.I;
dP.dim = [P.dim(1,1) P.dim(1,2); P.dim(2,1) 2*P.dim(2,2)];
ss = P.var1; st = P.var2;
dP.P = zeros(size(P.P)); 
dP.Q1 = zeros(dP.dim(1,1),dP.dim(2,2));
dP.Q2 = diff(P.Q2,ss); 
dP.R.R0 = [diff(P.R.R0,ss)+subs(P.R.R1-P.R.R2,st,ss) P.R.R0];
dP.R.R1 = [diff(P.R.R1,ss) zeros(size(P.R.R1))];
dP.R.R2 = [diff(P.R.R2,ss) zeros(size(P.R.R2))];