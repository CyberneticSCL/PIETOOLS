function [P] = subs_op(P1, old, new) 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% [P] = subs_op(P1, old, new) function performs addition of two operators P: R^p x L2^q to R^m x L2^n
% Date: 7/1/20
% Version: 1.0
% 
% INPUT
% P1: opvar class objects
% old: pvar object to be replaced
% new: value to be substituted
% 
% OUTPUT
% P: returns P after subtitution
%
% NOTES:
% For support, contact M. Peet, Arizona State University at mpeet@asu.edu
% or S. Shivakumar at sshivak8@asu.edu

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PIETools - subs
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
% Initial coding MMP, SS  - 7_1_2020
%


if ~isa(P1,'opvar')
    error('Input must be an opvar variable');
end

opvar P; P.dim = P1.dim; P.var1 = P1.var1; P.var2 = P1.var2;
P.I = P1.I;

fields = {'P', 'Q1', 'Q2', 'R0', 'R1', 'R2'};
for i=fields
if ~contains(i,'R')
if isa(P1.(i{:}),'polynomial')
    P.(i{:}) = subs(P1.(i{:}),old,new);
end
else
if isa(P1.R.(i{:}),'polynomial')
    P.R.(i{:}) = subs(P1.R.(i{:}),old,new);
end
end
end

end