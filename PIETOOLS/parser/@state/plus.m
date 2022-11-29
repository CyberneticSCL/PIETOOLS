function sumTerms = plus(objA,objB)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% State-class method that performs objA+objB
% Input: 
% objA, objB - state class objects
% Output:
% sumTerms - A*[objA;objB], where A=[1,1]
%
% For support, contact M. Peet, Arizona State University at mpeet@asu.edu
% or D. Jagt at djagt@asu.edu

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
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
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% If you modify this code, document all changes carefully and include date
% authorship, and a brief description of modifications
if ~isa(objA,'state')||~isa(objB,'state')
    error('Only state type objects can be added together');
end

s.type = '.'; s.subs = 'veclength';

if sum(subsref(objA,s))~=sum(subsref(objB,s))
    error('States of unequal length cannot be added');
end

isdot_A = isdot(objA); isdot_B = isdot(objB); isout_A=isout(objA); isout_B = isout(objB);

if any((isdot_A&isout_B)|(isout_A&isdot_B)|(isout_A&isout_B)|(isdot_A&isdot_B))
    error("Any linear combinations of outputs and time-derivative of state are not allowed");
end

[objC,permMatsA,permMatsB] = combine(objA,objB); % objA = permMats{1}*objC and objB = permMats{2}*objC

opvar T; T.R.R0 = permMatsA+permMatsB;
sumTerms = terms(T,objC);
end