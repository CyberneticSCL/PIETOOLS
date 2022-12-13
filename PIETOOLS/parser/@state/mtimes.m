function prodTerms = mtimes(K,objA)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% State-class method that performs K*objA
% Input: 
% objA - state class objects
% K - polynomial or double
% Output:
% prodTerms - A*[objA], where A=K
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

if isa(objA,'state')&&isa(K,'state')
    error('Two state type objects cannot be multiplied');
end

if isa(K,'state')
    error('Left multiplication by a state object is not supported');
end

isdot_A = isdot(objA); isout_A=isout(objA); 
% isdot_A = []; isout_A=[];
% for i=1:length(objA)
%     isdot_A = [isdot_A; objA(i).diff_order(1)*ones(subsref(objA(i),s),1)];
%     isout_A = [isout_A; strcmp(objA(i).type,'out')*ones(subsref(objA(i),s),1)];
% end
if any((isdot_A|isout_A))
    error("Multiplication involving vectors with outputs or time-derivative of state is not allowed");
end
s.type = '.'; s.subs = 'veclength';

if numel(K)~=1 && (size(K,2)~=sum(subsref(objA,s)))&& (sum(subsref(objA,s))~=1)
    error('Dimensions of multiplier and state vector do not match. Cannot be multiplied');
end

opvar T;
if numel(K)==1
    s.type = '.'; s.subs = 'veclength';
    T.R.R0 = K*eye(sum(subsref(objA,s)));
else
    T.R.R0 = K;
end
prodTerms = terms(T,objA);
end