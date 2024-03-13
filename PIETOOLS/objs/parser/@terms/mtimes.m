function prodTerms = mtimes(K,objA)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% terms-class method that performs K*objA
% Input: 
% objA - terms class objects
% K - polynomial or double matrix
% Output:
% prodTerms - terms class object K*objA
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
if isa(K,'terms')
    error('Left multiplication by state/terms object is not supported');
end
if isa(objA,'terms')&&isa(K,'terms')
    error('Two terms type objects cannot be multiplied');
end

isdot_A = isdot(objA.statevec); isout_A=isout(objA.statevec); 
% isdot_A = []; isout_A=[]; 
% for i=1:length(objA)
%     isdot_A = [isdot_A; objA(i).diff_order(1)*ones(subsref(objA(i),s),1)];
%     isout_A = [isout_A; strcmp(objA(i).type,'out')*ones(subsref(objA(i),s),1)];
% end
if any((isdot_A|isout_A))
    error("Multiplication involving vectors with outputs or time-derivative of state is not allowed");
end

if numel(K)~=1 && size(K,2)~=length(objA)
    error('Dimensions of multiplier and terms object do not match. Cannot be multiplied');
end

opvar T; 
if numel(K)==1
    T.R.R0 = K*eye(length(objA));
else
    T.R.R0 = K;
end

prodTerms = terms(T*objA.operator,objA.statevec);
end