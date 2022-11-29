function sumTerms = plus(objA,objB)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% terms-class method that performs objA+objB
% Input: 
% objA, objB - terms class objects
% Output:
% sumTerms - terms class object objA+objB
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
if isa(objA,'state')
    objA = state2terms(objA);
elseif ~isa(objA,'terms')
    error('Only state/terms type objects can be added');
end
if isa(objB,'state')
    objB = state2terms(objB);
elseif ~isa(objB,'terms')
    error('Only state/terms type objects can be added');
end

if length(objA)~=length(objB)
    error('Terms of unequal length cannot be added');
end

tempoperator = [objA.operator objB.operator];
tempstatevec = vertcat(objA.statevec, objB.statevec);

% <<<<<<< Updated upstream
isdot_A = isdot(tempstatevec); isout_A=isout(tempstatevec); 
% isdot_A = []; isout_A=[]; 
% for i=1:length(tempstatevec)
%     isdot_A = [isdot_A; tempstatevec(i).diff_order(1)*ones(tempstatevec(i).veclength,1)];
%     isout_A = [isout_A; strcmp(tempstatevec(i).type,'out')*ones(tempstatevec(i).veclength,1)];
% end
% =======
% isdot_A = []; isout_A=[]; 
% for i=1:length(tempstatevec)
%     isdot_A = [isdot_A; tempstatevec(i).diff_order(1)*ones(tempstatevec(i).veclength,1)];
%     isout_A = [isout_A; strcmp(tempstatevec(i).type,'out')*ones(tempstatevec(i).veclength,1)];
% end
% isdot_A = boolean(isdot_A); isout_A = boolean(isout_A);
% >>>>>>> Stashed changes

for i=1:sum(tempoperator.dim(:,1))
   rowval = tempoperator.R.R0(i,:);
   nonzerocols = polynomial([rowval(isdot_A),rowval(isout_A)]);
   if sum(~isequal(nonzerocols,zeros(size(nonzerocols))))>1
       error('Linear combination of time derivatives or outputs is not allowed');
   end
end

sumTerms = terms(tempoperator,tempstatevec);
end