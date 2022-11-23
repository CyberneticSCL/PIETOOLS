function sumTerms = plus(objA,objB)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% terms-class method that performs objA+objB
% Input: 
% objA, objB - terms class objects
% Output:
% sumTerms - terms class object objA+objB

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