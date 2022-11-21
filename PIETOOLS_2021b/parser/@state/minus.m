function sumTerms = minus(objA,objB)

if ~isa(objA,'state')||~isa(objB,'state')
    error('Only state type objects can be added together');
end

s.type = '.'; s.subs = 'veclength';

if sum(subsref(objA,s))~=sum(subsref(objB,s))
    error('States of unequal length cannot be added');
end

% <<<<<<< Updated upstream
isdot_A = isdot(objA); isdot_B = isdot(objB); isout_A=isout(objA); isout_B = isout(objB);
% for i=1:length(objA)
%     isdot_A = [isdot_A; objA(i).diff_order(1)*ones(subsref(objA(i),s),1)];
%     isout_A = [isout_A; strcmp(objA(i).type,'out')*ones(subsref(objA(i),s),1)];
% end
% for i=1:length(objB)
%     isdot_B = [isdot_B; objB(i).diff_order(1)*ones(subsref(objB(i),s),1)];
%     isout_B = [isout_B; strcmp(objB(i).type,'out')*ones(subsref(objB(i),s),1)];
% end
% =======
s.type = '.'; s.subs = 'veclength';
% isdot_A = []; isdot_B = []; isout_A=[]; isout_B = [];
% for i=1:length(objA)
%     isdot_A = [isdot_A; objA(i).diff_order(1)*ones(subsref(objA(i),s),1)];
%     isout_A = [isout_A; strcmp(objA(i).type,'out')*ones(subsref(objA(i),s),1)];
% end
% for i=1:length(objB)
%     isdot_B = [isdot_B; objB(i).diff_order(1)*ones(subsref(objB(i),s),1)];
%     isout_B = [isout_B; strcmp(objB(i).type,'out')*ones(subsref(objB(i),s),1)];
% end
% >>>>>>> Stashed changes
if any((isdot_A&isout_B)|(isout_A&isdot_B)|(isout_A&isout_B)|(isdot_A&isdot_B))
    error("Any linear combinations of outputs and time-derivative of state are not allowed");
end

[objC,permMatsA,permMatsB] = combine(objA,objB); % objA = permMats{1}*objC and objB = permMats{2}*objC

% objC = state2terms(objC);

opvar T; 
T.R.R0 = permMatsA-permMatsB;

sumTerms = terms(T,objC);
end