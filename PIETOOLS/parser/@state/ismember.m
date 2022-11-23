function [locA, idx] = ismember(objA, objB)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% State-class method that specifies if objA is in objB.
% Input: 
% objA, objB - state class objects
% Output:
% locA - logical array specifying if vector, objA, has components in objB
% idx - array specifying location of objA(i) in objB if objA(i) is present in objB

locA = zeros(length(objA),1); idx = locA;
for i=1:length(objA)
    for j=1:length(objB)
        if (isequal(objA(i),objB(j)))
            locA(i,1) = 1;
            idx(i,1) = j;
            break;
        end
    end
end
end