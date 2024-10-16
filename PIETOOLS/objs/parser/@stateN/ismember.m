function [locA, idx] = ismember(objA, objB)
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