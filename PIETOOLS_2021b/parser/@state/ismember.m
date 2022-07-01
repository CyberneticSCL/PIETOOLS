function locA = ismember(objA, objB)
locA = zeros(length(objA),1); 
for i=1:length(objA)
    for j=1:length(objB)
        if (objA(i) == objB(j))
            locA(i,1) = 1;
            break;
        end
    end
end


end