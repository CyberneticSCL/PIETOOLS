function len = stateveclength(obj)
len = zeros(length(obj),1);
for i=1:length(obj)
    len(i,1) = obj(i).veclength;
end
end