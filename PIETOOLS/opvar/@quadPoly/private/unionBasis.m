% @quadPoly/private/unionBasis.m
function [ZU, mapA, mapB] = unionBasis(ZA, ZB)
keyA = keys(ZA);
keyB = keys(ZB);

ZU   = ZA;
mapA = (1:size(ZA,1)).';
mapB = zeros(size(ZB,1),1);

key2idx = containers.Map('KeyType','char','ValueType','double');
for i = 1:size(ZA,1)
    key2idx(keyA{i}) = i;
end

KU = size(ZU,1);
for j = 1:size(ZB,1)
    k = keyB{j};
    if isKey(key2idx, k)
        mapB(j) = key2idx(k);
    else
        KU = KU + 1;
        ZU(KU,:) = ZB(j,:);
        key2idx(k) = KU;
        mapB(j) = KU;
    end
end
end

function k = keys(Z)
n = size(Z,1);
k = cell(n,1);
for i = 1:n
    k{i} = sprintf('%d,', Z(i,:));
end
end
