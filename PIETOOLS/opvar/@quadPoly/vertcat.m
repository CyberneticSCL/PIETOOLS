% @quadPoly/vertcat.m
function H = vertcat(varargin)

A = varargin;
if numel(A) == 1
    H = A{1};
    return;
end

n = A{1}.dim(2);
ns = A{1}.ns;
nt = A{1}.nt;

ZsU = A{1}.Zs;
ZtU = A{1}.Zt;

mapS = cell(numel(A),1);
mapT = cell(numel(A),1);
mapS{1} = (1:size(A{1}.Zs,1)).';
mapT{1} = (1:size(A{1}.Zt,1)).';

mTotal = A{1}.dim(1);

for i = 2:numel(A)
    if A{i}.dim(2) ~= n
        error('vertcat: column dimensions must match.');
    end
    mTotal = mTotal + A{i}.dim(1);

    [ZsU, ~, mapS{i}] = unionBasis(ZsU, A{i}.Zs);
    [ZtU, ~, mapT{i}] = unionBasis(ZtU, A{i}.Zt);
end

dsU = size(ZsU,1);
dtU = size(ZtU,1);

% Fast path if bases already identical
sameBases = true;
for i = 1:numel(A)
    if ~isequal(A{i}.Zs, ZsU) || ~isequal(A{i}.Zt, ZtU)
        sameBases = false;
        break;
    end
end

Ccells = cell(numel(A),1);
for i = 1:numel(A)
    mi = A{i}.dim(1);

    if sameBases
        Ccells{i} = A{i}.C;
    else
        ds = size(A{i}.Zs,1);
        dt = size(A{i}.Zt,1);
        Ccells{i} = lift(A{i}.C, mi, n, ds, dt, dsU, dtU, mapS{i}, mapT{i});
    end
end

CH = cat(1, Ccells{:});
H  = quadPoly(CH, ZsU, ZtU, [mTotal n], ns, nt);

end
