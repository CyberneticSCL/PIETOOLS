function [varsC,ZC,CC] = leftShiftMonomials_SD( ...
    varsA,ZA,CA,varsB,ZB,CB)

[nI,nJ] = size(CA);

NB = prod(cellfun(@numel,ZB));
if isempty(ZB)
    NB = 1;
end

if numel(CB) ~= nI
    error('leftShiftMonomials_SD:cellDimensions', ...
        'CB must contain one cell for each row of CA.');
end

if mod(size(CB{1},1),NB) ~= 0
    error('leftShiftMonomials_SD:CBDimensions', ...
        'CB dimensions are inconsistent with ZB.');
end

q = size(CB{1},1)/NB;
nCoef = numel(CA{1,1}.A);

if mod(nCoef,q) ~= 0
    error('leftShiftMonomials_SD:CADimensions', ...
        'CA dimensions are inconsistent with CB.');
end

mA = nCoef/q;
CAconstant = cell(nI,nJ);

for i = 1:nI
    for j = 1:nJ
        param = CA{i,j};

        if numel(param.A) ~= nCoef || size(param.B,2) ~= nCoef
            error('leftShiftMonomials_SD:affineDimensions', ...
                'Affine coefficient dimensions are inconsistent.');
        end

        CAconstant{i,j} = reshape(param.A,mA,q);
    end
end

[varsC,ZC,CCconstant,Tmap] = leftShiftMonomials_SS( ...
    varsA,ZA,CAconstant,varsB,ZB,CB);

CC = cell(nI,nJ);

for i = 1:nI
    T = Tmap{i};

    for j = 1:nJ
        CC{i,j} = struct( ...
            'A',CCconstant{i,j}(:), ...
            'B',CA{i,j}.B*T.');
    end
end
end