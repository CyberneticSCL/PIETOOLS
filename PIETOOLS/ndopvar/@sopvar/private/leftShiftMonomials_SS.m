function [varsC,ZC,CC] = leftShiftMonomials_SS(varsA,ZA,CA,varsB,ZB,CB)
% leftShiftMonomials_SS
%
% Computes
%
%   ((I_p \otimes ZA')*CA{I,J})*((I_q \otimes ZB')*CB{I})
%       = (I_p \otimes ZC')*CC{I,J}
%
% for each cell index I,J.
%
% INPUTS
%   varsA : variables for ZA
%   ZA    : 1-by-nA cell array of exponent vectors
%   CA    : nI-by-nJ cell array, CA{I,J} has size p*NA by q
%   varsB : variables for ZB
%   ZB    : 1-by-nB cell array of exponent vectors
%   CB    : 1-by-nI or nI-by-1 cell array, CB{I} has size q*NB by r
%
% OUTPUTS
%   varsC : merged variable list
%   ZC    : merged/product monomial basis
%   CC    : nI-by-nJ cell array, CC{I,J} has size p*NC by r

varsA = varsA(:).';
varsB = varsB(:).';
ZA = ZA(:).';
ZB = ZB(:).';

NA = prod(cellfun(@numel,ZA));
NB = prod(cellfun(@numel,ZB));

if isempty(ZA)
    NA = 1;
end
if isempty(ZB)
    NB = 1;
end

[nI,nJ] = size(CA);

if numel(CB) ~= nI
    error('leftShiftMonomials_SS: CB must have one cell for each row of CA.');
end

[mA,q] = size(CA{1,1});
[mB,r] = size(CB{1});

if mod(mA,NA) ~= 0
    error('leftShiftMonomials_SS: size(CA{1},1) must be divisible by NA.');
end

if mB ~= q*NB
    error('leftShiftMonomials_SS: size(CB{1},1) must equal size(CA{1},2)*NB.');
end

p = mA/NA;

% Merge variables. This uses MATLAB sorted union, matching current convention.
varsC = union(varsA,varsB);
nC = numel(varsC);

% Build output basis.
ZC = cell(1,nC);

for t = 1:nC
    ia = find(strcmp(varsA,varsC{t}),1);
    ib = find(strcmp(varsB,varsC{t}),1);

    if isempty(ia)
        ZC{t} = ZB{ib}(:);
    elseif isempty(ib)
        ZC{t} = ZA{ia}(:);
    else
        ea = ZA{ia}(:);
        eb = ZB{ib}(:);
        ZC{t} = unique(reshape(ea + eb.',[],1));
    end
end

NC = prod(cellfun(@numel,ZC));
if isempty(ZC)
    NC = 1;
end

% Expand ZA exponents into varsC coordinates.
EAc = zeros(NA,nC);

for t = 1:numel(varsA)
    k = find(strcmp(varsC,varsA{t}),1);

    lenA = cellfun(@numel,ZA);
    left  = prod(lenA(1:t-1));
    right = prod(lenA(t+1:end));

    EAc(:,k) = kron(ones(left,1),kron(ZA{t}(:),ones(right,1)));
end

% Expand ZB exponents into varsC coordinates.
EBc = zeros(NB,nC);

for t = 1:numel(varsB)
    k = find(strcmp(varsC,varsB{t}),1);

    lenB = cellfun(@numel,ZB);
    left  = prod(lenB(1:t-1));
    right = prod(lenB(t+1:end));

    EBc(:,k) = kron(ones(left,1),kron(ZB{t}(:),ones(right,1)));
end

% Expand ZC exponents into varsC coordinates.
EC = zeros(NC,nC);

for t = 1:nC
    lenC = cellfun(@numel,ZC);
    left  = prod(lenC(1:t-1));
    right = prod(lenC(t+1:end));

    EC(:,t) = kron(ones(left,1),kron(ZC{t}(:),ones(right,1)));
end

CC = cell(nI,nJ);

for I = 1:nI

    Bcoef = CB{I};

    if ~isequal(size(Bcoef),[q*NB,r])
        error('leftShiftMonomials_SS: inconsistent CB cell dimensions.');
    end

    for J = 1:nJ

        Acoef = CA{I,J};

        if ~isequal(size(Acoef),[p*NA,q])
            error('leftShiftMonomials_SS: inconsistent CA cell dimensions.');
        end

        G = zeros(p,r,NC);

        for ia = 1:NA

            rowsA = ia:NA:(p*NA);
            Ai = Acoef(rowsA,:);

            ea = EAc(ia,:);

            for ib = 1:NB

                rowsB = ib:NB:(q*NB);
                Bj = Bcoef(rowsB,:);

                e = ea + EBc(ib,:);

                k = find(all(EC == e,2),1);

                if isempty(k)
                    error('leftShiftMonomials_SS: internal monomial basis mismatch.');
                end

                G(:,:,k) = G(:,:,k) + Ai*Bj;
            end
        end

        CC{I,J} = reshape(permute(G,[3 1 2]),[NC*p,r]);
    end
end

end