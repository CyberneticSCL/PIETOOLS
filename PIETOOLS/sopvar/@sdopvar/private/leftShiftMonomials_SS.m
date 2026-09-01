function [varsC,ZC,CC,Tmap] = leftShiftMonomials_SS(varsA,ZA,CA,varsB,ZB,CB)
% leftShiftMonomials_SS
%
% Computes
%
%   ((I_p \otimes ZA')*CA{I,J})*((I_q \otimes ZB')*CB{I})
%       = (I_p \otimes ZC')*CC{I,J}
%
% for each cell index I,J.
%
% Optional 4th output:
%
%   vec(CC{I,J}) = Tmap{I} * vec(CA{I,J})
%
% This is useful when CA{I,J} is not numeric but represented through
% decision-parameter coefficients.

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

varsC = union(varsA,varsB);
nC = numel(varsC);

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

EAc = zeros(NA,nC);

for t = 1:numel(varsA)
    k = find(strcmp(varsC,varsA{t}),1);

    lenA = cellfun(@numel,ZA);
    left  = prod(lenA(1:t-1));
    right = prod(lenA(t+1:end));

    EAc(:,k) = kron(ones(left,1),kron(ZA{t}(:),ones(right,1)));
end

EBc = zeros(NB,nC);

for t = 1:numel(varsB)
    k = find(strcmp(varsC,varsB{t}),1);

    lenB = cellfun(@numel,ZB);
    left  = prod(lenB(1:t-1));
    right = prod(lenB(t+1:end));

    EBc(:,k) = kron(ones(left,1),kron(ZB{t}(:),ones(right,1)));
end

EC = zeros(NC,nC);

for t = 1:nC
    lenC = cellfun(@numel,ZC);
    left  = prod(lenC(1:t-1));
    right = prod(lenC(t+1:end));

    EC(:,t) = kron(ones(left,1),kron(ZC{t}(:),ones(right,1)));
end

if nC == 0
    Qmat = sparse(1,1,1,1,1);
else
    EC3dt = reshape(EAc,1,[],size(EAc,2)) + reshape(EBc,[],1,size(EBc,2));
    EC2dt = reshape(EC3dt,[],size(EC3dt,3));

    [~,Idx_array] = ismember(EC2dt,EC,'rows');

    if any(Idx_array == 0)
        error('leftShiftMonomials_SS: failed to map product exponents into ZC.');
    end

    Qmat = sparse(Idx_array,1:(NA*NB),ones(size(Idx_array)),NC,NA*NB);
end

IpQ = kron(speye(p),Qmat);

CC = cell(nI,nJ);

if nargout >= 4
    Tmap = cell(nI,1);
else
    Tmap = [];
end

for I = 1:nI

    Bcoef = CB{I};

    if ~isequal(size(Bcoef),[q*NB,r])
        error('leftShiftMonomials_SS: inconsistent CB cell dimensions.');
    end

    if nargout >= 4
        T = sparse(p*NC*r,p*NA*q);

        for k = 1:NB
            Dk = IpQ(:,k:NB:end);
            Ek = Bcoef(k:NB:end,:);
            T = T + kron(Ek.',Dk);
        end

        Tmap{I} = T;
    end

    for J = 1:nJ

        Acoef = CA{I,J};

        if ~isequal(size(Acoef),[p*NA,q])
            error('leftShiftMonomials_SS: inconsistent CA cell dimensions.');
        end

        if nargout >= 4
            CC{I,J} = reshape(Tmap{I}*Acoef(:),p*NC,r);
        else
            CC{I,J} = IpQ*kron(Acoef,speye(NB))*Bcoef;
        end
    end
end

end