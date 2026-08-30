function [CC,ZdC] = plus_decparam_batch(C,Zd)
% Given N decision parameter structures
%
%   C{k} = C{k}.A + C{k}.B' * Zd{k}
%
% this function returns CC and ZdC such that
%
%   sum_k C{k} = CC.A + CC.B' * ZdC
%
% Inputs:
%   C  - cell array {C1,...,CN}, where each C{k} has fields A and B
%   Zd - cell array {Zd1,...,ZdN} containing the corresponding parameters
%
% Outputs:
%   CC  - structure with fields A and B
%   ZdC - combined unique decision-parameter vector

N = numel(C);

% Constant part
Aout = C{1}.A;
for k = 2:N
    Aout = Aout + C{k}.A;
end

% Construct global parameter vector

n = cellfun(@numel,Zd);
offset = [0; cumsum(n(:))];

ZdAll = vertcat(Zd{:});

[ZdC,~,ic] = unique(ZdAll);

nC = numel(ZdC);
ncol = size(C{1}.B,2);

% Count total number of nonzeros
nnzTot = 0;
for k = 1:N
    nnzTot = nnzTot + nnz(C{k}.B);
end

% Construct sparse output directly
%
% Store only nonzero entries. Duplicate (row,column) entries are
% automatically summed by sparse().

I = zeros(nnzTot,1);
J = zeros(nnzTot,1);
V = zeros(nnzTot,1);

pos = 0;

for k = 1:N
    [ii,jj,vv] = find(C{k}.B);

    nk = numel(vv);
    range = pos + (1:nk);

    % Local parameter index -> global parameter index
    map = ic(offset(k) + ii);

    I(range) = map;
    J(range) = jj;
    V(range) = vv;

    pos = pos + nk;
end

Bout = sparse(I,J,V,nC,ncol);

CC = struct('A',Aout,'B',Bout);

end