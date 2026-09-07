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
%   ZdC - combined unique decision-parameter vector, in first-occurrence
%         order, matching the contract of 'CombineDecisionBasis'
%
% MMP, 09/07/2026: 'unique' -> 'unique(...,''stable'')'. Plain 'unique'
% sorts, so the returned ZdC did not match the order every other routine in
% the class uses, and a caller would have had its decision variables
% silently permuted.
%
% NOTE: this routine is currently unreachable and uncalled. It sits in the
% class folder, so it is a method and dispatches on its first argument,
% which every intended caller passes as a cell -- the same trap that made
% 'lrmultiply' unreachable until it was moved to 'private'. Its work is also
% now done by 'sync_basis', which merges the decision lists for all N
% operands at once and leaves the summation to the caller.

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

% 'stable' keeps first-occurrence order, which is the contract              % MMP, 09/07/2026
% 'CombineDecisionBasis' establishes: Zd = [ZdA; setdiff(ZdB,ZdA,'stable')]. % MMP, 09/07/2026
% Plain 'unique' SORTS, so a caller would silently get a different          % MMP, 09/07/2026
% decision variable ordering than every other routine in the class.         % MMP, 09/07/2026
[ZdC,~,ic] = unique(ZdAll,'stable');                                        % MMP, 09/07/2026
%[ZdC,~,ic] = unique(ZdAll);                                                % MMP, 09/07/2026 (was)

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