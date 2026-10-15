% Test_intop_cat.m -- Tests @intop/blkdiag.m, vertcat.m and horzcat.m.
% Oracle: the expected dimensions and populated blocks follow from the
% input sizes.
%
% DJ (CC), 10/15/2026: Initial coding;

% NOTES
% - Structural checks only (matdim of the resulting 'tensopmat', populated
%   blocks); the action of the result on a function is not checked.
% - Blocks of blkdiag/horzcat need disjoint variable names (b%d_s%d,
%   c%d_s%d), since common_vars rejects equal names on different random
%   domains.

%%%% blkdiag: no matching requirements between pieces
k = 3;
ms = randi(3,[1,k]);    ns = randi(3,[1,k]);
Ns = randi(3,[1,k]);
Aops = cell(1,k);
for i=1:k
    % Distinct variable-name prefix per block
    vars_i = arrayfun(@(v) sprintf('b%d_s%d',i,v), 1:Ns(i), 'UniformOutput', false);
    dom_i = randi([-2,2],[Ns(i),1]) +[zeros(Ns(i),1),2*rand([Ns(i),1])];
    degs_i = randi(3,[1,Ns(i)]);
    Aops{i} = rand_intop([ms(i),ns(i)],vars_i,dom_i,degs_i,3);
end
Pop = blkdiag(Aops{:});
if ~isa(Pop,'tensopmat')
    error("blkdiag: output should be of type 'tensopmat'.")
elseif ~isequal(size(Pop),[sum(ms),sum(ns)])
    error("blkdiag: matrix dimensions of output do not match expected [sum(ms),sum(ns)].")
elseif ~isequal(size(Pop.ops),[k,k])
    error("blkdiag: number of blocks in output does not match number of inputs.")
end
for i=1:k
    for j=1:k
        if i==j
            if isempty(Pop.ops{i,j})
                error("blkdiag: diagonal block (%d,%d) should not be empty.",i,i)
            end
        elseif ~isempty(Pop.ops{i,j})
            error("blkdiag: off-diagonal block (%d,%d) should be empty.",i,j)
        end
    end
end

%%%% vertcat: pieces must share the same input (variables, domain, degrees)
N = 2;
dom = randi([-2,2],[N,1]) +[zeros(N,1),2*rand([N,1])];
degs = randi(3,[1,N]);
ms = randi(3,[1,k]);
n = randi(3);
Bops = cell(1,k);
for i=1:k
    Bops{i} = rand_intop([ms(i),n],N,dom,degs,3);
end
Pop = vertcat(Bops{:});
if ~isa(Pop,'tensopmat')
    error("vertcat: output should be of type 'tensopmat'.")
elseif ~isequal(size(Pop),[sum(ms),n])
    error("vertcat: matrix dimensions of output do not match expected [sum(ms),n].")
end

%%%% horzcat: pieces must share the same output row-count (m)
m = randi(3);
ns = randi(3,[1,k]);
Cops = cell(1,k);
for i=1:k
    Ni = randi(3);
    % Distinct variable-name prefix per block
    vars_i = arrayfun(@(v) sprintf('c%d_s%d',i,v), 1:Ni, 'UniformOutput', false);
    dom_i = randi([-2,2],[Ni,1]) +[zeros(Ni,1),2*rand([Ni,1])];
    degs_i = randi(3,[1,Ni]);
    Cops{i} = rand_intop([m,ns(i)],vars_i,dom_i,degs_i,3);
end
Pop = horzcat(Cops{:});
if ~isa(Pop,'tensopmat')
    error("horzcat: output should be of type 'tensopmat'.")
elseif ~isequal(size(Pop),[m,sum(ns)])
    error("horzcat: matrix dimensions of output do not match expected [m,sum(ns)].")
end

disp('Test_intop_cat passed: blkdiag, vertcat and horzcat all produced correctly-sized tensopmat objects.')
