function At = ctranspose(A)
Aparams = A.params;
keyA = keys(Aparams);

paramsAt = cell(1, numel(keyA));

valA = Aparams(keyA);
varsAll = vars(A);

% in ctranspose, before the loop:
varObj    = cell(size(varsAll));
dumVarObj = cell(size(varsAll));
for j = 1:numel(varsAll)
    varObj{j}    = pvar(varsAll{j});
    dumVarObj{j} = pvar([varsAll{j}, '_dum']);
end


% transpose keys
[keyAt, keyAidx] = transposeKeys(keyA, varsAll, A.vars_in, A.vars_out);
% transpose parameters
for i=1:length(keyA)                        % this can be parallelized; parfor loops?
    paramsAt{i} = transposeParams(keyAidx(i,:), varObj, dumVarObj, valA{i});
end

At = nDopvar(A.vars_out, A.vars_in, A.dims', dictionary(keyAt,paramsAt));
end

function [keyAt, i_multi, i_multiFlip] = transposeKeys(keyA, vars, vars_in, vars_out)

n = numel(vars);

i_multi = nDopvar.keys2index(keyA,n);  % keys to multi-index

i_multiFlip = i_multi; 

% lower and upper integrals flip
i_multiFlip(i_multi==1) = 2;
i_multiFlip(i_multi==2) = 1;

vars_in_idx = ismember(vars,vars_in);
vars_out_idx = ismember(vars,vars_out);

i_multiFlip(i_multi(:,vars_in_idx&~vars_out_idx)==3) = 0;
i_multiFlip(i_multi(:,vars_out_idx&~vars_in_idx)==0) = 3;

keyAt = nDopvar.index2keys(i_multiFlip,n);
end


function Zt = transposeParams(keyidx,varObj, dumvarObj,Z)
Zt=Z';

swapPos = find(keyidx == 1 | keyidx == 2);

for p=1:numel(swapPos)
    i = swapPos(p);
    Zt = var_swap(Zt, varObj{i}, dumvarObj{i});
end
end