function C = mtimes(A,B)
Aparams = A.params;
Bparams = B.params;
keyA = keys(Aparams);
keyB = keys(Bparams);

maxkeyC = 3^(length(A.var_in));

keyC = 1:maxkeyC;
paramsC = cell(1,maxkeyC);

for i=1:numel(keyA)
    for j=1:numel(keyB)
        [keysList, paramsList] = compose(Aparams{i}, Bparams{j}, A.vars, B.vars, keyA(i), keyB(j));
        for k=1:numel(keysList)
            paramsC{keyList(k)} = paramsC{keyList(k)}+paramsList{k};
        end
    end
end

% prune empty parameters
emptyIdx = cellfun('isempty', paramsC);
keyC = keyC(~emptyIdx);
paramsC = paramsC(~emptyIdx);

C = nDopvar(B.vars_in, A.vars_out, [A.dims(1), B.dims(2)], dictionary(keyC,paramsC));
end

function [keys,params] = compose(Aparam, Bparam, Avars, Bvars)
keys = [];
params = {};
end