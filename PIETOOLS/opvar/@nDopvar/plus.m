function C = plus(A,B)
Aparams = A.params;
Bparams = B.params;
keyA = keys(Aparams)';
keyB = keys(Bparams)';

keyC = union(keyA,keyB);
paramsC = cell(1,length(keyC));

for i=1:length(keyC)
    if ~ismember(keyC(i),keyA) 
        paramsC{i} = Bparams{keyC(i)};
    elseif ~ismember(keyC(i),keyB) 
        paramsC{i} = Aparams{keyC(i)};
    else
        paramsC{i} = Aparams{keyC(i)}+Bparams{keyC(i)}; 
    end
end

C = nDopvar(A.vars_in, A.vars_out, A.dims, dictionary(keyC, paramsC));
end