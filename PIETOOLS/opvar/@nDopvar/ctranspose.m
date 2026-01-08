function At = ctranspose(A)
Aparams = A.params;
keyA = keys(Aparams);
keyAt = zeros(1,length(keyA));
paramsAt = cell(1, length(keyAt));
for i=1:length(keyA)
    keyAt(i) = transposeKey(keyA(i),A.vars);
    paramsAt{i} = transposeParams(keyA(i),A.vars,Aparams{keyA(i)});
end

At = nDopvar(A.vars_out, A.vars_in, A.dims', dictionary(keyAt,paramsAt));
end
function k = transposeKey(i,vars)
n = length(vars.vars_in);
hi = cell(1,n);

if n==1
    hi{1} = i;
else
    [hi{:}] = ind2sub(repmat(3,1,n),i);
end

hi_flipped = hi;

hi_flipped(cellfun(@(x) x==2, hi)) = {3};
hi_flipped(cellfun(@(x) x==3, hi)) = {2};

if n==1
    k = hi_flipped{1};
else
    k = sub2ind(repmat(3,1,n),hi_flipped{:});
end
end

function Zt = transposeParams(i,vars,Z)
Zt=Z';

n = length(vars.vars_in);
hi = cell(1,n);

if n==1
    hi{1} = i;
else
    [hi{:}] = ind2sub(repmat(3,1,n),i);
end
    
for i=1:length(hi)
    if (hi{i}==2) || (hi{i}==3)
        var = vars.vars_in(i);
        dumvar = var; dumvar.varname{1} = [var.varname{1},'_dum'];
        Zt = var_swap(Zt, var, dumvar);
    end
end
end
