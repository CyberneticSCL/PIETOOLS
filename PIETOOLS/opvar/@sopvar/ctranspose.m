function At = ctranspose(A)
At = A;
At.vars_out = A.vars_in;
At.vars_in = A.vars_out;
At.params = cell(size(A.params,2), size(A.params,1));

varsMain = union(A.vars_in,A.vars_out);
varsDummy = cellfun(@(x) [x,'_dum'], varsMain,UniformOutput=false);

nvars = size(A.params);

for i=1:numel(nvars)
    for j=1:nvars(j)
        At.params{j} = A.params{j}';
        if j~=1  % if integral terms then swap main and dummy variables
            At.params{j} = var_swap(At.params{j}, varsMain{i}, varsDummy{i});
        end        
    end
end
end
