function At = ctranspose(A)
% initialize the transpose object
At = A;  
At.vars_out = A.vars_in;  % input vars become output vars
At.vars_in = A.vars_out;  % output vars become input vars

At.dims = A.dims';       % matrix dimension is transposed

% collect all vars
varsMain = union(A.vars_in,A.vars_out); 
% create dummy vars corresponding to vars
varsDummy = cellfun(@(x) [x,'_dum'], varsMain,UniformOutput=false);

% total vars
nvars = numel(varsMain);
% parameter cell dimension
cellsize = size(A.params);

% for each parameter in cell structure repeat
for i=1:numel(A.params)
    % find multiindex from linear index
    Aidx = cell(1,numel(cellsize));
    [Aidx{:}] = ind2sub(cellsize,i);
    
    % create transposed index to place in correct location
    % 2 and 3 are swapped
    Atidx = cellfun(@(x) x+(x==2)-(x==3), Aidx, UniformOutput=false);
    % note, if x==1, then x + (x==2) - (x==3) = 1
    %          x==2, then x + (x==2) - (x==3) = 3
    %          x==3, then x + (x==2) - (x==3) = 2
    
    % transpose Parameter
    tmp = A.params{i}';
    
    % for each spatial variable si, check if we need to swap
    % si and si_dum 
    for j=1:nvars
        % dimension along si is 1, 
        % must be multiplier or full integral, swap (si,si_dum)
        if cellsize(j)==1   
            tmp = var_swap(tmp,varsMain{j},varsDummy{j});
        elseif any(Atidx{j}==[2,3]) % must be semisep term, swap (si,si_dum)
            tmp = var_swap(tmp,varsMain{j},varsDummy{j});
        else 
            % multiplier term along si, appearing in both domain and range
            % do not swap
        end
    end
    % place the parameter in the transposed location
    At.param{Atidx{:}} = tmp;
end
end
