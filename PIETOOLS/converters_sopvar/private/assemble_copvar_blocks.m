function P = assemble_copvar_blocks(entries,rmasks,rdims,cmasks,cdims,vars,dom)
% P = ASSEMBLE_COPVAR_BLOCKS(ENTRIES,RMASKS,RDIMS,CMASKS,CDIMS,VARS,DOM)
% assembles a 'copvar' from polynomial kernels, for any number of spatial
% variables.
%
% INPUTS
% - entries:    Struct array with fields
%                   r, c:   binary indices of the output and input space
%                           (sum_j 2^(j-1)*[space contains variable j]);
%                   rows, cols: rows/columns of the block (output and input
%                           space), in which the kernel is placed;
%                   sub:    1 x nvars array, the parameter index (1: multiplier,
%                           2: int_a^s, 3: int_s^b) along each variable; only
%                           the entries for variables shared by both spaces
%                           are used;
%                   C:      'polynomial' matrix, the kernel, using the primary
%                           variables for variables only in the output space,
%                           both primary and dummy for shared variables, and
%                           dummy variables for variables only in the input;
%               Kernels of entries with the same r, c, and sub are summed;
% - rmasks, rdims:  Binary indices and dimensions of the output spaces;
% - cmasks, cdims:  Binary indices and dimensions of the input spaces;
% - vars, dom:      nvars x 2 arrays of the primary and dummy variables, and
%                   nvars x 2 array of the domain of each variable;
%
% OUTPUT
% - P:      'copvar' with one row per output space and one column per input
%           space, in the order given. Blocks without any entries are
%           structurally zero.
%
% PIETOOLS - assemble_copvar_blocks
% If you modify this code, document all changes carefully and include date
% authorship, and a brief description of modifications
%
% Initial coding: Claude 10/28/2026

nvars = size(vars,1);
names = pvar2varname(vars(:,1));    names = names(:)';
[~,order] = sort(names);
nr = numel(rmasks);     nc = numel(cmasks);
params = cell(nr,nc);
for ie=1:numel(entries)
    e = entries(ie);
    i = find(rmasks==e.r);      j = find(cmasks==e.c);
    inr = logical(bitget(e.r,1:nvars));     inc = logical(bitget(e.c,1:nvars));
    S3 = order(inr(order) & inc(order));
    n3 = numel(S3);
    if isempty(params{i,j})
        if n3==0,       sz = [1,1];
        elseif n3==1,   sz = [3,1];
        else,           sz = 3*ones(1,n3);
        end
        K = cell(sz);
        for k=1:numel(K),   K{k} = polynomial(zeros(rdims(i),cdims(j)));    end
        params{i,j} = K;
    end
    K = params{i,j};
    if n3<=1
        lin = e.sub(S3(1:n3));      lin = max([lin,1]);
    else
        subc = num2cell(e.sub(S3));
        lin = sub2ind(3*ones(1,n3),subc{:});
    end
    K{lin}(e.rows,e.cols) = K{lin}(e.rows,e.cols) + e.C;
    params{i,j} = K;
end

reg = names(order);     regdom = dom(order,:);
space_r = false(nr,nvars);      space_c = false(nc,nvars);
for i=1:nr,     sp = logical(bitget(rmasks(i),1:nvars));    space_r(i,:) = sp(order);   end
for j=1:nc,     sp = logical(bitget(cmasks(j),1:nvars));    space_c(j,:) = sp(order);   end
blocks = cell(nr,nc);
for i=1:nr
    for j=1:nc
        if isempty(params{i,j}),    continue,   end
        inr = logical(bitget(rmasks(i),1:nvars));   inc = logical(bitget(cmasks(j),1:nvars));
        S3 = order(inr(order) & inc(order));
        S2 = order(inr(order) & ~inc(order));
        S1 = order(~inr(order) & inc(order));
        blocks{i,j} = sopvar_from_kernels(params{i,j},S2,S3,S1,vars,dom,[rdims(i),cdims(j)]);
    end
end
meta = struct('vars',{reg},'dom',regdom,'space_out',space_r,'space_in',space_c, ...
                'dim_out',rdims(:),'dim_in',cdims(:));
P = copvar(blocks,meta);
end
