function At = ctranspose(A)
% Pt = ctranspose(P) transposes a sdopvar operator P: L_2^p[S1,S3] to L_2^q[S2,S3]
% Date: 8/17/26
% Compute the number of variables S3 that appear in both the input and output spaces.
nvars3 = numel(intersect(A.vars.in,A.vars.out));
% each common variable has three kernel types: multiplier, lower integral, upper integral.
%1 = multiplier
%2 = lower integral
%3 = upper integral
%Under adjoint/transposition, lower and upper integral kernels swap. Multiplier terms stay multiplier
idx = repmat({[1 3 2]},1,nvars3);
params = A.params; 
% swap the parameters 2 and 3 only if there is at least one common variable and the 
% parameter array has the expected full 3^nvars3 number of cells.
if nvars3>0 && numel(params.A)==3^nvars3
    params.A = params.A(idx{:});
    params.B = params.B(idx{:});
end

function K = transposeVecMap(m,n)
old_idx = reshape(1:m*n,m,n);
new_idx = reshape(1:m*n,n,m).';
K = sparse(new_idx(:),old_idx(:),1,m*n,m*n);
end
old_left_size  = A.dims(1) * prod(cellfun(@numel,A.ZL));
old_right_size = A.dims(2) * prod(cellfun(@numel,A.ZR));
K = transposeVecMap(old_left_size,old_right_size);
% loop over every kernel cell.
for ii=1:numel(params.A)
    % K maps vec(C) to vec(C'). Since B stores vec(C_i)' by rows, multiply
    % by K' on the right to transpose every decision coefficient at once.
    params.A{ii} = K*params.A{ii};
    params.B{ii} = params.B{ii}*K.';
end
At_vars = struct('in',{A.vars.out},'out',{A.vars.in});
At_dom = struct('in',A.dom.out,'out',A.dom.in);
At = sdopvar(params,At_vars,A.Zd,A.ZR,A.ZL,At_dom,[A.dims(2),A.dims(1)]);
end
