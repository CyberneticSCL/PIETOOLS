function P_reshaped = reshape(P1, M, N)
% reshape  reshape of P1 object.
%   
%   B = RESHAPE(A,M,N) or RESHAPE(A,[M,N]) returns the M-by-N quadpoly whose
%   elements are taken columnwise from A. A must have M*N elements.
% Input P1 quadPoly obj (n times n)
% Output sum_i obj(i, i)
 

if isa(P1, 'double')
    P_reshaped = reshape(P1, M, N);
    return
end

if ~isa(P1, 'quadPoly')
    error('quadPoly:reshape:badType', 'reshape supports quadPoly or numeric inputs.');
end


sz = P1.dim;

if sz(1)*sz(2) ~= M*N
    error('quadPoly:reshape:dimMismatch', 'Number of elements must not change.');
end

% extract monomials and coefficients.
Zt = P1.Zt;
Zs = P1.Zs;
nt = P1.nt;
ns = P1.ns;
C = P1.C;


n_Zt = 1;% length of right monomials
for idx = 1:size(Zt,2)
    n_Zt = n_Zt*size(Zt{idx}, 1);
end

n_Zs = 1;% length of left monomials
for idx = 1:size(Zs,2)
    n_Zs = n_Zs*size(Zs{idx}, 1);
end

[i, j, v] = find(C);   % find nonzero entries

% Within-block indices and block indices
a  = mod(i-1, n_Zs) + 1;        % 1..n_Zs
ib = floor((i-1) / n_Zs);       % 0..m-1

b  = mod(j-1, n_Zt) + 1;        % 1..n_Zt
jb = floor((j-1) / n_Zt);       % 0..n-1

full_index = ib + jb*sz(1); % numbering of blocks

iblc2 = mod(full_index, M); %new row indx of blocks
jblc2 = floor(full_index/M); %new col indx of blocks

i2 = iblc2*n_Zs + a;             % rows: n blocks of size ds
j2 = jblc2*n_Zt + b;             % cols: m blocks of size dt

Cnew = sparse(i2, j2, v, M*n_Zs, N*n_Zt);

% % uncomment for testing Cnew
% % create 4D array of C indeces
% C_ND_array = reshape(full(C), n_Zs, sz(1), n_Zt, sz(2)); 
% % C_ND_array(:, i, :, j) is coefficients for P1(i, j).C
% 
% C_ND_array = permute(C_ND_array, [1, 3, 2, 4]); %change dimensions now
% % C_ND_array(:, :, i, j) is coefficients for P1(i, j).C
% 
% %Now reshape coefficients
% C_ND_array = reshape(C_ND_array, n_Zs, n_Zt, M, N);
% 
% %now permute indeces back
% C_ND_array = ipermute(C_ND_array, [1, 3, 2, 4]);
% 
% %finally, compute the final C matrix
% 
% C_reshaped = reshape(C_ND_array, n_Zs*M, n_Zt*N);
% 
% % C_reshaped is block reshape of C matrix without using sparce str

% define reshaped quadpoly 
P_reshaped = quadPoly(Cnew, Zs, Zt, [M, N], ns, nt);

end

