% Test_tensopmat_times.m -- verify @tensopmat/times.m: scalar.*tensopmat
% (Case 1, via apply_op), the block/range index arithmetic (Case 2), and
% matrix.*tensopmat for native blocks (Case 3).

% DJ (CC), 10/15/2026: Initial coding;

% NOTES
% - matrix.*tensopmat with 'tensopvar' blocks is not tested: it ends in
%   Aij.*sopvar via @tensopvar/times.m, and '@sopvar' has no 'times.m'
%   (MATLAB's "Invalid data type" error). The same holds for
%   matrix.*tensopvar in @tensopvar/times.m. Case 3 uses native
%   'polynomial' blocks instead.
% - Case 2 replicates times.m's index arithmetic by hand on a random
%   object; it does not call times.

tol = 1e-8;

p = randi([2,4]);  d = randi([2,4]);
N = randi([2,3]);  M = randi([2,3]);
dom = [0,1];

[Cop,info] = rand_tensopmat(p,d,N,M,dom,0.7);
mdim = size(Cop);

X = cell(1,d);
for j=1:d
    X{j} = rand_poly([info.dims_in(j),1],info.vars_in_pool(info.sel_in{j}),1,5);
end
Y = apply_op(Cop,X);

% Case 1: scalar .* tensopmat (the 'all(Adim==1)' branch)
a = 3.25;
Y1 = apply_op(a.*Cop,X);
err1 = cleanpoly(Y1-a*Y,1e-10);
if max(max(abs(err1.C)))>tol
    error("Test_tensopmat_times failed: scalar.*tensopmat mismatch.")
end

% Case 2 (structural): every populated block is visited exactly once and
% mapped to the correct row/column range
nnr_arr = [0;cumsum(Cop.dims.out)];
nnc_arr = [0;cumsum(Cop.dims.in)];
nvisited = 0;
for k=1:numel(Cop.ops)
    [i,j] = ind2sub(size(Cop.ops),k);
    if isempty(Cop.ops{i,j})
        continue
    end
    nvisited = nvisited+1;
    expected_size = [nnr_arr(i+1)-nnr_arr(i), nnc_arr(j+1)-nnc_arr(j)];
    if ~isequal(expected_size,size(Cop.ops{i,j}))
        error("Test_tensopmat_times failed: block (%d,%d) size does not match its row/column range.",i,j)
    end
end
if nvisited~=nnz(~cellfun(@isempty,Cop.ops))
    error("Test_tensopmat_times failed: not every populated block was visited.")
end

% Case 3: matrix .* tensopmat, for objects of blocks of type 'polynomial' acting on R (see `tensopmat`), with
% empty (zero) blocks and blocks of several rows and columns, also if the object has a single row and/or
% column, which is repeated to the size of the matrix (with the dimensions of the object updated). The expected
% value is computed from the action of the object on the unit vectors, which does not use times:
%   ((A.*C)*x)(r) = sum_c A(r,c)*(C*(x_c e_c'))(r'), with c' (r') the column (row) of C, 1 if C has a single one.
names = {'s1','s2'};
for rep_idx=1:30
    p = randi([1,3]);  d = randi([1,3]);
    dout = randi([1,3],p,1);    din = randi([1,3],d,1);
    Cfull = gen_native(dout,din,names,0.5);             % with empty blocks
    Crow = gen_native(1,din,names,0.7);                 % a single row of dimension 1
    Ccol = gen_native(dout,1,names,0.7);                % a single column of dimension 1
    C11 = gen_native(1,1,names,1);
    m2 = randi([2,3]);      n2 = randi([2,3]);
    cases = {randn(size(Cfull)),            Cfull;  ...     % same size, with empty blocks
             randn(m2,size(Crow,2)),        Crow;   ...     % one row, repeated
             randn(size(Ccol,1),n2),        Ccol;   ...     % one column, repeated
             randn(m2,n2),                  C11};           % one row and column, repeated
    for k=1:size(cases,1)
        Amat = cases{k,1};      Cb = cases{k,2};
        Cres = Amat.*Cb;
        [lv,msg] = isvalid(Cres);
        if ~lv || ~isequal(size(Cres),size(Amat))
            error("Test_tensopmat_times failed: matrix.*tensopmat (case 3.%d) is not a valid object of the size of the matrix: %s",k,msg)
        end
        [Y,Yexp] = apply_times_case(Amat,Cb,Cres);
        err = cleanpoly(Y-Yexp,1e-10);
        if max(max(abs(err.C)))>tol
            error("Test_tensopmat_times failed: matrix.*tensopmat (case 3.%d) does not match the action of the object on unit vectors.",k)
        end
    end
end

disp('Test_tensopmat_times passed.')


function C = gen_native(dout,din,names,dens)
% Random 'tensopmat' with blocks of type 'polynomial' in the variables names, of which a fraction dens is
% nonempty, with at least one block in each row and column
p = numel(dout);    d = numel(din);
mask = rand(p,d)<dens;
for i=1:p,  if ~any(mask(i,:)),     mask(i,randi(d)) = true;    end,    end
for j=1:d,  if ~any(mask(:,j)),     mask(randi(p),j) = true;    end,    end
E = cell(p,d);
for i=1:p,  for j=1:d
    if mask(i,j),   E{i,j} = rand_poly([dout(i),din(j)],names(:),1,3) + randn(dout(i),din(j));     end
end,    end
vars = struct('in',{cell(1,0)},'out',{names});
dm = struct('in',zeros(0,2),'out',repmat([0,1],numel(names),1));
C = tensopmat(E,vars,dm);
end

function [Y,Yexp] = apply_times_case(A,C,Cres)
% Applies Cres = A.*C to a random input x in R^n, and computes the expected value from the action of C on unit vectors
[m,n] = size(A);        [mc,nc] = size(C);
din = C.dims.in;        off = [0;cumsum(din)];
x = randn(n,1);
cb = 1:n;       if nc==1,   cb = ones(1,n);     end
rb = 1:m;       if mc==1,   rb = ones(1,m);     end
Yexp = 0;
for c=1:n
    j = find(cb(c)<=off(2:end),1);      e = zeros(din(j),1);        e(cb(c)-off(j)) = 1;
    Xc = cell(1,numel(din));
    for jj=1:numel(din),    Xc{jj} = zeros(din(jj),1);     end
    Xc{j} = e*x(c);
    Yc = apply_op(C,Xc);
    Yexp = Yexp + A(:,c).*Yc(rb(:));
end
% input of the product: blocks of the object, or a single block if the columns were repeated
if nc==1
    Xres = {x};
else
    Xres = cell(1,numel(din));
    for j=1:numel(din),     Xres{j} = x(off(j)+1:off(j+1));     end
end
Y = apply_op(Cres,Xres);
end