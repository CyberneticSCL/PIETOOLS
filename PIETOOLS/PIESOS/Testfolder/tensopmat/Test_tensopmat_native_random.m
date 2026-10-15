% Test_tensopmat_native_random.m -- randomized value test of the
% operations on a 'tensopmat' whose blocks are matrices of type 'double',
% 'polynomial' or 'dpvar' acting on R (see @tensopmat/tensopmat.m), for
% K=1,...,4 spatial variables (with different domains), p x d blocks of
% random sizes, empty blocks, rows that map to R (no variables), and rows
% in different subsets of the variables.
% The expected values are computed from the explicit matrices E{i,j}
% defining the blocks (blockwise matrix algebra, or the sum of E{i,j}*x_j
% for apply_op), not from the class itself, for
%   (1) construction, isvalid, dims, domains; (2) apply_op;
%   (3) uminus, minus, plus; (4) times (scalar and elementwise);
%   (5) mtimes (matrix times object, and object times object with
%   constant blocks); (6) horzcat, vertcat, blkdiag; (7) repmat, repelem
%   (each row and column of the matrix repeated); (8) subsref of matrix
%   indices, subsasgn of a block; (9) clean_op, is_id_op.
% The seed is printed, and part of any error message, so that a failure
% can be reproduced.
%
% DJ (CC), 10/15/2026: Initial coding;

nreps = 100;
tol = 1e-9;
rng('shuffle');         seed = randi(1e6);      rng(seed);
dpvar a b

for rep=1:nreps
    K = randi([1,4]);
    names = arrayfun(@(k) sprintf('x%d',k),1:K,'UniformOutput',false);
    dom = [-(1:K)',(1:K)'+1];                       % a different domain for each variable
    vars = struct('in',{cell(1,0)},'out',{names});
    dm = struct('in',zeros(0,2),'out',dom);
    p = randi([1,3]);   d = randi([1,3]);
    dout = randi([1,3],p,1);        din = randi([1,3],d,1);
    lbl = sprintf('rep %d (K=%d, p=%d, d=%d, seed %d)',rep,K,p,d,seed);
    [E1,C1] = gen_obj(names,vars,dm,dout,din,a,b);

    %% (1) construction
    check_obj(C1,E1,dout,din,names,dom,tol,['(1) ' lbl]);

    %% (2) apply_op, for constant inputs
    X = cell(1,d);      for j=1:d,  X{j} = randn(din(j),1);     end
    Yexp = cell(p,1);
    for i=1:p
        Yi = zeros(dout(i),1);
        for j=1:d,  if ~isempty(E1{i,j}),   Yi = Yi+E1{i,j}*X{j};   end,   end
        Yexp{i} = Yi;
    end
    assert_equal(apply_op(C1,X),vertcat(Yexp{:}),tol,['(2) apply_op ' lbl],true);

    %% (3) uminus, minus, plus, with an object with the same block dimensions
    [E2,C2] = gen_obj(names,vars,dm,dout,din,a,b);
    check_obj(-C1,map1(E1,@(x) -x),dout,din,names,dom,tol,['(3) uminus ' lbl]);
    check_obj(C1+C2,map2(E1,E2,@(x,y) x+y,1),dout,din,names,dom,tol,['(3) plus ' lbl]);
    check_obj(C1-C2,map2(E1,E2,@(x,y) x-y,-1),dout,din,names,dom,tol,['(3) minus ' lbl]);
    check_obj(C2+C1,map2(E2,E1,@(x,y) x+y,1),dout,din,names,dom,tol,['(3) plus, reversed ' lbl]);

    %% (4) times: scalar, and elementwise with a matrix (a function of the variables of the object, or a constant)
    vo = C1.vars.out;
    ri = [0;cumsum(dout)];      cj = [0;cumsum(din)];
    if isempty(vo),     Sc = randn;     Am = randn(sum(dout),sum(din));
    else
        Sc = rand_poly([1,1],vo(:),1,2);        Am = rand_poly([sum(dout),sum(din)],vo(:),1,3);
    end
    check_obj(Sc*C1,map1(E1,@(x) Sc*x),dout,din,names,dom,tol,['(4) scalar times ' lbl]);
    Ea = cell(p,d);
    for i=1:p,  for j=1:d
        if ~isempty(E1{i,j}),   Ea{i,j} = Am(ri(i)+1:ri(i+1),cj(j)+1:cj(j+1)).*E1{i,j};     end
    end,    end
    check_obj(Am.*C1,Ea,dout,din,names,dom,tol,['(4) elementwise times ' lbl]);
    % a single row and/or column, repeated to the size of the matrix: rows/columns of the matrix (dimension 1 blocks)
    [E11,C11] = gen_obj(names,vars,dm,1,1,a,b,true);                % 1 x 1 object, a nonempty block
    m2 = randi([2,3]);      n2 = randi([2,3]);      A2 = randn(m2,n2);
    check_obj(A2.*C11,{A2.*E11{1}},m2,n2,names,dom,tol,['(4) matrix times 1x1 object ' lbl]);

    %% (5) mtimes: matrix times object, and object times object with constant blocks
    r = randi([1,3]);       Dm = randn(r,sum(dout));
    Em = cell(1,d);
    for j=1:d
        Mj = zeros(r,din(j));
        for i=1:p,  if ~isempty(E1{i,j}),   Mj = Mj+Dm(:,ri(i)+1:ri(i+1))*E1{i,j};     end,   end
        Em{1,j} = Mj;
    end
    check_obj(Dm*C1,Em,r,din,names,dom,tol,['(5) matrix times object ' lbl]);
    d2 = randi([1,3]);      din2 = randi([1,3],d2,1);
    Kc = cell(d,d2);        for j=1:d,  for k=1:d2,     Kc{j,k} = randn(din(j),din2(k));    end,    end
    Ep = cell(p,d2);
    for i=1:p,  for k=1:d2
        Pik = zeros(dout(i),din2(k));
        for j=1:d,  if ~isempty(E1{i,j}),   Pik = Pik+E1{i,j}*Kc{j,k};     end,   end
        Ep{i,k} = Pik;
    end,    end
    check_obj(C1*tensopmat(Kc),Ep,dout,din2,names,dom,tol,['(5) object times constants ' lbl]);

    %% (6) horzcat, vertcat, blkdiag
    din3 = randi([1,3],randi([1,2]),1);
    [E3,C3] = gen_obj(names,vars,dm,dout,din3,a,b);                 % same rows as C1
    check_obj([C1,C3],[E1,E3],dout,[din;din3],names,dom,tol,['(6) horzcat ' lbl]);
    dout4 = randi([1,3],randi([1,2]),1);
    [E4,C4] = gen_obj(names,vars,dm,dout4,din,a,b);                 % same columns as C1
    check_obj([C1;C4],[E1;E4],[dout;dout4],din,names,dom,tol,['(6) vertcat ' lbl]);
    dout5 = randi([1,3],randi([1,2]),1);        din5 = randi([1,3],randi([1,2]),1);
    [E5,C5] = gen_obj(names,vars,dm,dout5,din5,a,b);
    Ebd = [E1,cell(p,numel(din5));cell(numel(dout5),d),E5];
    check_obj(blkdiag(C1,C5),Ebd,[dout;dout5],[din;din5],names,dom,tol,['(6) blkdiag ' lbl]);

    %% (7) repmat, repelem
    r1 = randi([1,3]);      r2 = randi([1,3]);
    check_obj(repmat(C1,r1,r2),repmat(E1,r1,r2),repmat(dout,r1,1),repmat(din,r2,1),names,dom,tol,['(7) repmat ' lbl]);
    Erp = cell(p,d);                                                % repeat each row and column of each block
    for i=1:p,  for j=1:d
        if ~isempty(E1{i,j}),   Erp{i,j} = E1{i,j}(repelem(1:dout(i),r1),repelem(1:din(j),r2));     end
    end,    end
    check_obj(repelem(C1,r1,r2),Erp,r1*dout,r2*din,names,dom,tol,['(7) repelem ' lbl]);

    %% (8) subsref of matrix indices (within each block: ascending, as a block is selected by its indices), subsasgn
    nR = sum(dout);     nC = sum(din);
    ridx = find(rand(1,nR)<0.7);        if isempty(ridx),   ridx = randi(nR);   end
    cidx = find(rand(1,nC)<0.7);        if isempty(cidx),   cidx = randi(nC);   end
    Cs = C1(ridx,cidx);
    [Es,dos,dis] = index_blocks(E1,dout,din,ridx,cidx);
    check_obj(Cs,Es,dos,dis,names,dom,tol,['(8) subsref ' lbl],true);       % the indices may only select empty blocks
    [i0,j0] = find(~cellfun('isempty',E1),1);                       % set a nonempty block to a new matrix
    if isempty(C1.vars.out) || ~any(C1.depmat.out(i0,:))
        Pnew = randn(dout(i0),din(j0));
    else
        Pnew = rand_poly([dout(i0),din(j0)],C1.vars.out(logical(C1.depmat.out(i0,:)))',1,3);
    end
    Ca = C1;        Ca{i0,j0} = Pnew;
    Ea = E1;        Ea{i0,j0} = Pnew;
    check_obj(Ca,Ea,dout,din,names,dom,tol,['(8) subsasgn ' lbl]);

    %% (9) clean_op does not change the values, and the identity is recognized
    Cc = clean_op(C1,1e-12);
    check_obj(Cc,E1,dout,din,names,dom,tol,['(9) clean_op ' lbl]);
    n1 = randi([1,3]);
    Cid = tensopmat({eye(n1)});
    if ~is_id_op(Cid) || is_id_op(Cid+Cid) || ~is_id_op(blkdiag(Cid,Cid))
        error("Test_tensopmat_native_random failed (9): identity of %s.",lbl)
    end
end

fprintf('Test_tensopmat_native_random passed (%d repetitions, seed %d).\n',nreps,seed)


%% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %%
function [E,C] = gen_obj(names,vars,dm,dout,din,a,b,isfull)
% Random p x d array E of blocks (matrices of type double, polynomial, dpvar, or empty), and the 'tensopmat' C
% defined by it. Each row of blocks has its own random subset of the variables; some rows map to R
p = numel(dout);    d = numel(din);     K = numel(names);
if nargin<8,    isfull = false;     end
mask = rand(p,d)>0.25 | isfull;
for i=1:p,  if ~any(mask(i,:)),     mask(i,randi(d)) = true;    end,    end
for j=1:d,  if ~any(mask(:,j)),     mask(randi(p),j) = true;    end,    end
E = cell(p,d);
for i=1:p
    isR = rand<0.15;
    rv = names(rand(1,K)<0.6);
    if isempty(rv),     rv = names(randi(K));   end
    for j=1:d
        if ~mask(i,j),  continue,   end
        if isR
            E{i,j} = randn(dout(i),din(j));
        else
            P = rand_poly([dout(i),din(j)],rv(:),2,4) + randn(dout(i),din(j));
            if rand<0.3
                P = P + a*rand_poly([dout(i),din(j)],rv(:),1,3) + b*randn(dout(i),din(j));
            end
            E{i,j} = P;
        end
    end
end
C = tensopmat(E,vars,dm);
end

function F = map1(E,fun)
F = cell(size(E));
for k=1:numel(E),   if ~isempty(E{k}),  F{k} = fun(E{k});   end,   end
end

function F = map2(E1,E2,fun,sgn)
% blockwise fun of two arrays of blocks, where an empty block is zero
F = cell(size(E1));
for k=1:numel(E1)
    if isempty(E1{k}) && isempty(E2{k}),    continue
    elseif isempty(E2{k}),                  F{k} = E1{k};
    elseif isempty(E1{k}),                  F{k} = sgn*E2{k};
    else,                                   F{k} = fun(E1{k},E2{k});
    end
end
end

function [Es,dos,dis] = index_blocks(E,dout,din,ridx,cidx)
% Blocks of the matrix with rows ridx, columns cidx of the block array E (consecutive indices in the same
% block form one block)
bo = cumsum(dout);      bi = cumsum(din);
[rb,rl] = locate(ridx,bo);      [cb,cl] = locate(cidx,bi);
rst = [1,find(diff(rb)~=0)+1];  rend = [rst(2:end)-1,numel(ridx)];
cst = [1,find(diff(cb)~=0)+1];  cend = [cst(2:end)-1,numel(cidx)];
Es = cell(numel(rst),numel(cst));
for i=1:numel(rst)
    for j=1:numel(cst)
        Eij = E{rb(rst(i)),cb(cst(j))};
        if ~isempty(Eij),   Es{i,j} = Eij(rl(rst(i):rend(i)),cl(cst(j):cend(j)));   end
    end
end
dos = (rend-rst+1)';    dis = (cend-cst+1)';
end

function [blk,loc] = locate(idx,bounds)
bounds = bounds(:)';
blk = zeros(size(idx));     loc = zeros(size(idx));     lo = [0,bounds(1:end-1)];
for k=1:numel(idx)
    blk(k) = find(idx(k)<=bounds,1);        loc(k) = idx(k)-lo(blk(k));
end
end

function check_obj(C,E,dout,din,names,dom,tol,label,allow_vacuous)
% C is a valid 'tensopmat' with blocks E (empty blocks of E are zero), dims dout, din, no input variables, and
% output variables among 'names' with their domains; unless allow_vacuous, some block of E must be nonzero
if nargin<9,    allow_vacuous = false;  end
[lv,msg] = isvalid(C);
if ~lv
    error("Test_tensopmat_native_random failed %s: isvalid: %s",label,msg)
end
if ~isequal(C.dims.out(:),dout(:)) || ~isequal(C.dims.in(:),din(:)) || ~isequal(size(C.ops),size(E))
    error("Test_tensopmat_native_random failed %s: dimensions.",label)
end
if ~isempty(C.vars.in) || ~isempty(C.depmat.in) && any(C.depmat.in(:))
    error("Test_tensopmat_native_random failed %s: input variables.",label)
end
[tf,loc] = ismember(C.vars.out,names);
if ~all(tf) || ~isequal(C.dom.out,dom(loc,:)) || ~isequal(C.vars.out,sort(C.vars.out))
    error("Test_tensopmat_native_random failed %s: output variables or domains.",label)
end
% the row of blocks maps to R (depmat.out zero) if no block in it depends on a variable
for i=1:size(E,1)
    if any(C.depmat.out(i,:)>1)
        error("Test_tensopmat_native_random failed %s: depmat.out>1.",label)
    end
end
is_vacuous = true;
for k=1:numel(E)
    if isempty(E{k})
        if ~isempty(C.ops{k}),  assert_equal(C.ops{k},0*C.ops{k},tol,label,true);   end
    else
        if isempty(C.ops{k})
            error("Test_tensopmat_native_random failed %s: block %d is empty.",label,k)
        end
        assert_equal(C.ops{k},E{k},tol,label,true);
        is_vacuous = is_vacuous && max_coeff(E{k})==0;
    end
end
if is_vacuous && ~allow_vacuous
    error("Test_tensopmat_native_random %s: all expected blocks are zero; test would be vacuous.",label)
end
end

function assert_equal(Y,Yref,tol,label,allowzero)
% Y and Yref of type double, polynomial or dpvar of the same size, equal up to tol
if ~isequal(size(Y),size(Yref))
    error("Test_tensopmat_native_random failed %s: size %s, expected %s.",label,mat2str(size(Y)),mat2str(size(Yref)))
end
D = Y-Yref;
err = max_coeff(D);
if err>tol
    error("Test_tensopmat_native_random failed %s: residual %.3e.",label,err)
end
if ~allowzero && max_coeff(Yref)==0
    error("Test_tensopmat_native_random %s: reference value is zero; test would be vacuous.",label)
end
end

function m = max_coeff(D)
if isa(D,'double'),         c = D;
elseif isa(D,'dpvar'),      c = D.C;
else,                       c = D.coefficient;      end
m = full(max([0;abs(nonzeros(c))]));
end
