% Test_tensopvar_class.m -- tests the methods of the 'tensopvar' class
% that have no test of their own: the constructor, initialize,
% size/matdim/deps, num_factors, is_id_op, common_vars, repmat, repelem,
% apply_op (not nested), and the semantics of mtimes, plus, minus and
% uminus.
%
% The expected values are computed independently of the class, from
% operators acting on constant inputs: factors that are multipliers from
% R, Rij = Mij(s) (explicit polynomial matrices, no input variables),
% applied to constant vectors xj give the polynomials Mij*xj, combined as
%   - Kronecker rows (otype(1)=2):   Y = sum_i kron(Mi1*x1,...,Mid*xd),
%   - pointwise rows (otype(1)=1):   Y = sum_i (Mi1*x1).*...*(Mid*xd),
% (apply_op renames the variable s of the output of a factor to s_k, k the
% degree of s up to that factor, for Kronecker rows; removed here by
% stripping '_1', or checked explicitly). Operators on a function space
% (multipliers on L2) are checked against products of polynomials. repmat
% and repelem are checked against the definitions kron(ones(r1,r2),T) and
% kron(T,ones(r1,r2)), acting on inputs kron(v,x1) and kron(xd,v), for
% which the output is sum(v)*kron(ones(r1,1),Y) and
% sum(v)*kron(Y,ones(r1,1)).
% Parts:
%   (1) constructor;
%   (2) initialize: derived fields and error conditions;
%   (3) size, matdim, deps;
%   (4) num_factors;
%   (5) is_id_op;
%   (6) common_vars;
%   (7) apply_op;
%   (8) repmat and repelem;
%   (9) mtimes;
%   (10) plus, minus, uminus.
%
% DJ (CC), 10/15/2026: Initial coding;

tol = 1e-9;
pvar s1 s2 s3 s1_1 s1_2
PV = {s1,s2,s3};        PN = {'s1','s2','s3'};
dom1 = [0,1];
% multiplier operator between spaces with the variables vin and vout (cellstr), of the matrix M
mkop = @(M,vin,vout) mat2tensopvar(M,struct('in',{vin},'out',{vout}),struct('in',repmat(dom1,numel(vin),1),'out',repmat(dom1,numel(vout),1)));
mkR = @(M,vout) mkop(M,cell(1,0),vout);                 % from R to L2[0,1]^numel(vout)
mkF = @(M,v) mkop(M,v,v);                               % multiplier on L2[0,1]^numel(v)
sop = @(T) T.ops{1,1};                                  % the 'sopvar' defining a single-factor operator

%% (1) constructor
T = tensopvar();
chk(isa(T,'tensopvar') && isempty(T.ops) && isequal(T.otype,[2,2]) && isequal(size(T),[0,0]),'(1) empty constructor');
S = mat2sopvar([1+s1;2],struct('in',{cell(1,0)},'out',{{'s1'}}),struct('in',zeros(0,2),'out',dom1));
T = tensopvar(S);
chk(isa(T,'tensopvar') && isequal(size(T),[2,1]) && isequal(size(T.ops),[1,1]) && isa(T.ops{1,1},'sopvar'),'(1) from sopvar');
T = tensopvar({S});
chk(isa(T,'tensopvar') && isequal(size(T),[2,1]),'(1) from cell');
expect_error(@() tensopvar(5),'Input must');
expect_error(@() tensopvar(S,S),'Too many');
tensopvar TA TB             % command syntax declares tensopvar objects
chk(isa(TA,'tensopvar') && isa(TB,'tensopvar') && isempty(TA.ops),'(1) command syntax: variables should be tensopvar objects');

%% (2) initialize: derived fields
for rep=1:12
    rng(rep);
    d = 1+mod(rep,3);       N = 1+mod(rep,2);       M = 1+mod(floor(rep/2),3);      nterms = 1+mod(rep,2);
    type = [1+mod(rep,2),2];
    dims = randi(3,d,2);    if type(1)==1,  dims(:,1) = dims(1,1);  end
    vars = struct('in',{arrayfun(@(i) sprintf('x%d',i),1:N,'UniformOutput',false)},'out',{arrayfun(@(i) sprintf('y%d',i),1:M,'UniformOutput',false)});
    dom = struct('in',[0,1;-1,2;0,3],'out',[0,2;1,4;-1,0]);        dom.in = dom.in(1:N,:);     dom.out = dom.out(1:M,:);
    depmat = struct('in',rand(d,N)>0.4,'out',rand(d,M)>0.4);
    depmat.in(1,:) = true;      depmat.out(1,:) = true;         % every variable appears in some factor
    Top = rand_tensopvar(dims,vars,dom,depmat,nterms,type,2,0.5);
    chk(isequal(Top.vars.in,vars.in) && isequal(Top.vars.out,vars.out),sprintf('(2.%d) vars',rep));
    chk(isequal(Top.dom.in,dom.in) && isequal(Top.dom.out,dom.out),sprintf('(2.%d) dom',rep));
    chk(isequal(Top.dims,dims),sprintf('(2.%d) dims',rep));
    chk(isequal(Top.depmat.in,double(depmat.in)) && isequal(Top.depmat.out,double(depmat.out)),sprintf('(2.%d) depmat',rep));
    Top2 = initialize(Top);
    chk(isequal(Top2.dims,Top.dims) && isequal(Top2.vars,Top.vars) && isequal(Top2.depmat,Top.depmat),sprintf('(2.%d) initialize idempotent',rep));
end
% error conditions
expect_error(@() initialize(set_ops({5},[2,2])),'must be of type');
expect_error(@() initialize(set_ops({sop(mkop(eye(2),{'s1'},{'s1'}));sop(mkop(eye(3),{'s1'},{'s1'}))},[2,2])),'Matrix dimensions');
expect_error(@() initialize(set_ops({sop(mkR(ones(2,1),{'s1'}));sop(mkR(ones(3,1),{'s1'}))},[2,2])),'Matrix dimensions');     % rows only
expect_error(@() initialize(set_ops({sop(mkR(ones(1,2),{'s1'}));sop(mkR(ones(1,3),{'s1'}))},[2,2])),'Matrix dimensions');     % columns only
expect_error(@() initialize(set_ops({sop(mkop(1,{'s1'},{'s1'}));sop(mkop(1,{'s2'},{'s2'}))},[2,2])),'variables of each term');
Sb = mat2tensopvar(1,struct('in',{{'s1'}},'out',{{'s1'}}),struct('in',[0,2],'out',[0,2]));       % s1 on [0,2]
expect_error(@() initialize(set_ops({sop(mkop(1,{'s1'},{'s1'})),sop(Sb)},[2,2])),'same domains');
expect_error(@() initialize(set_ops({sop(mkR(ones(2,1),{'s1'})),sop(mkR(ones(3,1),{'s1'}))},[1,2])),'Row dimensions');
expect_error(@() initialize(set_ops({sop(mkR(ones(1,2),{'s1'})),sop(mkR(ones(1,3),{'s1'}))},[2,1])),'Column dimensions');
% the same factor on different variables in different terms
A1 = sop(mkR(1,{'s1'}));        A2 = sop(mkR(1,{'s2'}));
expect_error(@() initialize(set_ops({A1,A2;A2,A1},[2,2])),'same factor in different terms');
ok_T = initialize(set_ops({A1,A2;A1,A2},[2,2]));          % same structure in each term is fine
chk(isequal(size(ok_T.ops),[2,2]) && isequal(ok_T.vars.out,{'s1','s2'}),'(2) valid two-term operator');

%% (3) size, matdim, deps
for rep=1:16
    rng(100+rep);
    d = 1+mod(rep,3);       type = [1+mod(rep,2),2];
    dims = randi(3,d,2);    if type(1)==1,  dims(:,1) = dims(1,1);  end
    vars = struct('in',{{'x1','x2'}},'out',{{'y1','y2'}});
    depmat = struct('in',rand(d,2)>0.3,'out',rand(d,2)>0.3);
    depmat.in(1,:) = true;      depmat.out(1,:) = true;
    Top = rand_tensopvar(dims,vars,[0,1],depmat,1,type,2,0.5);
    if type(1)==1,  m = dims(1,1);  else,  m = prod(dims(:,1));     end
    n = prod(dims(:,2));
    [m1,n1] = size(Top);
    chk(isequal(size(Top),[m,n]) && isequal(Top.matdim,[m,n]) && m1==m && n1==n && size(Top,1)==m && size(Top,2)==n,sprintf('(3.%d) size',rep));
    if type(1)==1,  degout = max(depmat.out,[],1);  else,  degout = sum(depmat.out,1);  end
    chk(isequal(Top.deps.out,double(degout)) && isequal(Top.deps.in,double(sum(depmat.in,1))),sprintf('(3.%d) deps',rep));
end
expect_error(@() size(Top,3),'Dimension must be');

%% (4) num_factors
for d=1:4
    Top = set_ops(repmat({sop(mkR(1,{'s1'}))},1,d),[2,2]);      Top = initialize(Top);
    chk(num_factors(Top)==d,sprintf('(4) num_factors d=%d',d));
end
% (pointwise rows: the sum of operators with different dependence on the output variables requires equal deps.out)
Ta = initialize(set_ops({sop(mkR([1+s1;2],{'s1'})),sop(mkR([s2,3;1,s2],{'s2'}))},[1,2]));
Tb = initialize(set_ops({sop(mkR([s1;s1^2],{'s1','s2'})),sop(mkR([2,s1;s1,1],{'s1','s2'}))},[1,2]));
Tab = Ta+Tb;
chk(isequal(size(Tab.ops),[2,1]) && isa(Tab.ops{1,1},'tensopvar') && num_factors(Tab)==2,'(4) nested sum has 2 leaf factors');
Tabb = Tab+Ta;
chk(num_factors(Tabb)==2 && isa(Tabb.ops{1,1},'tensopvar'),'(4) nested of nested');

%% (5) is_id_op
id1 = mkF(eye(2),{'s1'});       id2 = mkF(1,{'s2'});        id3 = mkF(1,{'s3'});
chk(is_id_op(id1),'(5) identity on L2^2[0,1]');
chk(is_id_op(mkF(1,{'s1'})),'(5) scalar identity');
chk(~is_id_op(mkF(2*eye(2),{'s1'})),'(5) 2*I');
chk(~is_id_op(mkF(eye(2)+[0,s1;0,0],{'s1'})),'(5) I + multiplier');
chk(~is_id_op(mkF(zeros(2),{'s1'})),'(5) zero');
chk(~is_id_op(mkF(ones(2,3),{'s1'})),'(5) not square');
chk(~is_id_op(mkR(1,{'s1'})),'(5) from R: input space differs from the output space');
chk(is_id_op(mkF((1+1e-14)*eye(2),{'s1'})),'(5) identity up to the default tolerance');
chk(~is_id_op(mkF((1+1e-6)*eye(2),{'s1'})),'(5) not identity above the default tolerance');
chk(is_id_op(mkF((1+1e-6)*eye(2),{'s1'}),1e-3),'(5) identity with a larger tolerance');
Id12 = otimes(id1,id2,[2,2]);       Id123 = otimes(Id12,id3,[2,2]);
chk(is_id_op(Id12) && is_id_op(Id123),'(5) tensor product of identities');
chk(~is_id_op(otimes(id1,mkF(2,{'s2'}),[2,2])),'(5) tensor product with 2*I');
Idp = otimes(mkF(1,{'s1'}),mkF(1,{'s2'}),[1,2]);
chk(~is_id_op(Idp),'(5) pointwise rows and Kronecker columns are not an identity');
Z12 = otimes(mkF(zeros(2),{'s1'}),id2,[2,2]);
chk(is_id_op(Id12+Z12),'(5) identity plus a zero term: the zero term is dropped');
chk(~is_id_op(Id12+Id12),'(5) two nonzero terms');

%% (6) common_vars
A = otimes(mkF(1+s1,{'s1'}),mkF(2,{'s2'}),[2,2]);          B = otimes(mkF(s3,{'s3'}),mkF(1,{'s3'}),[2,2]);
[A2,B2] = common_vars(A,B);
chk(isequal(A2.vars.in,{'s1','s2','s3'}) && isequal(B2.vars.in,{'s1','s2','s3'}) && isequal(A2.vars.out,B2.vars.out),'(6) merged variables');
chk(isequal(A2.dom.in,B2.dom.in) && isequal(size(A2.dom.in),[3,2]),'(6) merged domains');
chk(isequal(A2.depmat.in,[A.depmat.in,zeros(2,1)]) && isequal(B2.depmat.in,[zeros(2,2),B.depmat.in]),'(6) depmat: disjoint variables');
chk(isequal(A2.depmat.out,A2.depmat.in) && isequal(size(B2.depmat.out),[2,3]),'(6) depmat out');
% overlapping variables
A = otimes(mkF(1,{'s1','s2'}),mkF(1,{'s2'}),[2,2]);         B = otimes(mkF(1,{'s2','s3'}),mkF(1,{'s3'}),[2,2]);
[A2,B2] = common_vars(A,B);
chk(isequal(A2.vars.out,{'s1','s2','s3'}) && isequal(A2.depmat.out,[1,1,0;0,1,0]) && isequal(B2.depmat.out,[0,1,1;0,0,1]),'(6) overlapping variables');
% identical variables
[A2,B2] = common_vars(A,A);
chk(isequal(A2.depmat,A.depmat) && isequal(B2.vars,A.vars),'(6) identical variables');
% same variable, different domain
expect_error(@() common_vars(mkF(1,{'s1'}),Sb),'same domain');

%% (7) apply_op
% (7a) operators from R, Kronecker rows, factors in different variables, random terms and dimensions
for rep=1:18
    rng(200+rep);
    d = 1+mod(rep,3);       K = 1+mod(floor(rep/3),3);
    nr = randi(3,d,1);      nc = randi(3,d,1);
    Mij = cell(K,d);
    for i=1:K
        for j=1:d,  Mij{i,j} = randpoly(nr(j),nc(j),PV(j),2);   end
    end
    Top = build_R(Mij,cellfun(@(c) {c},PN(1:d),'UniformOutput',false),[2,2]);
    xs = arrayfun(@(j) randn(nc(j),1),1:d,'UniformOutput',false);
    Y = strip(apply_op(Top,xs));
    close_to(Y,ref_R(Mij,xs,[2,2]),tol,sprintf('(7a.%d)',rep));
end
% (7b) pointwise rows: all factors with the same number of rows, random dependence on s1,s2,s3
for rep=1:12
    rng(300+rep);
    d = 1+mod(rep,3);       K = 1+mod(floor(rep/3),2);      m = randi(3);
    nc = randi(3,d,1);
    Mij = cell(K,d);        vars = cell(1,d);
    for j=1:d
        idx = find(rand(1,3)>0.4);     if isempty(idx),  idx = j;  end
        vars{j} = PN(idx);
        for i=1:K,  Mij{i,j} = randpoly(m,nc(j),PV(idx),2);   end
    end
    Top = build_R(Mij,vars,[1,2]);
    xs = arrayfun(@(j) randn(nc(j),1),1:d,'UniformOutput',false);
    Y = apply_op(Top,xs);
    close_to(Y,ref_R(Mij,xs,[1,2]),tol,sprintf('(7b.%d)',rep));
end
% (7c) Kronecker product of factors in the same variable: the variables of the outputs are s1_1, s1_2
Top = build_R({1+s1,s1^2},{{'s1'},{'s1'}},[2,2]);
x1 = 2;     x2 = 3;
close_to(apply_op(Top,{x1,x2}),(1+s1_1)*s1_2^2*x1*x2,tol,'(7c) same variable');
close_to(apply_op(Top,{polynomial(x1),polynomial(x2)}),(1+s1_1)*s1_2^2*x1*x2,tol,'(7c) polynomial inputs');
% (7d) a single factor takes the input directly, not in a cell
Top = build_R({[1+s1,s1^2;2,s1]},{{'s1'}},[2,2]);       x = [1;-2];
close_to(strip(apply_op(Top,x)),[1+s1,s1^2;2,s1]*x,tol,'(7d) bare input');
% (7e) operators on a function space: multipliers on L2
pv1 = polynomial({'s1'}');
M1 = [1+s1,s1^2;2,3-s1];
Top = mkF(M1,{'s1'});           x = rand_poly([2,1],pv1,2,5);
close_to(strip(apply_op(Top,{x})),M1*x,tol,'(7e) multiplier on L2');
% one factor from R, one multiplier on L2, Kronecker rows
Mr = [1+s1;2];      Mf = [s2,1;2,s2^2];
Top = initialize(set_ops({sop(mkR(Mr,{'s1'})),sop(mkF(Mf,{'s2'}))},[2,2]));
xr = 3;     xf = rand_poly([2,1],polynomial({'s2'}'),2,5);
close_to(strip(apply_op(Top,{xr,xf})),kron(Mr*xr,Mf*xf),tol,'(7e) R factor and L2 factor');
% the same with pointwise rows: both factors with 2 rows
Top = initialize(set_ops({sop(mkR(Mr,{'s1'})),sop(mkF(Mf,{'s2'}))},[1,2]));
close_to(apply_op(Top,{xr,xf}),(Mr*xr).*(Mf*xf),tol,'(7e) pointwise rows');
% (7f) an operator that has not been initialized
Top = set_ops({sop(mkR(Mr,{'s1'}))},[2,2]);
close_to(strip(apply_op(Top,xr)),Mr*xr,tol,'(7f) not initialized');
% (7g) errors
Top = build_R({1+s1,s2},{{'s1'},{'s2'}},[2,2]);
expect_error(@() apply_op(Top,1),'1xd cell');
expect_error(@() apply_op(Top,{1,1,1}),'1xd cell');
expect_error(@() apply_op(Top,{[1;2],1}),'Column dimension');
expect_error(@() apply_op(Top,{1,'a'}),'polynomial');
Tpc = Top;      Tpc.otype = [2,1];
expect_error(@() apply_op(Tpc,{1,1}),'Pointwise products');

%% (8) repmat and repelem
for rep=1:24
    rng(400+rep);
    type = [1+mod(rep,2),2];
    d = 1+mod(floor(rep/2),3);          K = 1+mod(floor(rep/6),2);
    r = randi(3,1,2);                   if rep<=6,  r = [1,1]+(rep>3)*[0,1];  end
    m = randi(3);       nr = randi(3,d,1);      if type(1)==1,  nr(:) = m;  end
    nc = randi(3,d,1);
    Mij = cell(K,d);        vars = cell(1,d);
    for j=1:d
        if type(1)==1,  idx = 1+mod(j+rep,3);   else,  idx = j;     end       % Kronecker: factor j on s_j only, so '_1' names
        vars{j} = PN(idx);
        for i=1:K,  Mij{i,j} = randpoly(nr(j),nc(j),PV(idx),2);   end
    end
    Top = build_R(Mij,vars,type);
    xs = arrayfun(@(j) randn(nc(j),1),1:d,'UniformOutput',false);
    Y = apply_op(Top,xs);
    v = randn(r(2),1);
    lbl = sprintf('(8.%d) type [%d,2], d=%d, K=%d, r=[%d,%d]',rep,type(1),d,K,r);
    % repmat: kron(ones(r1,r2),T) acts on kron(v,x1) o x2 o ... o xd
    Tr = repmat(Top,r);
    xr = xs;    xr{1} = kron(v,xs{1});
    close_to(strip_if(apply_op(Tr,xr),type),sum(v)*kron(ones(r(1),1),strip_if(Y,type)),tol,['repmat ' lbl]);
    chk(isequal(size(Tr),[r(1)*size(Top,1),r(2)*size(Top,2)]),['size of repmat ' lbl]);
    chk(isequal(size(repmat(Top,r(1),r(2))),size(Tr)),['repmat(A,r1,r2) ' lbl]);
    % repelem: kron(T,ones(r1,r2)) acts on x1 o ... o kron(xd,v)
    Te = repelem(Top,r);
    xe = xs;    xe{d} = kron(xs{d},v);
    close_to(strip_if(apply_op(Te,xe),type),sum(v)*kron(strip_if(Y,type),ones(r(1),1)),tol,['repelem ' lbl]);
    chk(isequal(size(Te),[r(1)*size(Top,1),r(2)*size(Top,2)]),['size of repelem ' lbl]);
    chk(isequal(size(repelem(Top,r(1),r(2))),size(Te)),['repelem(A,r1,r2) ' lbl]);
end
Top = build_R({1+s1,s2},{{'s1'},{'s2'}},[2,2]);
for fn = {@repmat,@repelem}
    f = fn{1};
    expect_error(@() f(Top),'Not enough');
    expect_error(@() f(Top,[1,2,3]),'1x2 vector');
    expect_error(@() f(Top,-1,2),'nonnegative');
    expect_error(@() f(Top,1.5,2),'nonnegative');
    expect_error(@() f(Top,1,2,3));
end

%% (9) mtimes
rng(500);
Mij = {randpoly(2,3,{s1},2),randpoly(1,2,{s2},2);randpoly(2,3,{s1},2),randpoly(1,2,{s2},2)};
Top = build_R(Mij,{{'s1'},{'s2'}},[2,2]);
xs = {randn(3,1),randn(2,1)};
Y0 = apply_op(Top,xs);
% scalar multipliers
close_to(apply_op(3*Top,xs),3*Y0,tol,'(9) 3*T');
close_to(apply_op(Top*3,xs),3*Y0,tol,'(9) T*3');
close_to(apply_op(1*Top,xs),Y0,tol,'(9) 1*T');
% non-identity 'polynomial' multipliers are not supported (@sopvar/mtimes takes numeric multipliers only)
expect_error(@() polynomial(3)*Top,'only supported if the object is the identity');
expect_error(@() (1+s1)*Top,'only supported if the object is the identity');
% 'dpvar' multipliers are supported (@sopvar/mtimes routes them to dpvar_op_copvar)
dpvar da
Td = da*Top;
chk(isa(Td,'tensopvar') && isequal(size(Td),size(Top)),'(9) dpvar multiplier');
close_to(apply_op(polynomial(1)*Top,xs),Y0,tol,'(9) polynomial(1)*T');
% matrix multipliers, single factor
Ts = build_R({randpoly(3,2,{s1},2)},{{'s1'}},[2,2]);
xs1 = randn(2,1);       Ys = apply_op(Ts,xs1);
Afun = @(v) [1+v, v^2, 2; 3, v, 1-v];
Amat = [2,1,3;1,-1,2];      % numeric matrix (a 'polynomial' matrix is not supported, see above)
expect_error(@() Afun(s1)*Ts,'only supported if the object is the identity');
chk(isequal(size(Amat*Ts),[2,2]) && isequal(size(Ts*xs1(1)),size(Ts)),'(9) size of matrix*T');
close_to(apply_op(Amat*Ts,xs1),Amat*Ys,tol,'(9) matrix*T, single factor');
close_to(apply_op(eye(3)*Ts,xs1),Ys,tol,'(9) eye*T');
close_to(apply_op(polynomial(eye(3))*Ts,xs1),Ys,tol,'(9) polynomial eye*T');
% identity operator on L2^3[0,1]
close_to(apply_op(mkF(eye(3),{'s1'})*Ts,xs1),Ys,tol,'(9) identity operator*T');
% composition: A multiplier on L2 in each factor, B from R, mixed product property
for KA=1:2
    for KB=1:2
        rng(510+10*KA+KB);
        MA = cell(KA,2);    MB = cell(KB,2);
        for i=1:KA,     MA{i,1} = randpoly(2,2,{s1},2);     MA{i,2} = randpoly(1,1,{s2},2);     end
        for i=1:KB,     MB{i,1} = randpoly(2,3,{s1},2);     MB{i,2} = randpoly(1,2,{s2},2);     end
        TA = set_ops(cell(KA,2),[2,2]);
        for i=1:KA,     TA.ops{i,1} = sop(mkF(MA{i,1},{'s1'}));     TA.ops{i,2} = sop(mkF(MA{i,2},{'s2'}));     end
        TA = initialize(TA);
        TB = build_R(MB,{{'s1'},{'s2'}},[2,2]);
        xs = {randn(3,1),randn(2,1)};
        Yref = 0;
        for i=1:KA
            for k=1:KB
                Yref = Yref + kron(MA{i,1}*MB{k,1}*xs{1},MA{i,2}*MB{k,2}*xs{2});
            end
        end
        C = TA*TB;
        chk(isequal(size(C.ops),[KA*KB,2]) && isequal(size(C),[2,6]),sprintf('(9) composition structure, KA=%d, KB=%d',KA,KB));
        close_to(strip(apply_op(C,xs)),Yref,tol,sprintf('(9) composition, KA=%d, KB=%d',KA,KB));
    end
end
% errors
expect_error(@() TA*build_R({randpoly(2,2,{s1},2)},{{'s1'}},[2,2]),'factor');     % 2 factors composed with 1
expect_error(@() 'a'*Top,'not supported');
expect_error(@() Top*'a','not supported');
expect_error(@() eye(3)*Top,'Inner dimensions');
expect_error(@() Afun(s1)*Top,'Inner dimensions');
Tm = build_R({randpoly(3,2,{s1},2),randpoly(3,2,{s2},2)},{{'s1'},{'s2'}},[2,2]);
expect_error(@() randn(2,9)*Tm,'Left-multiplication');          % non-scalar matrix against two factors
Tcomp1 = initialize(set_ops({sop(mkF(1+s1,{'s1'})),sop(mkF(s2,{'s2'}))},[2,2]));
Tcomp2 = build_R({1,1},{{'s1'},{'s2'}},[2,2]);          Tcomp2b = Tcomp2;       Tcomp2b.otype = [1,2];
expect_error(@() Tcomp1*Tcomp2b,'Kronecker');                  % pointwise rows in B cannot be composed factor by factor
Tcomp3 = initialize(set_ops({sop(mkF(1+s1,{'s1'}))},[2,2]));
expect_error(@() Tcomp3*Tcomp2,'factor');                     % 1 factor with 2 factors

%% (10) plus, minus, uminus
for rep=1:8
    rng(600+rep);
    K1 = 1+mod(rep,2);      K2 = 1+mod(floor(rep/2),2);
    sameStructure = mod(rep,4)<2;
    M1 = cell(K1,2);        M2 = cell(K2,2);
    for i=1:K1,  M1{i,1} = randpoly(2,2,{s1},2);     M1{i,2} = randpoly(2,2,{s2},2);     end
    for i=1:K2
        if sameStructure
            M2{i,1} = randpoly(2,2,{s1},2);     M2{i,2} = randpoly(2,2,{s2},2);
        else
            M2{i,1} = randpoly(2,2,{s1,s2},2);  M2{i,2} = randpoly(2,2,{s1,s2},2);          % other dependence: nested sum
        end
    end
    % nested sums need pointwise rows (equal deps.out for factors with different dependence), see (4)
    type = [2-(~sameStructure),2];
    T1 = build_R(M1,{{'s1'},{'s2'}},type);
    if sameStructure
        T2 = build_R(M2,{{'s1'},{'s2'}},type);
    else
        T2 = build_R(M2,{{'s1','s2'},{'s1','s2'}},type);
    end
    xs = {randn(2,1),randn(2,1)};
    Y1 = strip(apply_op(T1,xs));        Y2 = strip(apply_op(T2,xs));
    Tp = T1+T2;     Tm = T1-T2;     Tu = -T1;
    chk(isa(Tp.ops{1,1},'tensopvar')==~sameStructure,sprintf('(10.%d) structure of the sum',rep));
    close_to(strip(apply_op(Tp,xs)),Y1+Y2,tol,sprintf('(10.%d) T1+T2',rep));
    close_to(strip(apply_op(Tm,xs)),Y1-Y2,tol,sprintf('(10.%d) T1-T2',rep));
    close_to(strip(apply_op(Tu,xs)),-Y1,tol,sprintf('(10.%d) -T1',rep));
    close_to(strip(apply_op(T1+0,xs)),Y1,tol,sprintf('(10.%d) T1+0',rep));
    close_to(strip(apply_op(0+T1,xs)),Y1,tol,sprintf('(10.%d) 0+T1',rep));
end
expect_error(@() build_R({1},{{'s1'}},[2,2])+build_R({[1,2]},{{'s1'}},[2,2]),'dimensions');
expect_error(@() build_R({1},{{'s1'}},[2,2])+build_R({1},{{'s2'}},[2,2]),'variables');
expect_error(@() build_R({1},{{'s1'}},[2,2])+'a','not be converted');

fprintf('Test_tensopvar_class passed (%d cases).\n',counter(0))


%% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %%
function n = counter(inc)
% number of checks performed
persistent c
if isempty(c),  c = 0;  end
c = c+inc;
n = c;
end

function chk(cond,label)
counter(1);
if ~cond
    error("Test_tensopvar_class failed %s.",label)
end
end

function expect_error(fun,msg)
counter(1);
try
    out = fun();            % with an output, so that tensopvar(...) is not read as command syntax (nargout=0) %#ok<NASGU>
catch ME
    if nargin>=2 && ~contains(ME.message,msg)
        error("Test_tensopvar_class failed: error '%s' does not contain '%s'.",ME.message,msg)
    end
    return
end
error("Test_tensopvar_class failed: no error for %s.",func2str(fun))
end

function T = set_ops(ops,type)
% 'tensopvar' with the cell ops and otype, not initialized
T = tensopvar();
T.ops = ops;
T.otype = type;
end

function M = randpoly(m,n,vars,deg)
% random m x n polynomial matrix in the pvars in the cell vars, of degree deg
M = polynomial(randn(m,n));
for k=1:numel(vars)
    for p=1:deg
        M = M + randn(m,n)*vars{k}^p;
    end
end
end

function T = build_R(Mij,vars,type)
% K x d cell Mij of polynomial matrices to the 'tensopvar' sum_i Ri1 o ... o Rid of multipliers from R,
% with factor j in the variables vars{j} (cellstr)
[K,d] = size(Mij);
T = tensopvar();        T.otype = type;     T.ops = cell(K,d);
for i=1:K
    for j=1:d
        vj = vars{j};
        Tij = mat2tensopvar(Mij{i,j},struct('in',{cell(1,0)},'out',{vj}),struct('in',zeros(0,2),'out',repmat([0,1],numel(vj),1)));
        T.ops{i,j} = Tij.ops{1,1};
    end
end
T = initialize(T);
end

function Y = ref_R(Mij,xs,type)
% expected action sum_i kron_j(Mij*xj) (type(1)=2) or sum_i prod_j (Mij*xj) (type(1)=1)
[K,d] = size(Mij);
Y = 0;
for i=1:K
    Yi = 1;
    for j=1:d
        if type(1)==2,  Yi = kron(Yi,Mij{i,j}*xs{j});
        else,           Yi = Yi.*(Mij{i,j}*xs{j});
        end
    end
    Y = Y+Yi;
end
end

function Y = strip(Y)
% remove the '_1' that apply_op appends to the variables of the output of a factor
Y = polynomial(Y);
Y.varname = regexprep(Y.varname,'_1$','');
end

function Y = strip_if(Y,type)
% strip for Kronecker rows only (pointwise rows keep the names)
if type(1)==2,  Y = strip(Y);   else,   Y = polynomial(Y);  end
end

function close_to(Y,Yref,tol,label)
counter(1);
Yref = polynomial(Yref);        Y = polynomial(Y);
if ~isequal(size(Y),size(Yref))
    error("Test_tensopvar_class failed %s: size %s, expected %s.",label,mat2str(size(Y)),mat2str(size(Yref)))
end
D = Y-Yref;
err = full(max([0;abs(nonzeros(D.coefficient))]));
if err>tol
    error("Test_tensopvar_class failed %s: residual %.3e.",label,err)
end
if isempty(nonzeros(Yref.coefficient))
    error("Test_tensopvar_class %s: reference value is zero; test would be vacuous.",label)
end
end
