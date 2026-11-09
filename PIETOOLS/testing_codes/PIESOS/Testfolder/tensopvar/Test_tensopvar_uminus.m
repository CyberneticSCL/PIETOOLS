% Test_tensopvar_uminus.m -- tests @tensopvar/uminus.m. The oracle is
% apply_op(-A,v) against -apply_op(A,v), and apply_op(A+(-A),v) against
% 0, on random inputs v. otype(1) alternates over {1,2}, otype(2)=2.
%
% DJ (CC), 10/15/2026: Initial coding;

ntrials = 20;
for trial=1:ntrials
tol = 1e-8;

nterms = randi(3);
d = randi(5);
N = randi(3);
M = randi(3);
if mod(trial,2)==1
    type = [1,2];
else
    type = [2,2];
end
dnsty = 0.25;
p = 2;
pdegv = 2;

dims = randi(4,[d,2]);
if type(1)==1
    dims(:,1) = dims(1,1);
end
if type(2)==1
    dims(:,2) = dims(1,2);
end
vars = struct();
% disjoint 'x'/'y' name pools: input and output variables must not share a name (rand_sopvar rejects one name with two domains)
vars.in = arrayfun(@(i)sprintf('x%d',i),1:N,'UniformOutput',false);
vars.out = arrayfun(@(i)sprintf('y%d',i),1:M,'UniformOutput',false);
dom = struct();
dom.in = randi([-2,2],[N,1]) +[zeros(N,1),rand([N,1])];
dom.out = randi([-2,2],[M,1]) +[zeros(M,1),rand([M,1])];
degs = struct();
degs.in = randi(3,[1,N]);     degs.out = randi(3,[1,M]);
depmat = struct();
depmat.in = logical(randi([0,1],d,N));
depmat.out = logical(randi([0,1],d,M));

Aop = rand_tensopvar(dims,vars,dom,depmat,nterms,type,degs,dnsty);

v = cell(1,d);
for j=1:d
    pvars_j = polynomial(Aop.ops{1,j}.vars.in(:));
    n_j = size(Aop.ops{1,j},2);
    v{j} = rand_poly([n_j,p],pvars_j,pdegv,50);
end

Av = apply_op(Aop,v);

% -Aop applied should equal -(Aop applied)
Bop = -Aop;
Bv = apply_op(Bop,v);
err1 = cleanpoly(Bv - (-Av),1e-10);
if max(max(abs(err1.C)))>tol
    error("Trial %d (type=[%d,%d]): uminus: residual exceeds tolerance.",trial,type(1),type(2))
end

% Aop + (-Aop) applied should equal 0
Cop = Aop + (-Aop);
Cv = apply_op(Cop,v);
err2 = cleanpoly(Cv,1e-10);
if max(max(abs(err2.C)))>tol
    error("Trial %d (type=[%d,%d]): Aop+(-Aop): residual exceeds tolerance (should cancel to 0).",trial,type(1),type(2))
end
end

fprintf('Test_tensopvar_uminus passed (%d trials, otype(1) in {1,2} alternating).\n',ntrials)
