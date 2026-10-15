% Test_tensopvar_times.m -- tests @tensopvar/times.m (the '.*' operator)
% for a scalar multiplier. The oracle is apply_op(a.*B,v) against
% a*apply_op(B,v) on random inputs v. otype(1) alternates over {1,2},
% otype(2)=2.
%
% DJ (CC), 10/15/2026: Initial coding;

% NOTES
% - Elementwise multiplication by a non-scalar matrix is not tested:
%   times.m does 'A.*C.ops{i,1}' with a 'sopvar', which has no times
%   method (the builtin rejects the object).

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

Bop = rand_tensopvar(dims,vars,dom,depmat,nterms,type,degs,dnsty);
m = size(Bop,1);   n = size(Bop,2);

v = cell(1,d);
for j=1:d
    pvars_j = polynomial(Bop.ops{1,j}.vars.in(:));
    n_j = size(Bop.ops{1,j},2);
    v{j} = rand_poly([n_j,p],pvars_j,pdegv,50);
end
Bv = apply_op(Bop,v);

%%%% Case 1: scalar multiplier
a = randn;
C1 = a.*Bop;
C1v = apply_op(C1,v);
err1 = cleanpoly(C1v - a*Bv,1e-10);
if max(max(abs(err1.C)))>tol
    error("Trial %d (type=[%d,%d]): Case 1 (a.*Bop): residual exceeds tolerance.",trial,type(1),type(2))
end

end

fprintf('Test_tensopvar_times passed (%d trials, otype(1) in {1,2} alternating).\n',ntrials)
